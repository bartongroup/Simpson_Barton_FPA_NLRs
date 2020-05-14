import os
import copy
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from rpy2 import robjects as robj
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

import click

proteus = importr('proteus')
imputeLCMD = importr('imputeLCMD')


def to_r_dataframe(pd_df):
    with localconverter(robj.default_converter + pandas2ri.converter):
        r_from_pd_df = robj.conversion.py2rpy(pd_df)
    return r_from_pd_df


def to_r_matrix(pd_df):
    return robj.r('as.matrix')(to_r_dataframe(pd_df))


def to_pd_dataframe(sobj):
    with localconverter(robj.default_converter + pandas2ri.converter):
        pd_from_r_df = robj.conversion.rpy2py(robj.r('as.data.frame')(sobj))
    return pd_from_r_df


def build_protein_table(evi_fn, metadata):
    evi = proteus.readEvidenceFile(evi_fn)
    pepdat = proteus.makePeptideTable(evi, to_r_dataframe(metadata))
    pepdat_df = to_pd_dataframe(pepdat.rx2('tab'))
    prodat = proteus.makeProteinTable(pepdat)
    prodat = proteus.normalizeData(prodat)
    return prodat


def _mar_mnar_impute(group, k=5):
    mnar_mask = group.isnull().all(axis=1).values
    # uses QRILC for MNAR data and KNN for MAR data
    mar_imputed = 2 ** to_pd_dataframe(
        imputeLCMD.impute_wrapper_KNN(
            to_r_matrix(np.log2(group)), k
        )
    )
    mnar_imputed = 2 ** to_pd_dataframe(
        imputeLCMD.impute_QRILC(
            to_r_matrix(np.log2(group)),
        ).rx2(1)
    )
    mar_imputed[mnar_mask] = mnar_imputed[mnar_mask]
    return mnar_imputed


def mar_mnar_impute(prodat, metadata, k=5):
    prodat = copy.deepcopy(prodat)
    prodat_df = to_pd_dataframe(prodat.rx2('tab'))
    prodat_df = prodat_df.dropna(how='all', axis=0)
    # impute the missing values
    cond_map = metadata.set_index('sample').condition.to_dict()
    prodat_imp = (prodat_df.groupby(cond_map, axis=1)
                           .apply(_mar_mnar_impute, k=k))
    # stick the imputed vals back in the proteusData list and clean up
    prodat[prodat.names.index('tab')] = to_r_matrix(prodat_imp)
    prodat[prodat.names.index('stats')] = proteus.intensityStats(prodat)
    prodat[prodat.names.index('detect')] = proteus.goodData(prodat)
    return prodat


def proteus_de(prodat, cntrl_cond, treat_cond):
    res = proteus.limmaDE(
        prodat,
        conditions=robj.StrVector([cntrl_cond, treat_cond])
    )
    res = to_pd_dataframe(res).set_index('protein')
    res.columns = res.columns.str.replace('.', '_')
    return res


def combine_p_values(bootstrap_res):
    pval_res = []
    logfc_res = []
    av_exprs_res = []
    idx = bootstrap_res[0].index
    for res in bootstrap_res:
        res = res.reindex(idx)
        pval_res.append(res.P_Value)
        logfc_res.append(res.logFC.values)
        av_exprs_res.append(res.AveExpr.values)
    ensemble_res = bootstrap_res[0].index.to_frame()
    ensemble_res['pval'] = stats.hmean(pval_res, axis=0)
    _, ensemble_res['fdr'], *_ = multipletests(ensemble_res.pval, method='fdr_bh')
    (ensemble_res['lower_ci_logFC'],
     ensemble_res['median_logFC'],
     ensemble_res['upper_ci_logFC']) = np.percentile(
        logfc_res, q=(2.5, 50, 97.5),  axis=0)
    (ensemble_res['lower_ci_av_exprs'],
     ensemble_res['median_av_exprs'],
     ensemble_res['upper_ci_av_exprs']) = np.percentile(
        av_exprs_res, q=(2.5, 50, 97.5),  axis=0)
    return ensemble_res.set_index('protein')


def bootstrap_proteus(evi_fn, metadata, cntrl_cond, treat_cond,
                      nboots=100, k=5):
    # first get the protein data table
    prodat = build_protein_table(evi_fn, metadata)
    bootstrap_res = []
    # now bootstrap the imputation and DE analysis
    for _ in range(nboots):
        prodat_imp = mar_mnar_impute(prodat, metadata, k)
        res = proteus_de(prodat_imp, cntrl_cond, treat_cond)
        bootstrap_res.append(res)
    res = combine_p_values(bootstrap_res)
    n_detected = to_pd_dataframe(prodat.rx2('stats'))
    n_detected = n_detected.pivot(
        index='id', columns='condition', values='ngood'
    )
    res = res.join(n_detected, how='left')
    return res


def generate_metadata(cntrl_sn, treat_sn, cntrl, treat):
    cntrl_sn = cntrl_sn.split(',')
    treat_sn = treat_sn.split(',')
    sample_names = cntrl_sn + treat_sn
    conds = [cntrl for _ in cntrl_sn] + [treat for _ in treat_sn]
    reps = list(range(1, len(cntrl_sn) + 1)) + \
           list(range(1, len(treat_sn) + 1))
    metadata = pd.DataFrame.from_dict({
        'experiment': sample_names,
        'measure': ['Intensity' for _ in sample_names],
        'sample': sample_names,
        'condition': conds,
        'replicate': reps
    })
    return metadata


@click.command()
@click.option('-e', '--maxquant-evidence-fn', required=True)
@click.option('-o', '--output-fn', required=True)
@click.option('-c', '--cntrl-sample-names', required=True)
@click.option('-t', '--treat-sample-names', required=True)
@click.option('--cntrl-cond-name', required=False, default='cntrl')
@click.option('--treat-cond-name', required=False, default='treat')
@click.option('-n', '--nboots', required=False, default=999)
@click.option('-k', '--knn', default=5)
def proteus_cli(maxquant_evidence_fn, output_fn,
                cntrl_sample_names, treat_sample_names,
                cntrl_cond_name, treat_cond_name,
                nboots, knn):
    metadata = generate_metadata(
        cntrl_sample_names, treat_sample_names,
        cntrl_cond_name, treat_cond_name
    )
    res = bootstrap_proteus(
        maxquant_evidence_fn, metadata,
        cntrl_cond_name, treat_cond_name,
        nboots, knn
    )
    res.to_csv(output_fn, sep='\t')


if __name__ == '__main__':
    proteus_cli()