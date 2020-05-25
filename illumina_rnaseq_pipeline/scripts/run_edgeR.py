import re
import os
from glob import glob

import numpy as np
import pandas as pd

import click

from rpy2 import robjects as robj
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter


tximport = importr("tximport")
edgeR = importr('edgeR')
limma = importr('limma')


def to_r_dataframe(pd_df):
    with localconverter(robj.default_converter + pandas2ri.converter):
        r_from_pd_df = robj.conversion.py2rpy(pd_df)
    return r_from_pd_df


def to_pd_dataframe(sobj):
    with localconverter(robj.default_converter + pandas2ri.converter):
        pd_from_r_df = robj.conversion.rpy2py(robj.r('as.data.frame')(sobj))
    return pd_from_r_df


def get_sample_ids(conds, reps):
    sample_names = np.repeat(conds, reps).tolist()
    sample_nos = np.concatenate([np.arange(n) + 1 for n in reps])
    sample_ids = [f'{sn}_{n}' for sn, n in zip(sample_names, sample_nos)]
    return sample_names, sample_ids


def get_design_matrix(conds, sample_ids):
    dm = pd.get_dummies(conds, dtype='int64')
    dm.index = sample_ids
    dm = to_r_dataframe(dm)
    dm = robj.r['as.matrix'](dm)
    return dm


def read_salmon(fn, sample_ids):
    tx2gene = pd.read_table(fn[0], usecols=[0], names=['TXNAME'], header=0)
    tx2gene['GENENAME'] = tx2gene.TXNAME.str.split('.', expand=True)[0]
    tx2gene = to_r_dataframe(tx2gene)
    fn = robj.StrVector(np.array(fn))
    fn.names = sample_ids
    salmon_data = tximport.tximport(
        fn, type='salmon', tx2gene=tx2gene, dropInfReps=True)
    return salmon_data.rx2('counts'), salmon_data.rx2('length')


def create_dgelist(counts, lengths, design):
    _create_dgelist = robj.r('''
        function (counts, lengths) { 
            normMat <- lengths/exp(rowMeans(log(lengths)))
            o <- log(edgeR::calcNormFactors(counts/normMat)) + log(colSums(counts/normMat))
            y <- edgeR::DGEList(counts)
            y$offset <- t(t(log(normMat)) + o)
            keep <- rowSums(edgeR::cpm(y) > 1) >= 2
            y <- y[keep, keep.lib.sizes=FALSE]
            return (y)
        }
    ''')
    dgelist = _create_dgelist(counts, lengths)
    dgelist = edgeR.estimateDisp(dgelist, design=design)
    return dgelist


def get_de(dgelist, design, cntrl_cond, treat_cond):
    fit = edgeR.glmFit(dgelist, design)
    contrast = limma.makeContrasts(f'{treat_cond} - {cntrl_cond}', levels=design)
    c_fit = edgeR.glmLRT(fit, contrast=contrast)
    return c_fit


def cfit_to_df(dgelist, c_fit):
    cpm = edgeR.cpm(dgelist, log=True)
    cpm_df = to_pd_dataframe(cpm)
    tt = edgeR.topTags(c_fit, n=np.inf, adjust_method='BH', sort_by='none')
    tt_df = to_pd_dataframe(tt)
    tt_df = tt_df.join(cpm_df, how='outer')
    return tt_df


def run_edgeR(cntrl_fns, cntrl_name, treat_fns, treat_name):
    n_cntrl = len(cntrl_fns)
    n_treat = len(treat_fns)
    conds, sample_ids = get_sample_ids([cntrl_name, treat_name], [n_cntrl, n_treat])
    design = get_design_matrix(conds, sample_ids)
    counts, length = read_salmon([*cntrl_fns, *treat_fns], sample_ids)
    dgelist = create_dgelist(counts, length, design)
    c_fit = get_de(dgelist, design, cntrl_name, treat_name)
    return cfit_to_df(dgelist, c_fit)


def clean_name_syntax(name):
    if re.match('^\d', name):
        name = 'c' + name
    name = name.replace('-', '_')
    return name


@click.command()
@click.option('-cf', '--cntrl-quantsf', multiple=True, required=True)
@click.option('-tf', '--treat-quantsf', multiple=True, required=True)
@click.option('-cn', '--cntrl-name', required=True)
@click.option('-tn', '--treat-name', required=True)
@click.option('-o', '--output-fn', required=True)
def cli(cntrl_quantsf, treat_quantsf,
        cntrl_name, treat_name,
        output_fn):
    cntrl_name = clean_name_syntax(cntrl_name)
    treat_name = clean_name_syntax(treat_name)
    cntrl_quantsf = list(cntrl_quantsf)
    treat_quantsf = list(treat_quantsf)
    res = run_edgeR(
        cntrl_quantsf, cntrl_name,
        treat_quantsf, treat_name
    )
    res.to_csv(output_fn, sep='\t')


if __name__ == '__main__':
    cli()