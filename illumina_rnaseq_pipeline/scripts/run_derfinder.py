import re
import numpy as np
import pandas as pd
from rpy2 import robjects as robj
from rpy2.robjects.packages import importr
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import pandas2ri

import pyBigWig as pybw
import click

derfinder = importr('derfinder')
edgeR = importr('edgeR')
limma = importr('limma')


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


def get_chroms_list(bw_fn):
    with pybw.open(bw_fn) as bw:
        chroms = list(bw.chroms().keys())
    return chroms


def load_data(fns, sample_names, chrom, cutoff=10):
    fns = robj.StrVector(fns)
    fns.names = robj.StrVector(sample_names)
    chroms = robj.StrVector([chrom])
    full_cov = derfinder.fullCoverage(
        files=fns, chrs=chroms,
        cutoff=cutoff, verbose=False,
    )
    return full_cov


def derfinder_granges_to_pd_dataframe(granges):
    invs = to_pd_dataframe(granges)
    invs = invs[['seqnames', 'start', 'end', 'strand']]
    invs.columns = ['chrom', 'start', 'end', 'strand']
    invs['chrom'] = invs.chrom.str.replace('^chr', '')
    invs['start'] = invs['start'] - 1
    return invs
    

def derfinder_counts_to_pd_dataframe(counts):
    counts = to_pd_dataframe(counts)
    return counts.round()


def get_region_counts(fns, sample_names, read_length, cutoff=10):
    chroms = get_chroms_list(fns[0])
    region_counts = []
    for chrom in chroms:
        cov = load_data(
            fns, sample_names, chrom, cutoff=cutoff
        )
        regions = derfinder.regionMatrix(
            cov, cutoff=cutoff, L=read_length, verbose=False
        )
        assert len(regions.names) == 1
        regions = regions.rx2(regions.names[0])
        invs = derfinder_granges_to_pd_dataframe(regions.rx2('regions'))
        counts = derfinder_counts_to_pd_dataframe(regions.rx2('coverageMatrix'))
        cutoff_mask = (counts >= cutoff).any(1)
        region_counts.append(
            pd.concat([invs, counts], axis=1)[cutoff_mask]
        )
    return pd.concat(region_counts, axis=0).reset_index(drop=True)


def get_expressed_regions(fwd_bw_fns, rev_bw_fns, sample_names,
                          read_length, cutoff=10):
    regions = []
    for strand, fns in zip(['+', '-'], [fwd_bw_fns, rev_bw_fns]):
        r = get_region_counts(fns, sample_names, read_length, cutoff)
        r.strand = strand
        regions.append(r)
    regions = pd.concat(regions).sort_values(['chrom', 'start'])
    return regions


def get_design_matrix(conds, sample_names):
    dm = pd.get_dummies(conds)
    dm.index = sample_names
    dm = to_r_matrix(dm)
    return dm


def create_dgelist(counts, design):
    dgelist = edgeR.DGEList(counts)
    dgelist = edgeR.calcNormFactors(dgelist)
    dgelist = edgeR.estimateDisp(dgelist, design=design)
    return dgelist


def get_de(dgelist, design, cntrl_cond, treat_cond):
    fit = edgeR.glmFit(dgelist, design)
    contrast = limma.makeContrasts(f'{treat_cond} - {cntrl_cond}', levels=design)
    cfit = edgeR.glmLRT(fit, contrast=contrast)
    return cfit


def cfit_to_df(dgelist, cfit):
    cpm = edgeR.cpm(dgelist, log=True)
    cpm_df = to_pd_dataframe(cpm)
    tt = edgeR.topTags(cfit, n=np.inf, adjust_method='BH', sort_by='none')
    tt_df = to_pd_dataframe(tt)
    tt_df = tt_df.join(cpm_df, how='outer')
    return tt_df.reset_index(drop=True)


def run_edgeR(expressed_regions,
              sample_names, conds,
              cntrl_cond, treat_cond):
    invs = expressed_regions[['chrom', 'start', 'end', 'strand']]
    counts = expressed_regions[sample_names]

    dm = get_design_matrix(conds, sample_names)
    dgelist = create_dgelist(
        to_r_matrix(counts), dm
    )
    cfit = get_de(dgelist, dm, cntrl_cond, treat_cond)
    res = cfit_to_df(dgelist, cfit)
    return pd.concat([invs, res], axis=1)


def get_sample_names(c, reps):
    conds = np.repeat(c, reps).tolist()
    ns = np.concatenate([np.arange(n) + 1 for n in reps])
    sample_names = [f'{c}_{i}' for sn, n in zip(conds, n)]
    return conds, sample_names


def find_differentially_expressed_regions(cntrl_fwd_bw_fns, cntrl_rev_bw_fns,
                                          treat_fwd_bw_fns, treat_rev_bw_fns,
                                          cntrl_cond, treat_cond,
                                          read_length, cutoff):
    conds, sample_names = get_sample_names(
        [cntrl_cond, treat_cond],
        [len(cntrl_fwd_bw_fns), len(treat_fwd_bw_fns)]
    )
    expressed_regions = get_expressed_regions(
        list(cntrl_fwd_bw_fns) + list(treat_fwd_bw_fns),
        list(cntrl_rev_bw_fns) + list(treat_rev_bw_fns),
        sample_names, read_length, cutoff
    )
    res = run_edgeR(
        expressed_regions,
        sample_names, conds,
        cntrl_cond, treat_cond
    )
    return res


def clean_name_syntax(name):
    if re.match('^\d', name):
        name = 'p' + name
    name = name.replace('-', '_')
    return name


@click.command()
@click.option('-cf', '--cntrl-fwd-bws', multiple=True, required=True)
@click.option('-cr', '--cntrl-rev-bws', multiple=True, required=True)
@click.option('-tf', '--treat-fwd-bws', multiple=True, required=True)
@click.option('-tr', '--treat-rev-bws', multiple=True, required=True)
@click.option('-o', '--output-fn', required=True)
@click.option('-cn', '--cntrl-name', required=True)
@click.option('-tn', '--treat-name', required=True)
@click.option('-l', '--read-length', required=True)
@click.option('-c', '--min-reads-cutoff', required=False, default=10)
def cli(cntrl_fwd_bws, cntrl_rev_bws,
        treat_fwd_bws, treat_rev_bws,
        output_fn, cntrl_name, treat_name,
        read_length, min_reads_cutoff):
    cntrl_name = clean_name_syntax(cntrl_name)
    treat_name = clean_name_syntax(treat_name)
    res = find_differentially_expressed_regions(
        cntrl_fwd_bws, cntrl_rev_bws,
        treat_fwd_bws, treat_rev_bws,
        cntrl_name, treat_name,
        read_length, min_reads_cutoff
    )
    res.to_csv(output_fn, sep='\t')


if __name___ == '__main__':
    cli()