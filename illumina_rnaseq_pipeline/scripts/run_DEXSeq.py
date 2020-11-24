import re
import os
from glob import glob

import numpy as np
import pandas as pd

from rpy2 import robjects as robj
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

import click

base = importr('base')
dexseq = importr('DEXSeq')
biocparallel = importr('BiocParallel')


def as_r_dataframe(pd_df):
    with localconverter(robj.default_converter + pandas2ri.converter):
        r_from_pd_df = robj.conversion.py2rpy(pd_df)
    return r_from_pd_df


def as_r_matrix(pd_df):
    return base.as_matrix(as_r_dataframe(pd_df))


def as_pd_dataframe(r_df):
    with localconverter(robj.default_converter + pandas2ri.converter):
        pd_from_r_df = robj.conversion.rpy2py(base.as_data_frame(r_df))
    return pd_from_r_df


def get_sample_ids(conds, reps):
    sample_names = np.repeat(conds, reps).tolist()
    sample_nos = np.concatenate([np.arange(n) + 1 for n in reps])
    sample_ids = [f'{sn}_{n}' for sn, n in zip(sample_names, sample_nos)]
    return sample_names, sample_ids


def read_featurecounts_output(counts_fn):
    counts = pd.read_csv(
        counts_fn, sep='\t', comment='#',
        dtype={'Chr': str}
    )
    gene_id = counts.Geneid
    gene_id.name = 'gene_id'
    exon_id = counts.groupby('Geneid').cumcount()
    exon_id = exon_id
    exon_id = gene_id + '.' + exon_id.astype(str)
    invs = counts.iloc[:, :6]
    invs.index = exon_id
    invs.columns = ['gene_id', 'chrom', 'start', 'end', 'strand', 'length']
    invs = invs.drop(['gene_id', 'length'], axis=1)
    counts = counts.iloc[:, 6:]

    # hardcoded for this snakemake pipeline,
    # needs some work to make it more generalisable
    counts.columns = counts.columns.str.extract('^aligned_data/(.+).sorted.bam$')[0]
    counts.index = exon_id
    return gene_id.tolist(), exon_id.tolist(), invs, counts


def get_cond_data(counts, cntrl_cond, treat_cond):
    sample_names = counts.columns.to_series()
    sample_names.name = None
    conds = sample_names.str.split('_', expand=True)[0]
    conds.name = 'condition'
    assert all(conds.isin([treat_cond, cntrl_cond]))
    return pd.DataFrame(conds, index=sample_names)


def build_dexseq_dataset(gene_id, exon_id, counts, cntrl_cond, treat_cond):
    cond_data = get_cond_data(counts, cntrl_cond, treat_cond)
    r_cond_data = as_r_dataframe(cond_data)
    r_gene_id = robj.StrVector(gene_id)
    r_exon_id = robj.StrVector(exon_id)
    r_counts = as_r_matrix(counts)
    design = robj.Formula('~sample + exon + condition:exon')
    dxd = dexseq.DEXSeqDataSet(
        countData=r_counts,
        sampleData=r_cond_data,
        design=design,
        featureID=r_exon_id,
        groupID=r_gene_id,
    )
    return dxd


def convert_dexseq_results(dxd):

    def get_lfc_col(dxd):
        cols = list(dxd.slots['listData'].names)
        c = [c for c in cols if c.startswith('log2fold')]
        assert len(c) == 1
        return c[0]

    # convert to a pandas object and save
    dxd_listdata = dxd.slots['listData']
    logfc_col = get_lfc_col(dxd)
    dxd_res = pd.DataFrame.from_dict({
        'gene_id': list(dxd_listdata.rx2('groupID')),
        'exon_id': list(dxd_listdata.rx2('featureID')),
        'CPM': list(dxd_listdata.rx2('exonBaseMean')),
        logfc_col: list(dxd_listdata.rx2(logfc_col)),
        'p_val': list(dxd_listdata.rx2('pvalue')),
        'fdr': list(dxd_listdata.rx2('padj'))
    })
    dxd_res = dxd_res.set_index('exon_id')
    dxd_res.index.name = None
    return dxd_res


def get_gene_level_fdr(dxr):
    per_gene_fdr = dexseq.perGeneQValue(dxr)
    gene_ids = list(per_gene_fdr.names)
    per_gene_fdr = np.array(per_gene_fdr)
    return pd.Series(per_gene_fdr, index=gene_ids)


def test_differential_exon_usage(gene_id, exon_id, invs, counts,
                                 cntrl_cond, treat_cond, processes=1):
    dxd = build_dexseq_dataset(
        gene_id, exon_id, counts, cntrl_cond, treat_cond
    )
    if processes > 1:
        mp_param = biocparallel.MulticoreParam(workers=processes)
    else:
        mp_param = biocparallel.SerialParam()
    dxd = dexseq.estimateSizeFactors_DEXSeqDataSet(dxd)
    dxd = dexseq.estimateDispersions_DEXSeqDataSet(dxd, BPPARAM=mp_param)
    dxd = dexseq.testForDEU(dxd, BPPARAM=mp_param)
    dxd = dexseq.estimateExonFoldChanges(
        dxd,
        fitExpToVar="condition",
        denominator=cntrl_cond,
    )
    dxr = dexseq.DEXSeqResults(dxd)
    exon_level_results = convert_dexseq_results(dxr)
    exon_level_results = pd.concat([invs, exon_level_results], axis=1)
    per_gene_fdr = get_gene_level_fdr(dxr)
    exon_level_results['gene_level_fdr'] = exon_level_results.gene_id.map(per_gene_fdr)
    return exon_level_results


@click.command()
@click.option('-f', '--featurecounts-fn', required=True)
@click.option('-o', '--output-fn', required=True)
@click.option('-c', '--cntrl-cond', required=True)
@click.option('-t', '--treat-cond', required=True)
@click.option('-p', '--processes', default=1)
def main(featurecounts_fn, output_fn, cntrl_cond, treat_cond, processes):
    gene_id, exon_id, invs, counts = read_featurecounts_output(
        featurecounts_fn
    )
    res = test_differential_exon_usage(
        gene_id, exon_id, invs, counts,
        cntrl_cond, treat_cond, processes
    )
    with open(output_fn, 'w') as bed:
        for record in res.itertuples(index=False):
            (
                chrom, start, end, strand, gene_id,
                cpm, logfc, p, fdr, g_fdr
            ) = record

            bed.write(
                f'{chrom}\t{start - 1:d}\t{end:d}\t{gene_id}\t.\t{strand}\t'
                f'{cpm:.2f}\t{logfc:.4f}\t{p:.4g}\t{fdr:.4g}\t{g_fdr:.4g}\n'
            )


if __name__ == '__main__':
    main()