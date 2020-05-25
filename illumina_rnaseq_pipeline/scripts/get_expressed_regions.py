import re
import numpy as np
from scipy import signal, ndimage as ndi
import click

import pyBigWig as pybw


class MultiBigWigParser:

    def __init__(self, bw_fns):
        # attempt to infer if stranded, each bw_fn should be comma separated list
        bw_fns = [tuple(fn.split(',')) for fn in bw_fns]
        if all([len(fn) == 2 for fn in bw_fns]):
            stranded = True
        elif all([len(fn) == 1 for fn in bw_fns]):
            stranded = False
        else:
            raise ValueError('Please provide either single bw files or comma separated pos,neg bw files')
        if stranded:
            self.handles = {
                bw_fn: (pybw.open(bw_fn[0]), pybw.open(bw_fn[1])) for bw_fn in bw_fns
            }
        else:
            self.handles = {
                bw_fn: (pybw.open(bw_fn[0]),) for bw_fn in bw_fns
            }
        self.closed = False
        self.stranded = stranded

    def fetch(self, chrom, start, end, strand=None):
        if strand is not None and not self.stranded:
            raise ValueError('cannot specify strand on unstranded bigwigs')
        if self.stranded:
            if strand == '+':
                queries = [
                    pos_bw.values(chrom, start, end, numpy=True)
                    for pos_bw, neg_bw in self.handles.values()
                ]
            elif strand == '-':
                queries = [
                    neg_bw.values(chrom, start, end, numpy=True)
                    for pos_bw, neg_bw in self.handles.values()
                ]
            elif strand is None or strand == '.':
                queries = [
                    np.nansum([
                        pos_bw.values(chrom, start, end, numpy=True),
                        neg_bw.values(chrom, start, end, numpy=True)
                    ], axis=0)
                    for pos_bw, neg_bw in self.handles.values()
                ]
            else:
                raise ValueError(f'strand {strand} unrecognised')
        else:
            queries = [
                bw[0].values(chrom, start, end, numpy=True)
                for bw in self.handles
            ]
        queries = np.asarray(queries)
        queries[np.isnan(queries)] = 0
        return queries

    def close(self):
        for bws in self.handles.values():
            for bw in bws:
                bw.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()


def flatten_intervals(invs):
    flattened = []
    all_invs = iter(np.sort(invs, axis=0))
    inv_start, inv_end = next(all_invs)
    for start, end in all_invs:
        if start <= inv_end:
            inv_end = max(inv_end, end)
        else:
            flattened.append([inv_start, inv_end])
            inv_start, inv_end = start, end
    if not flattened or flattened[-1] != [inv_start, inv_end]:
        flattened.append([inv_start, inv_end])
    return np.array(flattened)


def filter_terminal_exons(invs, max_intron_size, min_exon_size):
    if len(invs) == 1:
        return invs
    else:
        l_ex = invs[0, 1] - invs[0, 0]
        l_in = invs[1, 0] - invs[0, 1]
        if (l_ex < min_exon_size) or (l_in >= max_intron_size):
            invs = invs[1:]
            if len(invs) == 1:
                return invs
        else:
            r_ex = invs[-1, 1] - invs[-1, 0]
            r_in = invs[-1, 0] - invs[-2, 1]
            if (r_ex < min_exon_size) or (r_in >= max_intron_size):
                invs = invs[:-1]
    return invs


def get_record_range(invs,
                     max_terminal_intron_size=None,
                     min_terminal_exon_size=None,
                     filter_=True):
    invs = flatten_intervals(invs)
    if filter_:
        invs = filter_terminal_exons(
            invs,
            max_terminal_intron_size,
            min_terminal_exon_size,
        )
    return invs[0, 0], invs[-1, 1]


def get_gtf_attribute(gtf_record, attribute):
    try:
        attr = re.search(f'{attribute} "(.+?)";', gtf_record[8]).group(1)
    except AttributeError:
        raise ValueError(
            f'Could not parse attribute {attribute} '
            f'from GTF with feature type {record[2]}'
        )
    return attr


def gtf_iterator(gtf_fn,
                 extend_gene_five_prime=0,
                 use_5utr=False,
                 extend_gene_three_prime=0,
                 by_locus=True,
                 max_terminal_intron_size=100_000,
                 min_terminal_exon_size=20):
    gtf_records = {}
    if by_locus:
        gene_to_locus_mapping = {}
    with open(gtf_fn) as gtf:
        for i, record in enumerate(gtf):
            record = record.split('\t')
            chrom, _, feat_type, start, end, _, strand = record[:7]
            start = int(start) - 1
            end = int(end)
            if feat_type == 'transcript' and by_locus:
                locus_id = get_gtf_attribute(record, 'locus')
                gene_id = get_gtf_attribute(record, 'gene_id')
                gene_to_locus_mapping[gene_id] = locus_id
            elif feat_type == 'CDS' or feat_type == 'exon':
                gene_id = get_gtf_attribute(record, 'gene_id')
                idx = (chrom, gene_id, strand)
                if idx not in gtf_records:
                    gtf_records[idx] = {}
                if feat_type not in gtf_records[idx]:
                    gtf_records[idx][feat_type] = []
                gtf_records[idx][feat_type].append((start, end))

    if by_locus:
        # regroup gene invs by locus id:
        gtf_records_by_locus = {}
        for (chrom, gene_id, strand), feat_invs in gtf_records.items():
            locus_id = gene_to_locus_mapping[gene_id]
            new_idx = (chrom, locus_id, strand)
            if new_idx not in gtf_records_by_locus:
                gtf_records_by_locus[new_idx] = {}
            for feat_type, invs in feat_invs.items():
                if feat_type not in gtf_records_by_locus[new_idx]:
                    gtf_records_by_locus[new_idx][feat_type] = []
                gtf_records_by_locus[new_idx][feat_type] += invs
        gtf_records = gtf_records_by_locus

    # once whole file is parsed yield the intervals
    for (chrom, gene_id, strand), feat_invs in gtf_records.items():
        exon_start, exon_end = get_record_range(
            feat_invs['exon'],
            max_terminal_intron_size,
            min_terminal_exon_size,
        )
        try:
            cds_start, cds_end = get_record_range(
                feat_invs['CDS'],
                filter_=False,
            )
        except KeyError:
            # non-coding RNA
            cds_start, cds_end = exon_start, exon_end

        # remove region corresponding to 5'UTR if necessary
        if use_5utr:
            gene_start = exon_start
            gene_end = exon_end
        else:
            gene_start = cds_start if strand == '+' else exon_start
            gene_end = exon_end if strand == '+' else cds_end
        

        # add extensions to 3' and 5' ends
        start_ext, end_ext = extend_gene_five_prime, extend_gene_three_prime
        if strand == '-':
            start_ext, end_ext = end_ext, start_ext
        gene_start = max(0, gene_start - start_ext)
        gene_end = gene_end + end_ext

        yield chrom, gene_start, gene_end, gene_id, strand


def norm_factors(profiles):
    p = np.sum(profiles, axis=1)
    n = p.mean()
    return n / p


def normalise_profiles(cntrl_arrs, treat_arrs):
    prof = np.concatenate([cntrl_arrs, treat_arrs])
    n = norm_factors(prof)
    prof = (prof.T * n).T
    prof[np.isnan(prof)] = 0
    return prof[:len(cntrl_arrs)], prof[len(cntrl_arrs):]


def get_expressed_regions(cntrl_profiles, treat_profiles,
                          threshold=1, min_size=50):
    cntrl_exprs = (cntrl_profiles >= threshold).all(0)
    treat_exprs = (treat_profiles >= threshold).all(0)
    exprs = (cntrl_exprs | treat_exprs).astype('int32')

    if not exprs.any():
        return np.empty(shape=(0, 2), dtype='int32')
    chgp = np.diff(exprs)

    er_starts = np.argwhere(chgp == 1).ravel() + 1
    er_ends = np.argwhere(chgp  == -1).ravel() + 1
    if er_starts.size == 0:
        er_starts = np.array([0,])
    if er_ends.size == 0: 
        er_ends = np.array([len(exprs)])

    if er_starts[0] > er_ends[0]:
        er_starts = np.insert(er_starts, 0, 0)
    if er_starts[-1] > er_ends[-1]:
        er_ends = np.insert(er_ends, len(er_ends), len(exprs))

    ers = []
    for start, end in zip(er_starts, er_ends):
        if (end - start) >= min_size:
            ers.append((start, end))

    return np.asarray(ers)


def segment_expressed_region(cntrl_profiles, treat_profiles,
                             medfilt_size=101, exprs_change=1,
                             window_size=20, min_segment_size=50):
    # can't segment if region is less than 2* the min seg length
    if cntrl_profiles.shape[1] < (min_segment_size * 2):
        return np.array([])
    prof_logfc = (np.mean(np.log2(treat_profiles + 0.5), axis=0) - 
                  np.mean(np.log2(cntrl_profiles + 0.5), axis=0))
    filter_ = np.ones(window_size * 2)
    filter_[window_size:] *= -1
    prof_change = np.convolve(prof_logfc, filter_, mode='same') / window_size
    peaks, _ = signal.find_peaks(
        np.abs(prof_change),
        height=exprs_change,
        distance=min_segment_size
    )
    # remove first and last peaks if they are too close to edges
    if peaks.size and peaks[0] < min_segment_size:
        peaks = peaks[1:]
    r_size = cntrl_profiles.shape[1]
    if peaks.size and (r_size - peaks[-1]) < min_segment_size:
        peaks = peaks[:-1]
    return peaks


def median_filter_second_axis(profiles, size):
    return np.asarray([
        ndi.median_filter(p, size=size, mode='nearest')
        for p in profiles
    ])


def find_expressed_segments(cntrl_profiles, treat_profiles,
                            expression_threshold=1,
                            cov_medfilt_size=101,
                            rel_exprs_change_for_breakpoint=1,
                            exprs_change_window_size=25,
                            min_segment_size=50):

    cntrl_profiles = median_filter_second_axis(
        cntrl_profiles, size=cov_medfilt_size
    )
    treat_profiles = median_filter_second_axis(
        treat_profiles, size=cov_medfilt_size
    )

    expressed_regions = get_expressed_regions(
        cntrl_profiles, treat_profiles,
        threshold=expression_threshold,
        min_size=min_segment_size,
    )
    if not expressed_regions.size:
        return expressed_regions
    expressed_segments = []
    for start, end in expressed_regions:
        seg_idx = segment_expressed_region(
            cntrl_profiles[:, start: end],
            treat_profiles[:, start: end],
            medfilt_size=cov_medfilt_size,
            exprs_change=rel_exprs_change_for_breakpoint,
            window_size=exprs_change_window_size,
            min_segment_size=min_segment_size,
        )
        segs = np.repeat(seg_idx + start, 2)
        segs = np.insert(segs, [0, len(segs)], [start, end])
        segs = segs.reshape(-1, 2)
        expressed_segments.append(segs)
    return np.concatenate(expressed_segments).astype('int32')


def get_segment_counts(expressed_segments,
                       cntrl_profiles,
                       treat_profiles):
    cntrl_counts = []
    treat_counts = []
    for start, end in expressed_segments:
        c = np.median(cntrl_profiles[:, start: end], axis=1)
        t = np.median(treat_profiles[:, start: end], axis=1)
        cntrl_counts.append(c)
        treat_counts.append(t)
    return np.asarray(cntrl_counts).T, np.asarray(treat_counts).T


@click.command()
@click.option('-c', '--cntrl-bw-fns', multiple=True, required=True)
@click.option('-t', '--treat-bw-fns', multiple=True, required=True)
@click.option('-g', '--gtf-fn', required=True)
@click.option('-o', '--output-gtf-fn', required=True)
@click.option('-t', '--treat-bw-fns', multiple=True, required=True)
@click.option('--extend-gene-five-prime', default=0)
@click.option('--use-5utr/--ignore-5utr', default=True)
@click.option('--extend-gene-three-prime', default=0)
@click.option('--use-locus-tag/--use-gene-id-tag', default=True)
@click.option('--max-terminal-intron-size', default=100_000)
@click.option('--min-terminal-exon-size', default=30)
@click.option('--expression-threshold', default=1)
@click.option('--cov-medfilt-size', default=101)
@click.option('--logfc-for-breakpoint', default=1)
@click.option('--exprs-change-window-size', default=25)
@click.option('--min-segment-size', default=50)
def main(cntrl_bw_fns, treat_bw_fns,
         gtf_fn, output_gtf_fn,
         extend_gene_five_prime=0, use_5utr=True,
         extend_gene_three_prime=0, use_locus_tag=True,
         max_terminal_intron_size=10_000, min_terminal_exon_size=30,
         expression_threshold=1, cov_medfilt_size=101,
         logfc_for_breakpoint=1, exprs_change_window_size=25,
         min_segment_size=50):

    gtf = gtf_iterator(
        gtf_fn,
        extend_gene_five_prime,
        use_5utr,
        extend_gene_three_prime,
        use_locus_tag,
        max_terminal_intron_size,
        min_terminal_exon_size,
    )
    with MultiBigWigParser(cntrl_bw_fns) as cntrl, \
         MultiBigWigParser(treat_bw_fns) as treat, \
         open(output_gtf_fn, 'w') as bed:

        for chrom, start, end, gene_id, strand in gtf:
            cntrl_profiles, treat_profiles = normalise_profiles(
                cntrl.fetch(chrom, start, end, strand),
                treat.fetch(chrom, start, end, strand),
            )
            expressed_segments = find_expressed_segments(
                cntrl_profiles, treat_profiles,
                expression_threshold,
                cov_medfilt_size,
                logfc_for_breakpoint,
                exprs_change_window_size,
                min_segment_size
            )
            if len(expressed_segments) > 1:
                expressed_segments += start
                idx = np.arange(len(expressed_segments)) + 1
                if strand == '-':
                    idx = idx[::-1]
                for i, (start, end) in zip(idx, expressed_segments):
                    bed.write(f'{chrom}\tger\texpressed_region\t'
                              f'{start + 1}\t{end}\t.\t{strand}\t.\t'
                              f'gene_id "{gene_id}"; exon_id "{gene_id}.{i}";\n')


if __name__ == '__main__':
    main()