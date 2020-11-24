import re
import numpy as np
import pandas as pd

import pysam
import click


def strand_filterer(strand, read_type='paired'):
    if read_type == 'single':
        if strand == '+':
            def _filt(query):
                for r in query:
                    if not r.is_reverse:
                        yield r
        else:
            def _filt(query):
                for r in query:
                    if r.is_reverse:
                        yield r
    elif read_type == 'paired':
        if strand == '+':
            def _filt(query):
                for r in query:
                    if r.is_read1 and r.is_reverse:
                        yield r
                    elif r.is_read2 and not r.is_reverse:
                        yield r
        else:
            def _filt(query):
                for r in query:
                    if r.is_read1 and not r.is_reverse:
                        yield r
                    elif r.is_read2 and r.is_reverse:
                        yield r
    else:
        raise ValueError()
    return _filt


class MultiBamParser:

    def __init__(self, bam_fns, read_type='paired'):
        self.handles = {
            bam_fn: pysam.AlignmentFile(bam_fn) for bam_fn in bam_fns
        }
        self.closed = False
        self.read_type = read_type

    def find_introns(self, chrom, start, end, strand=None, introns=None):
        q = (chrom, start, end)
        if strand is not None:
            filt = strand_filterer(strand, self.read_type)
        results = {}
        for bam_fn, bam in self.handles.items():
            bam_q_iter = bam.fetch(*q)
            if strand is not None:
                bam_q_iter = filt(bam_q_iter)
            found = bam.find_introns(bam_q_iter)
            if introns:
                found = {i: (found[i] if i in found else 0) for i in introns}
            results[bam_fn] = found
        results = pd.DataFrame(results).fillna(0)
        if len(results):
            results.index.names = ['Start', 'End']
            results = results.reset_index()
            results.insert(0, 'Chr', chrom)
            results.insert(3, 'Strand', strand)
            return results
        else:
            return None


    def close(self):
        for bam in self.handles.values():
            bam.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()


def get_record_range(invs):
    starts, ends = zip(*invs)
    return min(starts), max(ends)


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
                 by_locus=True):
    gtf_records = {}
    if by_locus:
        transcript_to_locus_mapping = {}
    with open(gtf_fn) as gtf:
        for i, record in enumerate(gtf):
            record = record.split('\t')
            chrom, _, feat_type, start, end, _, strand = record[:7]
            start = int(start) - 1
            end = int(end)
            if feat_type == 'transcript' and by_locus:
                locus_id = get_gtf_attribute(record, 'locus')
                transcript_id = get_gtf_attribute(record, 'transcript_id')
                transcript_to_locus_mapping[transcript_id] = locus_id
            elif feat_type == 'exon':
                transcript_id = get_gtf_attribute(record, 'transcript_id')
                idx = (chrom, transcript_id, strand)
                if idx not in gtf_records:
                    gtf_records[idx] = {}
                if feat_type not in gtf_records[idx]:
                    gtf_records[idx][feat_type] = []
                gtf_records[idx][feat_type].append((start, end))

    # get introns
    for idx, feat_invs in gtf_records.items():
        exons = feat_invs['exon']
        exons.sort()
        intron_starts = [e[1] for e in exons[:-1]]
        intron_ends = [e[0] for e in exons[1:]]
        gtf_records[idx]['intron'] = list(zip(intron_starts, intron_ends))
    if by_locus:
        # regroup gene invs by locus id:
        gtf_records_by_locus = {}
        for (chrom, transcript_id, strand), feat_invs in gtf_records.items():
            locus_id = transcript_to_locus_mapping[transcript_id]
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
        start, end = get_record_range(feat_invs['exon'])
        introns = set(feat_invs['intron'])
        yield chrom, start, end, gene_id, strand, introns


@click.command()
@click.option('-b', '--bam-fns', multiple=True, required=True)
@click.option('-g', '--gtf-fn', required=True)
@click.option('-o', '--output-counts-fn', required=True)
@click.option('--use-locus-tag/--use-gene-id-tag', default=True)
def main(bam_fns, gtf_fn, output_counts_fn, use_locus_tag=True):

    gtf = gtf_iterator(
        gtf_fn,
        use_locus_tag,
    )
    res = []
    with MultiBamParser(bam_fns) as mbam:

        for chrom, start, end, gene_id, strand, introns in gtf:
            introns = mbam.find_introns(chrom, start, end, strand, introns)
            if introns is not None:
                introns.insert(0, 'Geneid', gene_id)
                res.append(introns)
    results = pd.concat(res)
    results.insert(5, 'Length', 0)
    results.to_csv(output_counts_fn, sep='\t', index=False)
    


if __name__ == '__main__':
    main()