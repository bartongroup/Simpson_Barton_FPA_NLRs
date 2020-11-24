def sample_name_subset(cond, raw_data_template='raw_data/{sample_name}.1.fastq.gz'):
    sample_names = glob_wildcards(raw_data_template).sample_name
    cond_sample_names = [sn for sn in sample_names if sn.startswith(cond)]
    return cond_sample_names


def cntrl_vs_treat_names(wildcards, file_template, suffix='',
                         raw_data_template='raw_data/{sample_name}.1.fastq.gz'):
    return {
        f'cntrl{suffix}': expand(
            file_template,
            sample_name=sample_name_subset(wildcards.cntrl, raw_data_template)
        ),
        f'treat{suffix}': expand(
            file_template,
            sample_name=sample_name_subset(wildcards.treat, raw_data_template)
        )
    }


def edgeR_input(wildcards):
    return cntrl_vs_treat_names(
        wildcards,
        'quantification/{sample_name}/quant.sf'
    )


rule run_edgeR:
    input:
        unpack(edgeR_input)
    output:
        tsv='differential_expression/edgeR/{treat}_vs_{cntrl}.tsv',
    params:
        cntrl=lambda wc, input: ' '.join(['-cf {}'.format(i) for i in input.cntrl]),
        treat=lambda wc, input: ' '.join(['-tf {}'.format(i) for i in input.treat])
    conda:
        'env_yamls/rpy2_edger.yaml'
    shell:
        '''
        python ../scripts/run_edgeR.py \
          {params.cntrl} \
          {params.treat} \
          -cn {wildcards.cntrl} \
          -tn {wildcards.treat} \
          -o {output.tsv}
        '''


def expressed_regions_input(wildcards):
    fwd = cntrl_vs_treat_names(
        wildcards,
        'coverage_tracks/{sample_name}.fwd.bw',
        suffix='_fwd'
    )
    rev = cntrl_vs_treat_names(
        wildcards,
        'coverage_tracks/{sample_name}.rev.bw',
        suffix='_rev'
    )
    fwd.update(rev)
    return fwd


def get_bw_flags(wildcards, input):
    cntrl = [','.join([i, i.replace('fwd', 'rev')]) for i in input.cntrl_fwd]
    cntrl = ' '.join(['-c {}'.format(i) for i in cntrl])
    treat = [','.join([i, i.replace('fwd', 'rev')]) for i in input.treat_fwd]
    treat = ' '.join(['-t {}'.format(i) for i in treat])
    return {'cntrl': cntrl, 'treat': treat}


rule get_expressed_regions:
    input:
        unpack(expressed_regions_input),
        gtf=nanopore('assembly/merged_nanopore_assembly.gtf')
    output:
        gtf='differential_expression/dexseq/{treat}_vs_{cntrl}.expressed_regions.gtf'
    params:
        bigwigs=get_bw_flags
    conda:
        'env_yamls/expressed_regions.yaml'
    shell:
        '''
        python ../scripts/get_expressed_regions.py \
          {params.bigwigs[cntrl]} \
          {params.bigwigs[treat]} \
          -g {input.gtf} \
          -o {output.gtf}
        '''


rule run_dexseq:
    input:
        'quantification/dexseq/{treat}_vs_{cntrl}.{count_type}.tsv'
    output:
        tsv='differential_expression/dexseq_{count_type}/{treat}_vs_{cntrl}.tsv'
    conda:
        'env_yamls/rpy2_dexseq.yaml'
    threads: 12
    shell:
        '''
        python ../scripts/run_DEXSeq.py -p {threads} \
          -f {input} \
          -c {wildcards.cntrl} \
          -t {wildcards.treat} \
          -o {output.tsv}
        '''