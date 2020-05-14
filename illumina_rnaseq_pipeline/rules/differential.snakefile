def sample_name_subset(cond):
    sample_names = glob_wildcards(
        'raw_data/{sample_name}.1.fastq.gz'
    ).sample_name
    cond_sample_names = [sn for sn in sample_names if sn.startswith(cond)]
    return cond_sample_names


def cntrl_vs_treat_names(wildcards, file_template, suffix=''):
    return {
        f'cntrl{suffix}': expand(
             file_template, sample_name=sample_name_subset(wildcards.cntrl)
        ),
        f'treat{suffix}': expand(
             file_template, sample_name=sample_name_subset(wildcards.treat)
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


def derfinder_input(wildcards):
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


rule run_derfinder:
    input:
        unpack(derfinder_input)
    output:
        tsv='differential_expression/derfinder/{treat}_vs_{cntrl}.tsv'
    params:
        read_length=150,
        cntrl_fwd=lambda wc, input: ' '.join(['-cf {}'.format(i) for i in input.cntrl_fwd]),
        cntrl_rev=lambda wc, input: ' '.join(['-cr {}'.format(i) for i in input.cntrl_rev]),
        treat_fwd=lambda wc, input: ' '.join(['-tf {}'.format(i) for i in input.treat_fwd]),
        treat_rev=lambda wc, input: ' '.join(['-tr {}'.format(i) for i in input.treat_rev]),
    conda:
        'env_yamls/rpy2_derfinder.yaml'
    shell:
        '''
        python ../scripts/run_derfinder.py -l {params.read_length} \
          {params.cntrl_fwd} {params.cntrl_rev} \
          {params.treat_fwd} {params.treat_rev} \
          -cn {wildcards.cntrl} \
          -tn {wildcards.treat} \
          -o {output.tsv}
        '''