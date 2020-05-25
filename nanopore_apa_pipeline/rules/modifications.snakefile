def sample_name_subset(cond):
    sample_names = glob_wildcards(
        basecalling('raw_data/{sample_name}')
    ).sample_name
    cond_sample_names = [sn for sn in sample_names if sn.startswith(cond)]
    return cond_sample_names


def differr_input(wildcards):
    return {
        'cntrl_bams': expand(
             'aligned_data/{sample_name}.bam',
             sample_name=sample_name_subset(wildcards.cntrl)
        ),
        'treat_bams': expand(
             'aligned_data/{sample_name}.bam',
             sample_name=sample_name_subset(wildcards.treat)
        )
    }


rule der_analysis:
    input:
        unpack(differr_input)
    output:
        bed='modifications/{treat}_vs_{cntrl}.der_sites.bed',
        thresh='modifications/{treat}_vs_{cntrl}.der_sites_thresholded.bed',
    params:
        cntrl_flag = lambda wc, input: ' '.join([f'-a {fn}' for fn in input.cntrl_bams]),
        treat_flag = lambda wc, input: ' '.join([f'-b {fn}' for fn in input.treat_bams]),
        fasta = config['fasta_fn']
    threads: 20
    conda:
        'env_yamls/differr.yaml'
    shell:
        '''
        differr \
          -r {params.fasta} \
          {params.cntrl_flag} \
          {params.treat_flag} \
          -o {output.bed}
        awk '$9 > 5' {output.bed} > {output.thresh}
        '''
