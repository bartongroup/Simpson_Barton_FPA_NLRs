def sample_name_subset(cond):
    sample_names = glob_wildcards('aligned_data/{sample_name,[^./]+}.bam').sample_name
    cond_sample_names = [sn for sn in sample_names if sn.startswith(cond)]
    return cond_sample_names


def expand_with_sample_name_subset(file_pattern):
    def _expand(cond):
        return expand(
            file_pattern,
            sample_name=sample_name_subset(cond)
        )
    return _expand


def d3pendr_input(wildcards):
    bam_expand = expand_with_sample_name_subset(
        'aligned_data/{sample_name}.filtered.bam'
    )
    bai_expand = expand_with_sample_name_subset(
        'aligned_data/{sample_name}.filtered.bam.bai'
    )
    return {
        'cntrl_bams': bam_expand(wildcards.cntrl),
        'cntrl_bais': bai_expand(wildcards.cntrl),
        'treat_bams': bam_expand(wildcards.treat),
        'treat_bais': bai_expand(wildcards.treat),
        'gtf': nanopore('assembly/merged_nanopore_assembly.gtf')
    }


rule run_d3pendr:
    input:
        unpack(d3pendr_input)
    output:
        'apa_results/{treat}_vs_{cntrl}.apa_results.bed'
    params:
        cntrl_flag=lambda wc, input: ' '.join([f'-c {fn}' for fn in input.cntrl_bams]),
        treat_flag=lambda wc, input: ' '.join([f'-t {fn}' for fn in input.treat_bams]),
        output_prefix=lambda wc: f'apa_results/{wc.treat}_vs_{wc.cntrl}',
        bootstraps=config['d3pendr_parameters'].get('nboots', 999),
        min_read_overlap=config['d3pendr_parameters'].get('min_read_overlap', 0.2),
        extend_3p=config['d3pendr_parameters'].get('extend_three_prime', 200),
        use_model='--use-gamma-model' if config['d3pendr_parameters'].get('use_gamma_model', True) \
                                      else '--no-model',
        test_hom='--test-homogeneity' if config['d3pendr_parameters'].get('test_homogeneity', False) \
                                      else '--no-test-homogeneity'
    threads: 24
    conda:
        'env_yamls/d3pendr.yml'
    shell:
        '''
        d3pendr \
          {params.cntrl_flag} \
          {params.treat_flag} \
          -a {input.gtf} \
          -o {params.output_prefix} \
          -p {threads} \
          --read-strand opposite \
          --read-end 5 \
          --use-locus-tag \
          --max-terminal-intron-size 10000 \
          --bootstraps {params.bootstraps} \
          --min-read-overlap {params.min_read_overlap} \
          {params.use_model} {params.test_hom}
        '''
