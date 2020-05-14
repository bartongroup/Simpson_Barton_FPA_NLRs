rule filter_spurious_tpe_from_oversplitting:
    input:
        bam='aligned_data/{sample_name}.bam',
        seq_sum=basecalling('sequencing_summaries/{sample_name}_sequencing_summary.txt')
    output:
        bam='aligned_data/{sample_name,[^.]+}.filtered.bam'
    conda:
        'env_yamls/d3pendr.yaml'
    shell:
        '''
        filter_nanopore_oversplitting.py \
          -b {input.bam} \
          -s {input.seq_sum} \
          -o {output}
        '''


rule run_d3pendr:
    input:
        cntrl_bams = expand(
            'aligned_data/{sample_name}.filtered.bam',
            sample_name=config['control_sample_names']
        ),
        treat_bams = expand(
            'aligned_data/{sample_name}.filtered.bam',
            sample_name=config['treatment_sample_names']
        ),
        gtf=config['gtf_fn'],
    output:
        'apa_results/{comp}.apa_results.bed'
    params:
        cntrl_flag=lambda wc, input: ' '.join([f'-c {fn}' for fn in input.cntrl_bams]),
        treat_flag=lambda wc, input: ' '.join([f'-t {fn}' for fn in input.treat_bams]),
        output_prefix=lambda wc: f'apa_results/{wc.comp}',
        bootstraps=config['d3pendr_parameters'].get('nboots', 999),
        min_read_overlap=config['d3pendr_parameters'].get('min_read_overlap', 0.2),
        use_model='--use-gamma-model' if config['d3pendr_parameters'].get('use_gamma_model', True) \
                                      else '--no-model',
        test_hom='--test-homogeneity' if config['d3pendr_parameters'].get('test_homogeneity', False) \
                                      else '--no-test-homogeneity'
    threads: 24
    conda:
        'env_yamls/d3pendr.yaml'
    shell:
        '''
        d3pendr \
          {params.cntrl_flag} \
          {params.treat_flag} \
          -a {input.gtf} \
          -o {params.output_prefix} \
          -p {threads} \
          --bootstraps {params.bootstraps} \
          --min-read-overlap {params.min_read_overlap} \
          {params.use_model} {params.test_hom}
        '''
