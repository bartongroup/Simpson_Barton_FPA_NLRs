rule run_d3pendr:
    input:
        cntrl_bams = expand(
            'aligned_data/{sample_name}.filtered.bam.bai',
            sample_name=config['control_sample_names']
        ),
        treat_bams = expand(
            'aligned_data/{sample_name}.filtered.bam.bai',
            sample_name=config['treatment_sample_names']
        ),
        gtf=config['gtf_fn'],
    output:
        'apa_results/{comp}.apa_results.bed'
    params:
        cntrl_flag=lambda wc, input: ' '.join([f'-c {fn[:-4]}' for fn in input.cntrl_bams]),
        treat_flag=lambda wc, input: ' '.join([f'-t {fn[:-4]}' for fn in input.treat_bams]),
        output_prefix=lambda wc: f'apa_results/{wc.comp}',
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
          --bootstraps {params.bootstraps} \
          --min-read-overlap {params.min_read_overlap} \
          {params.use_model} {params.test_hom}
        '''
