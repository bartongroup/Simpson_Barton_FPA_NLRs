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


def sample_name_subset(cond):
    sample_names = glob_wildcards(
        basecalling('basecalled_data/{sample_name}.dna.fastq')
    ).sample_name
    cond_sample_names = [sn for sn in sample_names if sn.startswith(cond)]
    return cond_sample_names


def d3pendr_input(wildcards):
    return {
        'cntrl_bams': expand(
             'aligned_data/{sample_name}.filtered.bam',
             sample_name=sample_name_subset(wildcards.cntrl)
        ),
        'treat_bams': expand(
             'aligned_data/{sample_name}.filtered.bam',
             sample_name=sample_name_subset(wildcards.treat)
        )
    }


rule run_d3pendr:
    input:
        unpack(d3pendr_input),
        gtf='assembly/merged_nanopore_assembly.gtf',
    output:
        'apa_results/{treat}_vs_{cntrl}.apa_results.bed'
    params:
        cntrl_flag=lambda wc, input: ' '.join([f'-c {fn}' for fn in input.cntrl_bams]),
        treat_flag=lambda wc, input: ' '.join([f'-t {fn}' for fn in input.treat_bams]),
        output_prefix=lambda wc: f'apa_results/{wc.treat}_vs_{wc.cntrl}',
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
          --use-locus-tag \
          --bootstraps {params.bootstraps} \
          --min-read-overlap {params.min_read_overlap} \
          {params.use_model} {params.test_hom}
        '''
