rule fastqc:
    input:
        'raw_data/{sample_name}.{read}.fastq.gz'
    output:
        'qc/{sample_name}.{read}_fastqc.html'
    conda:
        'env_yamls/fastqc.yaml'
    shell:
        '''
        fastqc -o qc {input}
        '''


def multiqc_input(wildcards):

    sample_names = glob_wildcards(
        'raw_data/{sample_name}.1.fastq.gz'
    ).sample_name

    if wildcards.qc_type == 'qc':
        return expand(
            'qc/{sample_name}.{read}_fastqc.html',
            sample_name=sample_names,
            read=[1, 2]
        )
    elif wildcards.qc_type == 'aligned_data':
        return expand(
            'aligned_data/{sample_name}.sorted.bamstats',
            sample_name=sample_names,
        )
    elif wildcards.qc_type == 'quantification':
        return expand(
            'quantification/{sample_name}/quant.sf',
            sample_name=sample_names,
        )
    else:
        raise NotImplementedError()
        


rule multiqc:
    input:
        multiqc_input
    output:
        '{qc_type}/multiqc_report.html'
    conda:
        'env_yamls/multiqc.yaml'
    shell:
        '''
        multiqc -f -dd 2 -o {wildcards.qc_type} {wildcards.qc_type}
        '''