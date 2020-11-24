rule fastqc:
    input:
        'raw_data/{sample_name}.fastq.gz'
    output:
        'qc/{sample_name}_fastqc.html'
    conda:
        'env_yamls/fastqc.yaml'
    shell:
        '''
        fastqc -o qc {input}
        '''


def multiqc_input(wildcards):

    sample_names = glob_wildcards(
        'raw_data/{sample_name}.fastq.gz'
    ).sample_name

    if wildcards.qc_type == 'qc':
        return expand(
            'qc/{sample_name}_fastqc.html',
            sample_name=sample_names,
        )
    elif wildcards.qc_type == 'aligned_data':
        sample_names = set([sn.rsplit('.', 1)[0] for sn in sample_names])
        return expand(
            'aligned_data/{sample_name}.sorted.bamstats',
            sample_name=sample_names,
        )
    elif wildcards.qc_type == 'quantification':
        sample_names = set([sn.rsplit('.', 1)[0] for sn in sample_names])
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