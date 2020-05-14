rule samtools_sort:
    input:
        bam='aligned_data/{sample_name}.bam',
    output:
        bam='aligned_data/{sample_name}.sorted.bam',
    threads: 8
    conda:
        'env_yamls/samtools.yaml'
    shell:
        '''
        samtools sort -@ {threads} -m 1G -o - {input.bam} > {output.bam}
        '''


rule samtools_index:
    input:
        bam='aligned_data/{sample_name}.{bamtype}.bam',
    output:
        bam='aligned_data/{sample_name}.{bamtype}.bam.bai'
    conda:
        'env_yamls/samtools.yaml'
    shell:
        '''
        samtools index {input.bam}
        '''


rule dustmask:
    input:
        bam='aligned_data/{sample_name}.sorted.bam',
        bai='aligned_data/{sample_name}.sorted.bam.bai'
    output:
        bam='aligned_data/{sample_name}.dustmasked.bam',
    params:
        dust=config['dustmask_fn']
    conda:
        'env_yamls/bedtools.yaml'
    shell:
        '''
        bedtools intersect -v -f 0.2 \
          -a {input.bam} \
          -b {params.dust} > {output.bam}
        '''


rule filter_indels:
    input:
        bam='aligned_data/{sample_name}.dustmasked.bam',
        bai='aligned_data/{sample_name}.dustmasked.bam.bai'
    output:
        bam='aligned_data/{sample_name}.filtered.bam',
    params:
        m=config.get('max_allowed_indels', 4)
    conda:
        'env_yamls/samtools.yaml'
    shell:
        '''
        samtools view -h {input.bam} |
        awk -v OFS='\t' \
          '{{if (substr($1, 1, 1) == "@") {{print}} \
             else if ((gsub(/I/, "I", $6) + gsub(/D/, "D", $6)) <= {params.m}) {{print}}}}' |
        samtools view -bS - > {output.bam}
        samtools index {output.bam}
        '''
