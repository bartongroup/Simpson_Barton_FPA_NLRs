rule get_introns:
    input:
        gtf=config['gtf_fn']
    output:
        'annotation/annot_introns.bed'
    conda:
        'env_yamls/minimap2.yaml'
    shell:
        '''
        paftools.js gff2bed -j <(sed -e 's/gene_id/gene_name/' {input.gtf}) > {output}
        '''


rule map_with_minimap2:
    input:
        fastq=basecalling('basecalled_data/{sample_name}.dna.fastq'),
        introns='annotation/annot_introns.bed',
        reference=config['fasta_fn'],
    output:
        'aligned_data/{sample_name,[^./]+}.bam',
    threads: 12
    params:
        intron_size=config['minimap2_parameters'].get('max_intron_size', 20_000),
        junc_bonus=config['minimap2_parameters'].get('annot_intron_bonus', 12)
    conda:
        'env_yamls/minimap2.yaml'
    shell:
        '''
        ../../scripts/minimap2/minimap2 -t {threads} \
          -a -L --cs=short \
          -k15 -w5 --splice \
          -g2000 -G{params.intron_size} \
          -A1 -B2 -O2,32 -E1,0 -C9 \
          -z200 -uf --splice-flank=yes \
          --junc-bonus={params.junc_bonus} --junc-bed {input.introns} \
          {input.reference} {input.fastq} |
        samtools view -bS - |
        samtools sort -m 1G -@ {threads} -o - - > {output}
        samtools index {output}
        '''


def sample_name_subset(cond):
    sample_names = glob_wildcards(
        basecalling('basecalled_data/{sample_name}.dna.fastq.gz')
    ).sample_name
    cond_sample_names = [sn for sn in sample_names if sn.startswith(cond)]
    return cond_sample_names


def pool_input(wildcards):
    return expand(
        'aligned_data/{sample_name}.filtered.bam',
        sample_name=sample_name_subset(wildcards.cond)
    )


rule pool_bams:
    input:
        pool_input
    output:
        'aligned_data/pooled/{cond}.bam'
    threads: 12
    conda:
        'env_yamls/minimap2.yaml'
    shell:
        '''
        samtools merge -r -@ {threads} {output} {input}
        samtools index {output}
        '''