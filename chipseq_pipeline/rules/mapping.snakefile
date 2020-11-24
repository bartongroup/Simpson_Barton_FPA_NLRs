import os


def bowtie_input(wc):
    unpaired = f'raw_data/{wc.sample_name}.fastq.gz'
    if os.path.exists(unpaired):
        return {'read': unpaired}
    else:
        read = f'raw_data/{wc.sample_name}.1.fastq.gz'
        mate = f'raw_data/{wc.sample_name}.2.fastq.gz'
        assert os.path.exists(read) and os.path.exists(mate)
        return {'read': read, 'mate': mate}


def bowtie_input_flags(wc, input):
    if hasattr(input, 'mate'):
        return f'-1 {input.read} -2 {input.mate}'
    else:
        return f'-U {input.read}'


rule map_bowtie2:
    input:
        unpack(bowtie_input)
    output:
        bam='aligned_data/{sample_name}.bam',
        bai='aligned_data/{sample_name}.bam.bai'
    conda:
        'env_yamls/bowtie2.yaml'
    params:
        bt2_prefix=os.path.splitext(config['genome_fasta_fn'])[0],
        input_flags=bowtie_input_flags
    threads: 16
    shell:
        '''
        bowtie2 --threads {threads} --mm --very-sensitive \
          --maxins 800 --no-mixed --no-discordant \
          -x {params.bt2_prefix} {params.input_flags} |
        samtools view -bS - |
        samtools sort -o - - > {output.bam}
        samtools index {output.bam}
        '''


rule get_normalised_coverage:
    input:
        bam='aligned_data/{sample_name}.bam',
    output:
        bw='coverage_tracks/{sample_name}.bw'
    conda:
        'env_yamls/deeptools.yaml'
    threads: 12
    shell:
        '''
        bamCoverage --normalizeUsing CPM -p {threads} \
          --minMappingQuality 5 \
          --binSize 1 \
          -b {input.bam} \
          -o {output.bw}
        '''
