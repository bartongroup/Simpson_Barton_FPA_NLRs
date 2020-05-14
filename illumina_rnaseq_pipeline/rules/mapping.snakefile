rule build_STAR_index:
    '''Create the index required for alignment with STAR'''
    output:
        directory('STAR_index')
    threads: 28
    params:
        fasta_fn = config['genome_fasta_fn'],
        gtf_fn = config['gtf_fn'],
        overhang = 149
    conda:
        'env_yamls/star.yaml'
    shell:
        '''
        mkdir {output};
        STAR \
          --runThreadN {threads} \
          --runMode genomeGenerate \
          --genomeDir {output} \
          --genomeFastaFiles {params.fasta_fn} \
          --sjdbGTFfile {params.gtf_fn} \
          --sjdbOverhang {params.overhang}
        '''


rule map_with_STAR:
    '''map reads with STAR spliced aligner'''
    input:
        read='raw_data/{sample_name}.1.fastq.gz',
        mate='raw_data/{sample_name}.2.fastq.gz',
        index='STAR_index'
    output:
        bam='aligned_data/{sample_name}.sorted.bam',
        bai='aligned_data/{sample_name}.sorted.bam.bai',
        bamstats='aligned_data/{sample_name}.sorted.bamstats',
    threads: 28
    conda:
        'env_yamls/star.yaml'
    shell:
        '''
        TOPDIR=$(pwd)
        mkdir -p aligned_data/{wildcards.sample_name}.tmpdir
        cd aligned_data/{wildcards.sample_name}.tmpdir
        STAR \
          --runThreadN {threads} \
          --genomeDir $TOPDIR/{input.index} \
          --readFilesIn $TOPDIR/{input.read} $TOPDIR/{input.mate} \
          --readFilesCommand "zcat" \
          --outFilterMultimapNmax 5 \
          --alignSJoverhangMin 8 \
          --alignSJDBoverhangMin 3 \
          --outFilterMismatchNmax 5 \
          --alignIntronMin 60 \
          --alignIntronMax 20000 \
          --outSAMtype BAM Unsorted
        cd ../..
        samtools sort \
          -m 2G -@ {threads} \
          -o {output.bam} \
          aligned_data/{wildcards.sample_name}.tmpdir/Aligned.out.bam
        samtools index {output.bam}
        samtools flagstat {output.bam} > {output.bamstats}
        rm -rf aligned_data/{wildcards.sample_name}.tmpdir
        '''


rule split_strand:
    input:
        bam='aligned_data/{sample_name}.sorted.bam',
        bai='aligned_data/{sample_name}.sorted.bam.bai'
    output:
        bam='aligned_data/{sample_name}.sorted.{strand}.bam',
        bai='aligned_data/{sample_name}.sorted.{strand}.bam.bai'
    params:
        samflags_1=lambda wc: '-f 128 -F 16' if wc.strand == 'fwd' else '-f 144',
        samflags_2=lambda wc: '-f 80' if wc.strand == 'fwd' else '-f 64 -F 16'
    threads: 4
    conda:
        'env_yamls/star.yaml'
    shell:
        '''
        samtools view -@ {threads} -b {params.samflags_1} {input.bam} > {output.bam}.1.bam
        samtools index -@ {threads} {output.bam}.1.bam
        samtools view -@ {threads} -b {params.samflags_2} {input.bam} > {output.bam}.2.bam
        samtools index -@ {threads} {output.bam}.2.bam
        samtools merge -@ {threads} {output.bam} {output.bam}.1.bam {output.bam}.2.bam
        samtools index -@ {threads} {output.bam}
        rm {output.bam}.[12].bam
        rm {output.bam}.[12].bam.bai
        '''


rule genome_coverage:
    input:
        bam='aligned_data/{sample_name}.sorted.{strand}.bam',
        bai='aligned_data/{sample_name}.sorted.{strand}.bam.bai'
    output:
        'coverage_tracks/{sample_name}.{strand}.bw',
    params:
        genome=config['genome_fasta_fn'],
    conda:
        'env_yamls/star.yaml'
    shell:
        '''
        samtools depth -d0 {input.bam} | 
          awk -v OFS='\t' '{{print $1, $2-1, $2, $3}}' > {output}.tmp.bdg
        bedGraphToBigWig {output}.tmp.bdg <(cut -f-2 {params.genome}.fai) {output}
        rm {output}.tmp.bdg
        '''