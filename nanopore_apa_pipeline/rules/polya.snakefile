rule nanopolish_polya:
    input:
        fastq=basecalling('basecalled_data/{sample_name}.dna.fastq'),
        readdb=basecalling("basecalled_data/{sample_name}.dna.fastq.index.readdb"),
        bam='aligned_data/{sample_name,[^./]+}.bam',
        reference=config['fasta_fn'],
    output:
        polya='polya_tails/{sample_name}.polya_tail_lengths.tsv.gz'
    threads: 28
    conda:
        'env_yamls/nanopolish.yaml'
    resources:
        job_class='long'
    params:
        fastq_data=basecalling('basecalled_data/'),
        fast5_data=basecalling('raw_data/')
    shell:
        '''
        OUTPUT=$(readlink -f {output.polya})
        cp -L --parents `find {params.fastq_data} -name '{wildcards.sample_name}.dna.fastq*'` $TMPDIR
        cp -L --parents `find aligned_data/ -name '{wildcards.sample_name}.bam*'` $TMPDIR
        cp -L --parents `find {params.fast5_data}/{wildcards.sample_name}/ -name '*.fast5'` $TMPDIR
        cd $TMPDIR
        nanopolish polya -t {threads} \
          -r {input.fastq} \
          -b {input.bam} \
          -g {input.reference} |
        gzip > {output}
        '''