rule build_salmon_index:
    output:
        directory('salmon_index')
    params:
        fasta_fn = config['transcriptome_fasta_fn']
    conda:
        'env_yamls/salmon.yaml'
    threads: 8
    shell:
        '''
        cat {params.fasta_fn} > tmp.fa
        salmon index -p {threads} -i {output} -t tmp.fa
        rm tmp.fa
        '''


rule pseudoalign_with_salmon:
    input:
        read='raw_data/{sample_name}.1.fastq.gz',
        mate='raw_data/{sample_name}.2.fastq.gz',
        index='salmon_index'
    output:
        'quantification/{sample_name}/quant.sf'
    params:
        prefix=lambda wc, output: os.path.split(output[0])[0]
    conda:
        'env_yamls/salmon.yaml'
    threads: 8
    shell:
        '''
        salmon quant -l A -p {threads} \
          -i {input.index} \
          -1 {input.read} -2 {input.mate} \
          -o {params.prefix}
        '''