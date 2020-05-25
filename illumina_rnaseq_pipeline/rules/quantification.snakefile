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


def sample_name_subset(cond):
    sample_names = glob_wildcards(
        'raw_data/{sample_name}.1.fastq.gz'
    ).sample_name
    cond_sample_names = [sn for sn in sample_names if sn.startswith(cond)]
    return cond_sample_names


def cntrl_vs_treat_names(wildcards, file_template, suffix=''):
    return {
        f'cntrl{suffix}': expand(
             file_template, sample_name=sample_name_subset(wildcards.cntrl)
        ),
        f'treat{suffix}': expand(
             file_template, sample_name=sample_name_subset(wildcards.treat)
        )
    }


def featurecounts_input(wildcards):
    input_ = cntrl_vs_treat_names(
        wildcards,
        'aligned_data/{sample_name}.sorted.bam',
        suffix='_bams'
    )
    input_['gtf'] = f'differential_expression/dexseq/{wildcards.treat}_vs_{wildcards.cntrl}.expressed_regions.gtf'
    return input_


rule quantify_expressed_regions:
    input:
        unpack(featurecounts_input)
    output:
        'quantification/dexseq/{treat}_vs_{cntrl}.counts.tsv'
    threads: 20
    conda:
        'env_yamls/featurecounts.yaml'
    shell:
        '''
        featureCounts -T {threads} \
          -f -s 2 -O --primary -p -B -C \
          -a {input.gtf} \
          -o {output} \
          -t "expressed_region" \
          -g "gene_id" \
          {input.cntrl_bams} {input.treat_bams}
        '''