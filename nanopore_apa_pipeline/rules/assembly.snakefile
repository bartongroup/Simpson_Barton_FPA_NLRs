def merged_stringtie_input(wildcards):
    sample_names = glob_wildcards(
        basecalling('basecalled_data/{sample_name}.dna.fastq')
    ).sample_name
    assembly_conds = config['assembly_conds']
    sample_names = [
       sn for sn in sample_names
       if sn.rsplit('_', 1)[0] in assembly_conds
    ]
    return {
        'bams': expand(
            'aligned_data/{sample_name}.filtered.bam',
             sample_name=sample_names,
        ),
        'annot_gtf': config['gtf_fn']
    }


rule merged_stringtie_assembly:
    input:
        unpack(merged_stringtie_input)
    output:
        'assembly/merged_nanopore_assembly.gtf'
    threads: 28
    conda:
        'env_yamls/stringtie.yaml'
    shell:
        '''
        stringtie -p {threads} -t -L \
          -o {output}.tmp.gtf \
          {input.bams}
        cat {input.annot_gtf} {output}.tmp.gtf |
          gffread -T -M -K -Q > {output}
        rm {output}.tmp.gtf
        '''