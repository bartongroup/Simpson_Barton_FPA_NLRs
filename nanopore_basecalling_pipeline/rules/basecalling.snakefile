import os
from glob import glob

checkpoint guppy_basecall:
    '''
    rule for basecalling FAST5 data with guppy.
    It produces multiple fastq files whose names are hard to
    preempt so we need to use checkpoints. Also it is prone
    to crashing on the cluster so we need to be able to resume
    a previous basecalling run by having a dummy output file.
    '''
    input:
        fast5_root='raw_data/{sample_name}'
    params:
        flowcell=config['flowcell'],
        kit=config['kit'],
    output:
        "basecalling/{sample_name}/{sample_name}.complete"
    threads: 28
    shell:
        '''
        ../scripts/ont-guppy-cpu/bin/guppy_basecaller \
          --recursive \
          --flowcell {params.flowcell} \
          --kit {params.kit} \
          --num_callers {threads} \
          --cpu_threads_per_caller 4 \
          --records_per_fastq 0 \
          --reverse_sequence yes \
          --input_path {input.fast5_root} \
          --save_path basecalling/{wildcards.sample_name}
        touch {output}
        '''


def get_all_basecalled_fastqs(wildcards):
    checkpoint_output = checkpoints.guppy_basecall.get(**wildcards).output[0]
    checkpoint_output = os.path.split(checkpoint_output)[0]
    return glob(
        os.path.join(checkpoint_output, 'fastq_runid_*.fastq')
    )


rule concatenate_fastqs:
    input:
        get_all_basecalled_fastqs
    output:
        "basecalled_data/{sample_name}.rna.fastq"
    threads: 1
    shell:
        '''
        for FASTQ in {input}; 
        do
          cat $FASTQ >> {output}
        done
        '''


rule rna_to_dna:
    input:
        "basecalled_data/{sample_name}.rna.fastq"
    output:
        "basecalled_data/{sample_name,[^.]+}.dna.fastq"
    threads: 1
    conda:
        'env_yamls/seqkit.yaml'
    shell:
        "cat {input} | seqkit seq --rna2dna > {output}"


rule rename_sequencing_summary:
    input:
        "basecalling/{sample_name}/{sample_name}.complete"
    output:
        "sequencing_summaries/{sample_name}_sequencing_summary.txt"
    threads: 1
    shell:
        '''
        cp basecalling/{wildcards.sample_name}/sequencing_summary.txt {output}
        '''