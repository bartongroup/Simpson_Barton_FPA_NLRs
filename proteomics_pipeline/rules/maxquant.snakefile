rule create_template_xml:
    output:
        xml='mqvar.xml'
    conda:
        'env_yamls/beautifulsoup.yaml'
    script:
        '../scripts/create_maxquant_xml.py'


rule run_maxquant:
    input:
        xml='mqvar.xml'
    output:
        evi='maxquant_output/combined/txt/evidence.txt'
    threads: 28
    conda:
        'env_yamls/maxquant.yaml'
    shell:
        '''
        mkdir -p maxquant_output
        maxquant -n {input.xml} >> maxquant_output/maxquant_dryrun.log
        maxquant {input.xml}
        '''


rule run_proteus:
    input:
        'maxquant_output/combined/txt/evidence.txt'
    output:
        'proteus_output/{treat}_vs_{cntrl}.tsv'
    params:
        cntrl_sn=lambda wc: ','.join([sn for sn in config['sample_fns'] if sn.startswith(wc.cntrl)]),
        treat_sn=lambda wc: ','.join([sn for sn in config['sample_fns'] if sn.startswith(wc.treat)]),
        nboots=config['nboots']
    threads: 1
    conda:
        'env_yamls/proteus.yaml'
    shell:
        '''
        python ../../scripts/run_proteus.py \
          -e {input} \
          -o {output} \
          -c {params.cntrl_sn} \
          -t {params.treat_sn} \
          --cntrl-cond-name {wildcards.cntrl} \
          --treat-cond-name {wildcards.treat} \
          -n {params.nboots}
        '''