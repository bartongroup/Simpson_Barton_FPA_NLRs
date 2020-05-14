import os
import glob
import re
from bs4 import BeautifulSoup


def add_fasta_fn(mqvar_xml, fasta_fn, identifier_rule=None, description_rule=None):
    if identifier_rule is None:
        identifier_rule = '>([^\s]+)'
    if description_rule is None:
        description_rule = ''
    fasta_xml = mqvar_xml.find('fastaFiles')
    fasta_xml.fastaFilePath.string = os.path.abspath(fasta_fn)
    fasta_xml.identifierParseRule.string = identifier_rule
    fasta_xml.descriptionParseRule.string = description_rule
    return mqvar_xml


def add_sample_fns(mqvar_xml, sample_fns):
    filepaths_xml = mqvar_xml.find('filePaths')
    filepaths_xml.clear()
    
    sample_names_xml = mqvar_xml.find('experiments')
    sample_names_xml.clear()
    
    fractions_xml = mqvar_xml.find('fractions')
    fractions_xml.clear()
    
    ptms_xml = mqvar_xml.find('ptms')
    ptms_xml.clear()
    
    group_inds_xml = mqvar_xml.find('paramGroupIndices')
    group_inds_xml.clear()
    
    ref_channels_xml = mqvar_xml.find('referenceChannel')
    ref_channels_xml.clear()

    for sn, fns in sample_fns.items():
        for fn in fns:
            fn_tag = mqvar_xml.new_tag('string')
            fn_tag.string = os.path.abspath(fn)
            filepaths_xml.append(fn_tag)

            sn_tag = mqvar_xml.new_tag('string')
            sn_tag.string = sn
            sample_names_xml.append(sn_tag)

            fractions_tag = mqvar_xml.new_tag('short')
            fractions_tag.string = '32767'
            fractions_xml.append(fractions_tag)

            ptms_tag = mqvar_xml.new_tag('boolean')
            ptms_tag.string = 'False'
            ptms_xml.append(ptms_tag)

            group_inds_tag = mqvar_xml.new_tag('int')
            group_inds_tag.string = '0'
            group_inds_xml.append(group_inds_tag)

            ref_channels_tag = mqvar_xml.new_tag('string')
            ref_channels_tag.string = ''
            ref_channels_xml.append(ref_channels_tag)
    return mqvar_xml


def add_threads(mqvar_xml, threads):
    mqvar_xml.numThreads.string = str(threads)
    return mqvar_xml


def add_output_paths(mqvar_xml, output_xml_fn):
    mqvar_xml.fixedCombinedFolder.string = os.path.join(os.getcwd(), 'maxquant_output')
    mqvar_xml.tempFolder.string = os.path.join(os.getcwd(), 'maxquant_output/tmp')
    return mqvar_xml


# this code is run by snakemake
input_template_xml_fn = os.path.abspath(os.path.join(os.getcwd(), '../../scripts/template.xml'))
print(input_template_xml_fn)
output_xml_fn = snakemake.output.xml
with open(input_template_xml_fn) as template, open(output_xml_fn, 'wb') as out:
    mqvar_xml = BeautifulSoup(template.read().encode('utf-8'), 'xml')
    # bs4 automatically strips out end-tags when the text is empty... and maxQuant can't deal
    for tag in mqvar_xml.find_all():
        if len(tag.get_text(strip=True)) == 0 and len(list(tag.find_all())) == 0:
            # re-setting the string to an empty string fixes this
            tag.string = ''
    mqvar_xml = add_fasta_fn(mqvar_xml, snakemake.config['fasta_fn'])
    mqvar_xml = add_sample_fns(mqvar_xml, snakemake.config['sample_fns'])
    mqvar_xml = add_threads(mqvar_xml, snakemake.config['threads'])
    mqvar_xml = add_output_paths(mqvar_xml, output_xml_fn)
    out.write(mqvar_xml.encode('utf-8', formatter=None))