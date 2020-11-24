## Widespread premature transcription termination controls Arabidopsis thaliana NLR genes

This repository contains the pipelines and notebooks used to generate figures for the manuscript "Widespread premature transcription termination controls Arabidopsis thaliana NLR genes".

### Pipelines

* All pipelines should be run with `snakemake --use-conda`.
* The pipeline `nanopore_apa_pipeline` depends on the output of `nanopore_basecalling_pipeline`.
* The pipelines `illumina_rnaseq_pipeline` and `helicos_apa_pipeline` both depend on the output of `nanopore_apa_pipeline`.

### Data availability:

IVI-MS data is available from the ProteomeXchange Consortium via the PRIDE partner repository with the dataset identifier PXD022684. Col-0 nanopore DRS data is available from ENA accession PRJEB32782. fpa-8 and 35S::FPA:YFP nanopore DRS data is available from ENA accession PRJEB41451. Col-0, fpa-8 and 35S::FPA:YFP Helicos DRS data is available from ENA accession XXX. Col-0, fpa-8 and 35S::FPA:YFP Illumina RNA-Seq data is available from ENA accession PRJEB41455.