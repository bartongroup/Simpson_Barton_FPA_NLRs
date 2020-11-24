## Nanopore APA pipeline

This pipeline aligns, filters and assembles transcripts from the nanopore DRS data, as well as performing differential poly(A) site usage analysis using [d3pendr](https://github.com/bartongroup/d3pendr). It should be run before running the pipelines `helicos_apa_pipeline` and `illumina_rnaseq_pipeline`.