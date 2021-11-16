# Outputs

The output of each step is outputted in each respective folder, named after the software used.

- **fastp**: Fasq files containing reads which passed quality control.
- **kaiju**: Classification status of reads that were kept.
- **MultiQC**: Info about pipeline and software versions.

- **pipeline_info**: Contains general info about the pipeline, resources used and run.

- **seqtk**: Contains the reads extracted from the original input fastq.gz files listed in the `sample.csv`. This is what you are probably looking for. Files are named after *sample* field from the `sample.csv` file; *e.g.* `sample_1.qc.keep.R1.fq.gz`.

  





