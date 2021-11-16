# Usage information

A typical run looks like:

```
nextflow run ikmb/filter_with_kaiju \
--samples samples.csv \
--kaiju_db PATH_TO/kaijudb \
--library pe
--keep U
```

Here, reads left **u**nclassified by kaiju database will be extracted from the input fastq.gz paired-end (**pe**) library files listed in the `samples.csv`.

The file `samples.csv` is a unquoted, comma-separated file containing the information about samples' fastq.gz reads to be processed. It can have as many columns as wanted. But first three columns **must** have headers as following: 

sample,fastq_1,fastq_2,...

For instance: 

```
sample,fastq_1,fastq_2,strandedness,sample.id,name,species,replicate
halichondria.panicea_REP1,/PATH_TO/20Nov6_A01-L1_S31_L001_R1_001.fastq.gz,/PATH_TO/20Nov6_A01-L1_S31_L001_R2_001.fastq.gz,reverse,68_T0S1,20Nov6_A01,Halichondria panicea,1
halichondria.panicea_REP2,/PATH_TO/20Nov6_B01-L1_S33_L001_R1_001.fastq.gz,/PATH_TO/20Nov6_B01-L1_S33_L001_R2_001.fastq.gz,reverse,69_T0S2,20Nov6_B01,Halichondria panicea,2
halichondria.panicea_REP3,/PATH_TO/20Nov6_C01-L1_S35_L001_R1_001.fastq.gz,/PATH_TO/20Nov6_C01-L1_S35_L001_R2_001.fastq.gz,reverse,70_T0S3,20Nov6_C01,Halichondria panicea,3
```

If single end reads, fill fastq_2 column with NA. Do not leave it empty.

**Note:** In this pipeline, the final file names match with the sample ids, which should be unique.

## Choosing Kaiju database and what to keep

The choice of Kaiju database and the parameter what to keep is coupled and depends on what to you to get out of it. For instance, in the case of  [Schmittmann et al 2021](https://www.frontiersin.org/articles/10.3389/fimmu.2021.689051/full) the authors aimed to exclude putative microbial reads from the host transcriptome. In this pipeline, that would translate to using a database with microbial sequences, but not animal, such as (nr, see [list](https://github.com/bioinformatics-centre/kaiju)) and keep unclassified (`--keep U`) reads, *i.e.*, reads that were not classified as either Archaea, Bacteria or Virus. Kaiju also allow for using your own customized database. Please refer to their [manual](https://github.com/bioinformatics-centre/kaiju) for more info about their databases.











