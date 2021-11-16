# Pipeline structure



This pipeline retrieves fastq.gz files from a **sample file**, quality filter the reads using [fastp](https://academic.oup.com/bioinformatics/article/34/17/i884/5093234), classify them with [Kaiju](https://www.nature.com/articles/ncomms11257), and extract reads from the original fastq.gz files using [seqtk](https://github.com/lh3/seqtk) based on kaiju classification status, *i.e.*, whether reads were classified or unclassified.

A simplified workflow is shown below.

![](images/pipeline.png)

















