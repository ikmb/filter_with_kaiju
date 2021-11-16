![](images/ikmb_bfx_logo.png)

# Filter with Kaiju

This pipeline filter your input sequences based on Kaiju ([Menzel, P. et al. (2016)](http://www.nature.com/ncomms/2016/160413/ncomms11257/full/ncomms11257.html), [Github page](https://github.com/bioinformatics-centre/kaiju)) output. It was developed as part of the [CRC1182](https://www.metaorganism-research.com/) and inspired by the workflow used by [Schmittmann et al 2021](https://www.frontiersin.org/articles/10.3389/fimmu.2021.689051/full).

When working with metaorganisms, is often desired to filter out sequencing reads from one symbiosis entities. This is what this pipeline is designed for. For instance, [Schmittmann et al 2021](https://www.frontiersin.org/articles/10.3389/fimmu.2021.689051/full) used this approach to filter out putative microbial reads from the transcriptome of sponge *Halichondria panicea* before *de novo* assembly of the transcriptome.

This pipeline does not contain novel algorithms. These are based on the following softwares:

- [fastp](https://github.com/OpenGene/fastp)
- [kaiju](https://kaiju.binf.ku.dk/)
- [seqtk](https://github.com/lh3/seqtk)

## Documentation 

1. [What happens in this pipeline?](docs/pipeline.md)
2. [Installation and configuration](docs/installation.md)
3. [Running the pipeline](docs/usage.md)
4. [Output](docs/output.md)
5. [Troubleshooting](docs/troubleshooting.md)



