# moimix: an R package for evaluating multiplicity of infection in malaria parasites

## Features

* Estimate multiplicity of infection from massively parallel sequencing data
* Estimate heterzygosity and within-isolate diversity directly from read-counts
* Call major alleles within isolates from B-allele frequencies
* Prepare SNP barcode data for use in [COIL](http://www.broadinstitute.org/infect/malaria/coil/)
* Simulate single nucleotide variant data with known multiplicity of infection

## How do I install moimix?

There are plans to put _moimix_ on [Bioconductor](http://www.bioconductor.org/)
in the future, however it is currently only available to install as a development
version from Github:

```{r}
# install using devtools packages
devtools::install_github("sa-lee/moimix")
```

## What data input does _moimix_ require?

_moimix_ makes use of the Genomic Data Storage (GDS) format used
by the Bioconductor package [SeqArray](http://www.bioconductor.org/packages/release/bioc/html/SeqArray.html)
to provide fast access to VCF files in R.

To convert a VCF file to the GDS:
```{r}
library(SeqArray)
seqVCF2GDS("isolate_snps.vcf.gz", "isolate_snps.gds")
```

It is also possible to estimate MOI from a matrix of read counts where
the first column is the number of reads supporting the reference allele and the
second column is the number of reads supporting the alternate allele.

## How do I use _moimix_?
See the [introduction vignette](https://github.com/sa-lee/moimix/blob/master/vignettes/introduction.Rmd) for usage examples.

## How do I cite _moimix_?
Manuscript is currently in preparation. If you use moimix please cite the following

Lee S, Harrison A, Tessier N, Tavul L, Miotto O, Siba P, Kwiatkowski D, Müller I, Barry AE and Bahlo M, _Assessing clonality in malaria parasites using massively parallel sequencing data_, 2016, in preparation.







