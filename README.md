# epiAneufinder: identifying copy number variations from single-cell ATAC-seq data

epiAneufinder is an algorithm used for calling Copy Number Variations (CNVs) from single-cell ATAC (scATAC) data.
ingle-cell open chromatin profiling via the single-cell Assay for Transposase-Accessible Chromatin using sequencing (scATAC-seq) assay has become a mainstream measurement of open chromatin in single-cells.
Here we present a novel algorithm, epiAneufinder, which exploits the read count information from scATAC-seq data to extract genome-wide copy number variations (CNVs) for individual cells.
Using epiAneufinder it is possible to identify single cell copy number variations from every scATAC dataset, without the need of additional experiments.


## Description

The algorithm works in three steps:
1. Data preparation (binning, GC correction, removal of blacklisted regions)
2. Genome segmentation based on maximum Anderson-Darling distance
3. Gain/loss assignments

## Getting Started

epiAneufinder is an R package

### Dependencies

* devtools
* BiocManager
* GenomicAlignments 
* SummarizedExperiment 
* plyranges
* Rsamtools
* GenomeInfoDb
* BSgenome.Hsapiens.UCSC.hg38
* GenomicRanges
* Biostrings
* BiocGenerics
* S4Vectors
* GenomicFeatures
* ggdendro

### Installing

The easiest way to install epiAneufinder is through conda.
Installation example:
Optionally: Create a conda enviroment and install R and devtools
```
conda create -n epianeufinder r-base r-essentials r-devtools
```
1. Start R and install dependencies
```
install.packages(c("BiocManager", "ggdendro"))
```
```
BiocManager::install(c("GenomicAlignments", "SummarizedExperiment", "plyranges", "Rsamtools", "GenomeInfoDb", "BSgenome.Hsapiens.UCSC.hg38", "GenomicRanges", "Biostrings", "BiocGenerics", "S4Vectors", "GenomicFeatures"))
```
2. Load devtools and install the epiAneufinder package
```
library(devtools)
install_github("colomemaria/epiAneufinder")
```
The user should also download the blacklisted regions of the hg38 and update accordingly the parameter in the wrapper script.
By default, the genome version to be used is the hg38. If there is another genome, it should be installed by the user, along with the corresponding blacklisted regions.

### Executing program

We provide a wrapper script to be used.
Once the user has updated the parameters that wants to change, the wrapper can be called from console with the following command:

```
Rscript epiAneufinder_wrapper.R
```

## Help

Any advise for common problems or issues.

## Authors

Contributors names and contact info

Akshaya Ramakrishnan  
Aikaterini Symeonidi  
Patrick Hanel  
Michael Schubert  
Maria Colomé-Tatché

## Version History

* 0.1
    * Initial Release

## License

This project is licensed under the [NAME HERE] License - see the LICENSE.md file for details

## Acknowledgments


