# epiAneufinder: Identifying copy number variations from single-cell ATAC-seq data

epiAneufinder is an algorithm used for calling Copy Number Variations (CNVs) from single-cell ATAC (scATAC) data.
Single-cell open chromatin profiling via the single-cell Assay for Transposase-Accessible Chromatin using sequencing (scATAC-seq) assay has become a mainstream measurement of open chromatin in single-cells. epiAneufinder exploits the read count information from scATAC-seq data to extract genome-wide copy number variations (CNVs) for each individual cell. epiAneufinder allows the addition of single-cell CNV information to scATAC-seq data, without the need of additional experiments, unlocking a layer of genomic variation which is otherwise unexplored. 

### Description

The algorithm works in three steps:
1. Data preparation (binning, GC correction, removal of blacklisted regions)
2. Genome segmentation based on maximum Anderson-Darling distance
3. Gain/loss assignments

### Getting Started

epiAneufinder is an R package

### Installation

Installation example:  
*OPTIONAL: Create a conda enviroment and install R and devtools
```
conda create -n epianeufinder r-base r-essentials
```
1. Start R and install dependencies
```
install.packages(c("devtools", "BiocManager", "ggdendro"))
BiocManager::install(c("GenomicAlignments", "SummarizedExperiment", "plyranges", "Rsamtools", "GenomeInfoDb", "BSgenome.Hsapiens.UCSC.hg38", "GenomicRanges", "Biostrings", "BiocGenerics", "S4Vectors", "GenomicFeatures"))
```
2. Load devtools and install the epiAneufinder package
```
library(devtools)
install_github("colomemaria/epiAneufinder")
```
The user should also download the blacklisted regions of the hg38 genome and update the parameter accordingly.
By default, the genome version to be used is the hg38. If there is another genome/species that needs to be used, it should be installed by the user, along with the corresponding blacklisted regions.

### Executing program

```
library(epiAneufinder)
epiAneufinder(input="sample.tsv", #Enter path to your fragments.tsv file or the folder containing bam files
              outdir="epiAneufinder_results", #Path to the directory where results should be written 
              blacklist="blacklist.bed", #Path to bed file that contains the blacklisted regions of your genome
              windowSize=1e5, 
              genome="BSgenome.Hsapiens.UCSC.hg38", #Substitute with relevant BSgenome
              exclude=c('chrX','chrY','chrM'), 
              reuse.existing=TRUE,
              title_karyo="Karyogram of sample data", 
              ncores=4,
              minFrags=20000)
```

We also provide a wrapper script to be used (epiAneufinder_wrapper.R). 

The user needs to update the path/name of the input data folder, as well as the path/name of the blacklisted regions to use on the wrapper. Test data and the hg38 blacklisted regions can be found in the sample_data folder. 

The user may also change on the wrapper any of the parameters of the algorithm, if needed. 

Once the user has updated the parameters on the wrapper, the wrapper can be called from console with the following command:

```
Rscript epiAneufinder_wrapper.R
```

### Help



### Authors

Contributors names and contact info

Akshaya Ramakrishnan (akshaya4r@gmail.com)

Aikaterini Symeonidi (aikaterini.symeonidi@helmholtz-muenchen.de) 

Patrick Hanel (patrick.hanel@helmholtz-muenchen.de) 

Michael Schubert  

Maria Colomé-Tatché (maria.colome@helmholtz-muenchen.de)

### Version History

* 0.1
    * Initial Release

### License

This project is licensed under the GNU General Public License v3.0 - see the LICENSE.md file for details

### Acknowledgments


