# epiAneufinder: Identifying copy number alterations from single-cell ATAC-seq data

epiAneufinder is an algorithm used for calling Copy Number Variations (CNVs) from single-cell ATAC (scATAC) data.
Single-cell open chromatin profiling via the single-cell Assay for Transposase-Accessible Chromatin using sequencing (scATAC-seq) assay has become a mainstream measurement of open chromatin in single-cells. epiAneufinder exploits the read count information from scATAC-seq data to extract genome-wide copy number variations (CNVs) for each individual cell. epiAneufinder allows the addition of single-cell CNV information to scATAC-seq data, without the need of additional experiments, unlocking a layer of genomic variation which is otherwise unexplored. 

Ramakrishnan, A., Symeonidi, A., Hanel, P. et al. epiAneufinder identifies copy number alterations from single-cell ATAC-seq data. Nat Commun 14, 5846 (2023). https://doi.org/10.1038/s41467-023-41076-1

All additional scripts used for the publication can be found here https://github.com/colomemaria/epiAneufinder_analyses.git

**Python version:**  We are currently developing a python version of epiAneufinder, which is still in beta-testing. Feel free to explore it and report any improvement suggestions and issues on the related GitHub: https://github.com/colomemaria/pyEpiAneufinder

### Description

The algorithm works in three steps:
1. Data preparation (binning, GC correction, removal of blacklisted regions)
2. Genome segmentation based on maximum Anderson-Darling distance
3. Gain/loss assignments

**Remark:** We improved the implementation of the multi-threading, so that the runtime compared to the previous versions should improve now considerably!

### Getting Started

epiAneufinder is an R package

### Installation

Installation example:  
*OPTIONAL: Create a conda enviroment and install R and devtools
```bash
conda create -n epianeufinder r-base r-essentials
```
1. Start R and install dependencies
```bash
conda activate epianeufinder; R
```
In R:
```R
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

epiAneufinder has been extensovely tested for the following version of dependencies:
GenomicAlignments v1.28.0, SummarizedExperiment v1.22.0, plyranges v1.12.1, Rsamtools v2.8.0, GenomeInfoDb v1.28.4, BSgenome.Hsapiens.UCSC.hg38 v1.4.3, GenomicRanges v1.44.0, Biostrings v2.60.2, BiocGenerics v0.38.0, S4Vectors v0.30.0, GenomicFeatures v1.44.2, devtools v2.4.3, BiocManager v1.30.16 and ggdendro v0.1.22.


The installation process of epiAneufinder, in a prepared enviroment (all dependencies already installed), takes approximatelly 3 minutes.

The software has been tested for Linux and MacOS. 

### Executing program

This is a short introduction in the epiAneufinder functionality, more information can be found in the vignette [introduction-epiAneufinder](vignettes/introduction-epiAneufinder.html).

```
library(epiAneufinder)
epiAneufinder(input="sample.tsv", #Enter path to your fragments.tsv file or the folder containing bam files
              outdir="epiAneufinder_results", #Path to the directory where results should be written 
              blacklist="hg38-blacklist.v2.bed", #Path to bed file that contains the blacklisted regions of your genome
              windowSize=1e5, 
              genome="BSgenome.Hsapiens.UCSC.hg38", #Substitute with relevant BSgenome
              exclude=c('chrX','chrY','chrM'), 
              reuse.existing=TRUE,
              title_karyo="Karyogram of sample data", 
              ncores=4,
              minFrags=20000,
              minsizeCNV=0,
              k=4,
              plotKaryo=TRUE)
```

We also provide a wrapper script to be used (epiAneufinder_wrapper.R). 

The user needs to update the path/name of the input data folder, as well as the path/name of the blacklisted regions to use on the wrapper. Test data and the hg38 blacklisted regions can be found in the sample_data folder. 

The user may also change on the wrapper any of the parameters of the algorithm, if needed. 

Once the user has updated the parameters on the wrapper, the wrapper can be called from console with the following command:

```
Rscript epiAneufinder_wrapper.R
```

Depending on the available resources, running the demo shouldn't take more than 15 minutes. Tested in a typical desktop, running time was 7.34 minutes with 6Gb of RAM used. 

### Authors

Contributors names and contact info

Akshaya Ramakrishnan (akshaya4r@gmail.com)

Aikaterini Symeonidi (asymeonidi@bmc.med.lmu.de and ksymeonidh@gmail.com) 

Patrick Hanel (patrick.hanel@helmholtz-muenchen.de) 

Katharina Schmid (katharina.schmid@bmc.med.lmu.de)

Maria Richter (maria.richter@bmc.med.lmu.de)

Michael Schubert  

Maria Colomé-Tatché (maria.colome@helmholtz-muenchen.de)

### Version History

* 1.1.5
    * Adapted function to read also fragment files from the new 10X version (with a 6th column)

* 1.1.4
    * Added parameters to `karyotype_measures` to optional create a scatterplot of aneuploidy vs heterogeneity per chromosome
    * Calculate the CNV burden (=aneuploidy) per cell with the function `cnv_burden_per_cell`

* 1.1.3
    * Added the function `karyotype_measures` to calculate aneuploidy and heterogeneity scores of the result karyograms (see vignette).
    
* 1.1.2
    *  Added the boolean parameter `gc_correction` (default TRUE) whether GC correction should be performed. We recommend strongly to run the GC correction (improves performance considerable). This option is implemented for the cases when the input count matrix was already GC corrected before.

* 1.1.1
    *  Improved the `selected_cells` parameter to work also when the named cells were filtered out before.

* 1.1.0
    * Fixed problem with multi-threading that caused runtime problems in the GC correction and AD calculation.
    * Added two new filtering options: 1) providing a file with cell barcodes to give the user more flexibility in selecting cells (parameter for filename `selected_cells`) and 2) removing cells that have less or equal a certain percentage of non-zero windows (parameter `threshold_cells_bins`, default <= 5\% non-zero windows)
    * Switching the order, first filtering lowly covered regions before applying the GC correction.
    * Added the option to save regions that were filtered out in a file removed_regions.tsv (called `save_removed_regions`, default `FALSE`).

* 1.0.3
    * Added the option to start the algorithm directly from a count matrix (Warning: performance might drop compared to bam/fragment files, see vignette)
    * Extended the text output to show a few more infos, e.g. number cells and windows
    * Added the parameter mapqFilter for a flexible filter of the bam files
    * Added a function to split cells into sublcones
    * Added to additional plotting functions, to plot the read distribution of an individual cell and to plot additional annotations in a side bar next to the karyogram
    * Added a vignette to improve the documentation

* 1.0.2
    * Updated function documentations
    * Included a better explaining error message for mismatching chromosome versions
    * Corrected minor bug if one chromosome fits in one window
    * Made karyogram plotting optional (argument plotKaryo)
    * Added in the wrapper the default parameters as used in the paper

* 0.1
    * Initial Release


### License

This project is licensed under the GNU General Public License v3.0 - see the LICENSE.md file for details


