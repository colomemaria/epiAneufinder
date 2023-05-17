#!/bin/bash

#SBATCH -o BAFExtract_HCT_%J.txt
#SBATCH -o BAFExtract_HCT_%J.err
#SBATCH -J BAFExtract_HCT
#SBATCH -c 1
#SBATCH --mem=64G
#SBATCH -t 48:00:00
#SBATCH -p slim18

#Load the required moduls
module load ngs/samtools/1.9

#Run step 1
samtools view data/HCT116/HCT116_merged_sorted.bam | tools/BAFExtract/bin/BAFExtract -generate_compressed_pileup_per_SAM stdin data/annotations/genome_input_BAFExtract_hg38/hg38_genome_list HCT116 50 0
echo „Finished step 1 - generating compressed pileup per SAM“

#Run step 2
tools/BAFExtract/bin/BAFExtract -get_SNVs_per_pileup data/annotations/genome_input_BAFExtract_hg38/hg38_genome_list HCT116 data/annotations/genome_input_BAFExtract_hg38/hg38 20 4 0.1 HCT116.af
echo „Finished step 2 - getting SNVs per pileup“ 