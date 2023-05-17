#!/bin/bash

#SBATCH -o BAFExtract_SNU601_%J.txt
#SBATCH -o BAFExtract_SNU601_%J.err
#SBATCH -J BAFExtract_SNU601
#SBATCH -c 1
#SBATCH --mem=32G
#SBATCH -t 48:00:00
#SBATCH -p slim18

#Load the required moduls
module load ngs/samtools/1.9

#Create directory with BAFExtract bins if the dir doesn't exist yet (not don't by the tool itself)
mkdir -p SNU601_bins_BAFExtract

#Run step 1
samtools view data/SNU601/SNU601_mapped/outs/possorted_genome_bam.bam | tools/BAFExtract/bin/BAFExtract -generate_compressed_pileup_per_SAM stdin genome_input_BAFExtract/hg19_genome_list SNU601_bins_BAFExtract 50 0
echo „Finished step 1 - generating compressed pileup per SAM“

#Run step 2
tools/BAFExtract/bin/BAFExtract -get_SNVs_per_pileup genome_input_BAFExtract/hg19_genome_list SNU601_bins_BAFExtract genome_input_BAFExtract/hg19 20 4 0.1 SNU601_BAFExtract.af
echo „Finished step 2 - getting SNVs per pileup“ 


