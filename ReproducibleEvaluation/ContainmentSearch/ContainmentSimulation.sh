#!/bin/bash
set -e


source ~/.bashrc
conda activate OPAL
#################################################
#Human Genome
Name='Human'
wgsim -N 100000 -1 20000 -2 20000 -e 0.001 -r 0.001 -R 0.001 -X 0.0001 /Data/GCF_000001405.40_GRCh38.p14_genomic.fna /FASTQ-files/${Name}-HiFi1.fastq /FASTQ-files/${Name}-HiFi2.fastq


#################################################
#Largest genome
Name='Pinus_taeda'
wgsim -N 100000 -1 20000 -2 20000 -e 0.001 -r 0.001 -R 0.001 -X 0.0001 /Data/Pinus_taeda_loblolly-pine/GCA_001447015.2_Sugar_pine_JHU_assembly_genomic.fna /FASTQ-files/${Name}-HiFi1.fastq /FASTQ-files/${Name}-HiFi2.fastq


#################################################
#RefSeq-100GB
Name='RefSeq-100GB'
wgsim -N 100000 -1 20000 -2 20000 -e 0.001 -r 0.001 -R 0.001 -X 0.0001 /Data/all-organisms.fna /FASTQ-files/${Name}-HiFi1.fastq /FASTQ-files/${Name}-HiFi2.fastq


##################################################
#RefSeq-200GB
Name='RefSeq-200GB'
wgsim -N 10000 -1 20000 -2 20000 -e 0.001 -r 0.001 -R 0.001 -X 0.0001 /Data/all-organisms-200GB.fna /FASTQ-files/${Name}-HiFi1.fastq /FASTQ-files/${Name}-HiFi2.fastq

