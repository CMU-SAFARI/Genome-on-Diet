#!/bin/bash
set -e


source ~/.bashrc
conda activate OPAL
#################################################
#Human Genome
Name='Human'
cd /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/FASTQ-files
sbatch --exclusive -w kratos5 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/FASTQ-files/${Name}-command-out.txt" --wrap="wgsim -N 100000 -1 20000 -2 20000 -e 0.001 -r 0.001 -R 0.001 -X 0.0001 /home/alserm/minimap2-alser/Data/GCF_000001405.40_GRCh38.p14_genomic.fna /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/FASTQ-files/${Name}-HiFi1.fastq /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/FASTQ-files/${Name}-HiFi2.fastq"


#sbatch --exclusive -w kratos5 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/FASTQ-files/${Name}-command-out.txt" --wrap="/home/alserm/minimap2-alser/Use-cases/Tools/PBSIM-PacBio-Simulator/src/pbsim --data-type CCS --depth 2 --prefix ${Name}-HiFi --length-mean 9000 --length-sd 1000  --length-min 8000 --length-max 10000 --accuracy-min 0.9 --model_qc /home/alserm/minimap2-alser/Use-cases/Tools/PBSIM-PacBio-Simulator/data/model_qc_ccs /home/alserm/minimap2-alser/Data/GCF_000001405.40_GRCh38.p14_genomic.fna"
#

#################################################
#Largest genome
Name='Pinus_taeda'
cd /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/FASTQ-files
sbatch --exclusive -w kratos5 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/FASTQ-files/${Name}-command-out.txt" --wrap="wgsim -N 100000 -1 20000 -2 20000 -e 0.001 -r 0.001 -R 0.001 -X 0.0001 /mnt/batty-shared/species/Pinus_taeda_loblolly-pine/GCA_001447015.2_Sugar_pine_JHU_assembly_genomic.fna /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/FASTQ-files/${Name}-HiFi1.fastq /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/FASTQ-files/${Name}-HiFi2.fastq"


#Name='Pinus_taeda'
#cd /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/FASTQ-files
#sbatch --exclusive -w kratos5 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/FASTQ-files/${Name}-command-out.txt" --wrap="/home/alserm/minimap2-alser/Use-cases/Tools/PBSIM-PacBio-Simulator/src/pbsim --data-type CCS --depth 2 --prefix ${Name}-HiFi --length-mean 9000 --length-sd 1000  --length-min 8000 --length-max 10000 --accuracy-min 0.9 --model_qc /home/alserm/minimap2-alser/Use-cases/Tools/PBSIM-PacBio-Simulator/data/model_qc_ccs /mnt/batty-shared/species/Pinus_taeda_loblolly-pine/GCA_001447015.2_Sugar_pine_JHU_assembly_genomic.fna"

#################################################
#RefSeq-100GB
Name='RefSeq-100GB'
cd /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/FASTQ-files
sbatch --exclusive -w kratos4 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/FASTQ-files/${Name}-command-out.txt" --wrap="wgsim -N 100000 -1 20000 -2 20000 -e 0.001 -r 0.001 -R 0.001 -X 0.0001 /home/alserm/minimap2-alser/Data/all-organisms.fna /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/FASTQ-files/${Name}-HiFi1.fastq /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/FASTQ-files/${Name}-HiFi2.fastq"


#Name='RefSeq-100GB'
#cd /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/FASTQ-files
#sbatch --exclusive -w kratos4 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/FASTQ-files/${Name}-command-out.txt" --wrap="/home/alserm/minimap2-alser/Use-cases/Tools/PBSIM-PacBio-Simulator/src/pbsim --data-type CCS --depth 2 --prefix ${Name}-HiFi --length-mean 9000 --length-sd 1000  --length-min 8000 --length-max 10000 --accuracy-min 0.9 --model_qc /home/alserm/minimap2-alser/Use-cases/Tools/PBSIM-PacBio-Simulator/data/model_qc_ccs /home/alserm/minimap2-alser/Data/all-organisms.fna"


##################################################
##RefSeq-200GB
#Name='RefSeq-200GB'
#cd /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/FASTQ-files
#sbatch --exclusive -w kratos3 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/FASTQ-files/${Name}-command-out.txt" --wrap="wgsim -N 10000 -1 20000 -2 20000 -e 0.001 -r 0.001 -R 0.001 -X 0.0001 /home/alserm/minimap2-alser/Data/all-organisms-200GB.fna /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/FASTQ-files/${Name}-HiFi1.fastq /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/FASTQ-files/${Name}-HiFi2.fastq"


#Name='RefSeq-200GB'
#cd /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/FASTQ-files
#sbatch --exclusive -w kratos3 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/FASTQ-files/${Name}-command-out.txt" --wrap="/home/alserm/minimap2-alser/Use-cases/Tools/PBSIM-PacBio-Simulator/src/pbsim --data-type CCS --depth 2 --prefix ${Name}-HiFi --length-mean 9000 --length-sd 1000  --length-min 8000 --length-max 10000 --accuracy-min 0.9 --model_qc /home/alserm/minimap2-alser/Use-cases/Tools/PBSIM-PacBio-Simulator/data/model_qc_ccs /home/alserm/minimap2-alser/Data/all-organisms-200GB.fna"
