#!/bin/bash
#set -e



#source ~/.bashrc
#conda activate SNIFFLES

cd /home/alserm/minimap2-alser/Use-cases/Read-mapping/Read-quality

sbatch --exclusive -w kratos1 -c 40 -J bwa_hg38 --wrap="NanoPlot --fastq /home/alserm/minimap2-alser/Use-cases/Read-mapping/HG002_NA24385_son-GRCh38-17Sept2022/D1_S1_L001_R1_001-017.fastq --threads 40 --outdir /home/alserm/minimap2-alser/Use-cases/Read-mapping/Read-quality --prefix Illumina"

sbatch --exclusive -w kratos2 -c 40 -J bwa_hg38 --wrap="NanoPlot --fastq /home/alserm/minimap2-alser/Use-cases/Read-mapping/HG002_NA24385_son-GRCh38-17Sept2022/m64011_190830_220126.fastq --threads 40 --outdir /home/alserm/minimap2-alser/Use-cases/Read-mapping/Read-quality --prefix HiFi"

sbatch --exclusive -w kratos3 -c 40 -J bwa_hg38 --wrap="NanoPlot --fastq /home/alserm/minimap2-alser/Use-cases/Read-mapping/HG002_NA24385_son-GRCh38-17Sept2022/HG002_ONT-UL_GIAB_20200204_1000filtered_2Mreads.fastq --threads 40 --outdir /home/alserm/minimap2-alser/Use-cases/Read-mapping/Read-quality --prefix ONT"