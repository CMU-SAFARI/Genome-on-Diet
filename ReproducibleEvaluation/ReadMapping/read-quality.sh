#!/bin/bash
#set -e



#conda install -c bioconda nanoplot


NanoPlot --fastq D1_S1_L001_R1_001-017.fastq --threads 40 --outdir Read-quality --prefix Illumina

NanoPlot --fastq m64011_190830_220126.fastq --threads 40 --outdir Read-quality --prefix HiFi

NanoPlot --fastq HG002_ONT-UL_GIAB_20200204_1000filtered_2Mreads.fastq --threads 40 --outdir Read-quality --prefix ONT