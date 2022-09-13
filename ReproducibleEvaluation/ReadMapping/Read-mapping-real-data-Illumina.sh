#!/bin/bash
#set -e

cd /home/alserm/minimap2-alser/src2
#source ~/.bashrc
#conda activate Metalign

#FASTQpath='/home/alserm/Molecules2Variations/GM24149_1_2Mreads.fastq'
#FASTQpath='/home/alserm/Molecules2Variations/PBmixSequel729_1_A01_PBTH_30hours_19kbV2PD_70pM_HumanHG003.fastq'
FASTQpath='/home/alserm/Molecules2Variations/ERR240727_1.fastq'



for i in {9..21..1}
do


    sbatch --exclusive -w kratos5 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Genome-on-Diet-Chr1-Illumina-stats-command-out_k21w$i.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Genome-on-Diet-Chr1-Illumina-stats-time-memory_k21w$i.txt" /home/alserm/minimap2-alser/src2/minimap2_avx -t 40 -ax sr -Z 10 -W 2 -i 2 -k 21 -w $i -N 1 -r 0.05,120,200 -n 0.9,0 --AF_max_loc 10 --secondary=yes -a -o "/home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Illumina/Genome-on-Diet-Chr1-Illumina-stats_k21w$i.sam" /home/alserm/minimap2-alser/Data/GRCh38.p14_NC_000001.11.fna ${FASTQpath}"

    #sbatch --exclusive -w kratos5 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/minimap2-Chr1-Illumina-stats-command-out_k21w$i.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/minimap2-Chr1-Illumina-stats-time-memory_k21w$i.txt" /home/alserm/minimap2-alser/minimap2-original-7June2022/minimap2 -t 40 -ax sr -N 1 -a -k 21 -w $i --secondary=yes -o "/home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Illumina/minimap2-Chr1-Illumina-stats_k21w$i.sam" /home/alserm/minimap2-alser/Data/GRCh38.p14_NC_000001.11.fna ${FASTQpath}"

    sbatch --exclusive -w kratos5 -c 40 -J bwa_hg38 --wrap="samtools stats --threads 40 "/home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Illumina/Genome-on-Diet-Chr1-Illumina-stats_k21w$i.sam" | grep ^SN | cut -f 2- > "/home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Genome-on-Diet-Chr1-Illumina-stats_k21w$i.txt""

  
    #sbatch --exclusive -w kratos5 -c 40 -J bwa_hg38 --wrap="samtools stats --threads 40 "/home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Illumina/minimap2-Chr1-Illumina-stats_k21w$i.sam" | grep ^SN | cut -f 2- > "/home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/minimap2-Chr1-Illumina-stats_k21w$i.txt""

done

#samtools flagstat ../Use-cases/Read-mapping/Mapping-accuracy/ONT/Genome-on-Diet-Chr1-ONT-pbsim2fq.sam | awk -F "[(|%]" 'NR == 5 {print $0}'


#ssh alserm@safari-proxy
#wget https://sh.rustup.rs -O rustup.sh && sh rustup.sh
#source $HOME/.cargo/env
#cargo install samcomp
#
#samcomp -m prim /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/ONT/minimap2-Chr1-ONT-stats_k21w10.sam /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/ONT/Genome-on-Diet-Chr1-ONT-stats_k21w10.sam 


samcomp -m prim /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Illumina/minimap2-Chr1-Illumina-stats_k21w11.sam /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Illumina/Genome-on-Diet-Chr1-Illumina-stats_k21w11.sam 


grep 'seconds' /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/minimap2-Chr1-Illumina*
grep 'mapped:' /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/minimap2-Chr1-Illumina* | awk 'ORS=NR%3?FS:RS'
grep 'non-primary alignments' /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/minimap2-Chr1-Illumina*
grep 'Maximum resident set size' /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/minimap2-Chr1-Illumina*


grep 'seconds' /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Genome-on-Diet-Chr1-Illumina*
grep 'mapped:' /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Genome-on-Diet-Chr1-Illumina* | awk 'ORS=NR%3?FS:RS'
grep 'non-primary alignments' /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Genome-on-Diet-Chr1-Illumina*
grep 'Maximum resident set size' /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Genome-on-Diet-Chr1-Illumina*


