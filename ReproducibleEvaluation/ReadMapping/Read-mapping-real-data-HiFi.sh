#!/bin/bash
#set -e

cd /home/alserm/minimap2-alser/src2
#source ~/.bashrc
#conda activate Metalign

#FASTQpath='/home/alserm/Molecules2Variations/GM24149_1_2Mreads.fastq'
FASTQpath='/home/alserm/Molecules2Variations/PBmixSequel729_1_A01_PBTH_30hours_19kbV2PD_70pM_HumanHG003.fastq'
#FASTQpath='/home/alserm/Molecules2Variations/ERR240727_1.fastq'


for i in {16..25..1}
do
    sbatch --exclusive -w kratos5 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Genome-on-Diet-Chr1-HiFi-stats-command-out_k19w$i.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Genome-on-Diet-Chr1-HiFi-stats-time-memory_k19w$i.txt"  /home/alserm/minimap2-alser/src2/minimap2_avx -t 40 -ax map-hifi -Z 10 -W 2 -i 0.2 -k 19 -w $i -N 1 -r 0.04,400,800 -n 0.8,0.005  --AF_max_loc 10 --sort=merge --frag=no -F200,1 --secondary=yes -a -o "/home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/HiFi/Genome-on-Diet-Chr1-HiFi-stats_k19w$i.sam" /home/alserm/minimap2-alser/Data/GRCh38.p14_NC_000001.11.fna ${FASTQpath}"

    #sbatch --exclusive -w kratos2 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/minimap2-Chr1-HiFi-stats-command-out_k19w$i.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/minimap2-Chr1-HiFi-stats-time-memory_k19w$i.txt" /home/alserm/minimap2-alser/minimap2-original-7June2022/minimap2 -t 40 -ax map-hifi -N 1 -a -k 19 -w $i --secondary=yes -o "/home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/HiFi/minimap2-Chr1-HiFi-stats_k19w$i.sam" /home/alserm/minimap2-alser/Data/GRCh38.p14_NC_000001.11.fna ${FASTQpath}"


    sbatch --exclusive -w kratos5 -c 40 -J bwa_hg38 --wrap="samtools stats --threads 40 "/home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/HiFi/Genome-on-Diet-Chr1-HiFi-stats_k19w$i.sam" | grep ^SN | cut -f 2- > "/home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Genome-on-Diet-Chr1-HiFi-stats_k19w$i.txt""

  
    #sbatch --exclusive -w kratos2 -c 40 -J bwa_hg38 --wrap="samtools stats --threads 40 "/home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/HiFi/minimap2-Chr1-HiFi-stats_k19w$i.sam" | grep ^SN | cut -f 2- > "/home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/minimap2-Chr1-HiFi-stats_k19w$i.txt""

done

#samtools flagstat ../Use-cases/Read-mapping/Mapping-accuracy/ONT/Genome-on-Diet-Chr1-ONT-pbsim2fq.sam | awk -F "[(|%]" 'NR == 5 {print $0}'


#ssh alserm@safari-proxy
#wget https://sh.rustup.rs -O rustup.sh && sh rustup.sh
#source $HOME/.cargo/env
#cargo install samcomp
#
#samcomp -m prim /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/ONT/minimap2-Chr1-ONT-stats_k21w10.sam /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/ONT/Genome-on-Diet-Chr1-ONT-stats_k21w10.sam 




samcomp -m prim /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/HiFi/minimap2-Chr1-HiFi-stats_k19w19.sam /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/HiFi/Genome-on-Diet-Chr1-HiFi-stats_k19w19.sam 




grep 'seconds' /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/minimap2-Chr1-HiFi*
grep 'mapped:' /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/minimap2-Chr1-HiFi* | awk 'ORS=NR%3?FS:RS'
grep 'non-primary alignments' /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/minimap2-Chr1-HiFi*
grep 'Maximum resident set size' /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/minimap2-Chr1-HiFi*


grep 'seconds' /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Genome-on-Diet-Chr1-HiFi*
grep 'mapped:' /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Genome-on-Diet-Chr1-HiFi* | awk 'ORS=NR%3?FS:RS'
grep 'non-primary alignments' /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Genome-on-Diet-Chr1-HiFi*
grep 'Maximum resident set size' /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Genome-on-Diet-Chr1-HiFi*


