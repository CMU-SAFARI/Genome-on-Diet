#!/bin/bash
#set -e

#source ~/.bashrc
#conda activate Metalign

FASTQpath='m64011_190830_220126.fastq'


for i in {7..19..2}
do
    sbatch --exclusive -w kratos3 -c 40 -J bwa_hg38 -o "HiFi/Genome-on-Diet-GRCh38-HiFi-stats-command-out_INDEL_SNP_k19w$i.txt" --wrap="/usr/bin/time -v --output "HiFi/Genome-on-Diet-GRCh38-HiFi-stats-time-memory_INDEL_SNP_k19w$i.txt" ../../Genome-on-Diet-SNPs-Indels/GDiet_avx -t 40 --MD -ax map-hifi -Z 10 -W 2 -i 0.2 -k 19 -w $i -N 1 -r 0.04,400,800 -n 0.8,0.005 --AF_max_loc 10 --sort=merge --frag=no -F200,1 --secondary=yes -a -o "HiFi/Genome-on-Diet-GRCh38-HiFi-stats_INDEL_SNP_k19w$i.sam" GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta ${FASTQpath}"

    sbatch --exclusive -w kratos3 -c 40 -J bwa_hg38 --wrap="samtools stats --threads 40 "HiFi/Genome-on-Diet-GRCh38-HiFi-stats_INDEL_SNP_k19w$i.sam" | grep ^SN | cut -f 2- > "HiFi/Genome-on-Diet-GRCh38-HiFi-stats_INDEL_SNP_k19w$i.txt""

    sbatch --exclusive -w kratos5 -c 40 -J bwa_hg38 -o "HiFi/Genome-on-Diet-GRCh38-HiFi-stats-command-out_SV_k19w$i.txt" --wrap="/usr/bin/time -v --output "HiFi/Genome-on-Diet-GRCh38-HiFi-stats-time-memory_SV_k19w$i.txt" ../../Genome-on-Diet-SNPs-Indels-SVs/GDiet_avx -t 40 --MD -ax map-hifi -Z 10 -W 2 -i 0.2 -k 19 -w $i -N 1 -r 1000 --vt_dis=800 --vt_nb_loc=10 --vt_df1=0.011 --vt_df2=0.15 --max_min_gap=4000 --vt_f=0.04 --sort=merge --frag=no -F200,1 --secondary=yes -a -o "HiFi/Genome-on-Diet-GRCh38-HiFi-stats_SV_k19w$i.sam" GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta ${FASTQpath}"

    sbatch --exclusive -w kratos4 -c 40 -J bwa_hg38 -o "HiFi/minimap2-GRCh38-HiFi-stats-command-out_k19w$i.txt" --wrap="/usr/bin/time -v --output "HiFi/minimap2-GRCh38-HiFi-stats-time-memory_k19w$i.txt" ./minimap2 -t 40 --MD -ax map-hifi -N 1 -a -k 19 -w $i --secondary=yes -o "HiFi/minimap2-GRCh38-HiFi-stats_k19w$i.sam" GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta ${FASTQpath}"
  
    sbatch --exclusive -w kratos4 -c 40 -J bwa_hg38 --wrap="samtools stats --threads 40 "HiFi/minimap2-GRCh38-HiFi-stats_k19w$i.sam" | grep ^SN | cut -f 2- > "HiFi/minimap2-GRCh38-HiFi-stats_k19w$i.txt""

done

#samcomp -m prim /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/HiFi/minimap2-GRCh38-HiFi-stats_k19w19.sam /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/HiFi/Genome-on-Diet-GRCh38-HiFi-stats_k19w19.sam 


cd HiFi/

grep 'seconds' minimap2-GRCh38-HiFi*.txt
grep 'mapped:' minimap2-GRCh38-HiFi*.txt | awk 'ORS=NR%3?FS:RS'
grep 'non-primary alignments' minimap2-GRCh38-HiFi*.txt
grep 'Maximum resident set size' minimap2-GRCh38-HiFi*.txt
grep 'PROFILING' minimap2-GRCh38-HiFi*.txt


grep 'seconds' Genome-on-Diet-GRCh38-HiFi*.txt
grep 'mapped:' Genome-on-Diet-GRCh38-HiFi*.txt | awk 'ORS=NR%3?FS:RS'
grep 'non-primary alignments' Genome-on-Diet-GRCh38-HiFi*.txt
grep 'Maximum resident set size' Genome-on-Diet-GRCh38-HiFi*.txt
grep 'PROFILING' Genome-on-Diet-GRCh38-HiFi*.txt


