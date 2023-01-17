#!/bin/bash
#set -e

cd /home/alserm/minimap2-alser/
#source ~/.bashrc
#conda activate Metalign


FASTQpath= 'HG002_ONT-UL_GIAB_20200204_1000filtered_2Mreads_SplitedInto30K.fastq'
FASTQpath1='HG002_ONT-UL_GIAB_20200204_1000filtered_2Mreads.fastq'



for i in {6..18..2}
do
    /usr/bin/time -v python UltraLongReadSplitter.py ${FASTQpath1} 30000 > ${FASTQpath}
    
    sbatch --exclusive -w kratos1 -c 40 -J bwa_hg38 -o "ONT/Genome-on-Diet-GRCh38-ONT-stats-command-out_SV_k15w$i.txt" --wrap="/usr/bin/time -v --output "ONT/Genome-on-Diet-GRCh38-ONT-stats-time-memory_SV_k15w$i.txt" GDiet-LongReads/GDiet_avx -t 40 --MD -ax map-ont -Z 10 -W 2 -i 0.2 -k 15 -w $i -N 1 -r 1300 --vt_dis=1000 --vt_nb_loc=3 --vt_df1=0.007 --vt_df2=0.007 --max_min_gap=4000 --vt_f=0.04 -s 35000 --vt_cov 0.3 --sort=merge --frag=no -F200,1 --secondary=yes -o "ONT/Genome-on-Diet-GRCh38-ONT-stats_SV_k15w$i.sam" GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta ${FASTQpath}"
    
    ##The line below is when you want faster version with almost the same performance for variant detection of minimap2
    #sbatch --exclusive -w kratos1 -c 40 -J bwa_hg38 -o "ONT/Genome-on-Diet-GRCh38-ONT-stats-command-out_SV_k15w$i.txt" --wrap="/usr/bin/time -v --output "ONT/Genome-on-Diet-GRCh38-ONT-stats-time-memory_SV_k15w$i.txt" GDiet-LongReads/GDiet_avx -t 40 --MD -ax map-ont -Z 10 -W 2 -i 0.2 -k 15 -w $i -N 1 -r 1300 --vt_dis=1000 --vt_nb_loc=3 --vt_df1=0.016 --vt_df2=0.01 --max_min_gap=4000 --vt_f=0.04 -s 35000 --vt_cov 0.3 --sort=merge --frag=no -F200,1 --secondary=yes -o "ONT/Genome-on-Diet-GRCh38-ONT-stats_SV_k15w$i.sam" GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta ${FASTQpath}"

    sbatch --exclusive -w kratos1 -c 40 -J bwa_hg38 --wrap="samtools stats --threads 40 "ONT/Genome-on-Diet-GRCh38-ONT-stats_SV_k15w$i.sam" | grep ^SN | cut -f 2- > "ONT/Genome-on-Diet-GRCh38-ONT-stats_SV_k15w$i.txt""


    sbatch --exclusive -w kratos3 -c 40 -J bwa_hg38 -o "ONT/minimap2-GRCh38-ONT-stats-command-out_k15w$i.txt" --wrap="/usr/bin/time -v --output "ONT/minimap2-GRCh38-ONT-stats-time-memory_k15w$i.txt" ./minimap2 -t 40 --MD -ax map-ont -N 1 -a -k 15 -w $i --secondary=yes -o "ONT/minimap2-GRCh38-ONT-stats_k15w$i.sam" GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta ${FASTQpath1}"

    sbatch --exclusive -w kratos3 -c 40 -J bwa_hg38 --wrap="samtools stats --threads 40 "ONT/minimap2-GRCh38-ONT-stats_k15w$i.sam" | grep ^SN | cut -f 2- > "ONT/minimap2-GRCh38-ONT-stats_k15w$i.txt""

done


cd ONT/
grep 'seconds' minimap2-GRCh38-ONT*.txt
grep 'mapped:' minimap2-GRCh38-ONT*.txt | awk 'ORS=NR%3?FS:RS'
grep 'non-primary alignments' minimap2-GRCh38-ONT*.txt
grep 'Maximum resident set size' minimap2-GRCh38-ONT*.txt
grep 'PROFILING' minimap2-GRCh38-ONT*.txt


grep 'seconds' Genome-on-Diet-GRCh38-ONT*.txt
grep 'mapped:' Genome-on-Diet-GRCh38-ONT*.txt | awk 'ORS=NR%3?FS:RS'
grep 'non-primary alignments' Genome-on-Diet-GRCh38-ONT*.txt
grep 'Maximum resident set size' Genome-on-Diet-GRCh38-ONT*.txt
grep 'PROFILING' Genome-on-Diet-GRCh38-ONT*.txt

