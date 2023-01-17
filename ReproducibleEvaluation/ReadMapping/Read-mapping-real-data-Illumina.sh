#!/bin/bash
#set -e

FASTQpath='D1_S1_L001_R1_001-017.fastq'


for i in {5..19..2}
do

    sbatch --exclusive -w kratos1 -c 40 -J bwa_hg38 -o "Illumina/Genome-on-Diet-GRCh38-Illumina-stats-command-out_INDEL_SNP_k21w$i.txt" --wrap="/usr/bin/time -v --output "Illumina/Genome-on-Diet-GRCh38-Illumina-stats-time-memory_INDEL_SNP_k21w$i.txt" GDiet-ShortReads/GDiet_avx --MD -t 40 -ax sr -Z 10 -W 2 -i 2 -k 21 -w $i -N 1 -r 0.05,150,200 -n 0.95,0.3 -s 100 --AF_max_loc 2 --secondary=yes -a -o "Illumina/Genome-on-Diet-GRCh38-Illumina-stats_INDEL_SNP_k21w$i.sam" GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta ${FASTQpath}"

    sbatch --exclusive -w kratos1 -c 40 -J bwa_hg38 --wrap="samtools stats --threads 40 "Illumina/Genome-on-Diet-GRCh38-Illumina-stats_INDEL_SNP_k21w$i.sam" | grep ^SN | cut -f 2- > "Illumina/Genome-on-Diet-GRCh38-Illumina-stats_INDEL_SNP_k21w$i.txt""

    sbatch --exclusive -w kratos1 -c 40 -J bwa_hg38 -o "Illumina/minimap2-GRCh38-Illumina-stats-command-out_k21w$i.txt" --wrap="/usr/bin/time -v --output "Illumina/minimap2-GRCh38-Illumina-stats-time-memory_k21w$i.txt" ./minimap2 -t 40 --MD -ax sr -N 1 -a -k 21 -w $i --secondary=yes -o "Illumina/minimap2-GRCh38-Illumina-stats_k21w$i.sam" GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta ${FASTQpath}"

    sbatch --exclusive -w kratos1 -c 40 -J bwa_hg38 --wrap="samtools stats --threads 40 "Illumina/minimap2-GRCh38-Illumina-stats_k21w$i.sam" | grep ^SN | cut -f 2- > "Illumina/minimap2-GRCh38-Illumina-stats_k21w$i.txt""

done


cd Illumina/
grep 'seconds' minimap2-GRCh38-Illumina*.txt
grep 'mapped:' minimap2-GRCh38-Illumina*.txt | awk 'ORS=NR%3?FS:RS'
grep 'non-primary alignments' minimap2-GRCh38-Illumina*.txt
grep 'Maximum resident set size' minimap2-GRCh38-Illumina*.txt
grep 'PROFILING' minimap2-GRCh38-Illumina*.txt


grep 'seconds' Genome-on-Diet-GRCh38-Illumina*.txt
grep 'mapped:' Genome-on-Diet-GRCh38-Illumina*.txt | awk 'ORS=NR%3?FS:RS'
grep 'non-primary alignments' Genome-on-Diet-GRCh38-Illumina*.txt
grep 'Maximum resident set size' Genome-on-Diet-GRCh38-Illumina*.txt
grep 'PROFILING' Genome-on-Diet-GRCh38-Illumina*.txt
