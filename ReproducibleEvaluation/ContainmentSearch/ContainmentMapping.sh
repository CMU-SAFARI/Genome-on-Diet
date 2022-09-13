#!/bin/bash
set -e

#/home/alserm/minimap2-alser/Use-cases/kmer-counting/Genome-on-Diet-16Aug2022-ConatinmentIndexing/minimap2_avx -t 40 -ax map-hifi -Z 10 -W 2 -i 0.2 -k 19 -w 16 -N 1 -r 0.04,400,800 -n 0.8,0 --AF_max_loc 1 --sort=merge --frag=no -F200,1 --secondary=no -a /home/alserm/minimap2-alser/Data/GRCh38.p14_NC_000001.11.fna /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/FASTQ-files/Human-HiFi1.fastq | awk -F'\t' '{print $1 "\t" $14}' | head -n 100


#################################################
#Human Genome
cd /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/
Name='Human'
sbatch --exclusive -w kratos3 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_minimap2_${Name}-command-out.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_minimap2_${Name}-time-memory.txt" /home/alserm/minimap2-alser/minimap2-original-7June2022/minimap2 -t 40 -x map-hifi -I 300G -k 19 -w 16 /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_minimap2_${Name}.sam /home/alserm/minimap2-alser/Data/GCF_000001405.40_GRCh38.p14_genomic.fna /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/FASTQ-files/Human-HiFi1.fastq"

Pattern='11'
PatternLength='2'
Name='Human'
sbatch --exclusive -w kratos5 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}-command-out_${Pattern}.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}-time-memory_${Pattern}.txt" /home/alserm/minimap2-alser/Use-cases/kmer-counting/Genome-on-Diet-16Aug2022-ConatinmentIndexing/minimap2_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 -i 0.2 -N 1 -r 0.04,400,800 -n 0.8,0 --AF_max_loc 1 --sort=merge --frag=no -F200,1 --secondary=no -a -o /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}_${Pattern}.sam /home/alserm/minimap2-alser/Data/GCF_000001405.40_GRCh38.p14_genomic.fna /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/FASTQ-files/Human-HiFi1.fastq"

Pattern='10'
PatternLength='2'
Name='Human'
sbatch --exclusive -w kratos5 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}-command-out_${Pattern}.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}-time-memory_${Pattern}.txt" /home/alserm/minimap2-alser/Use-cases/kmer-counting/Genome-on-Diet-16Aug2022-ConatinmentIndexing/minimap2_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 -i 0.2 -N 1 -r 0.04,400,800 -n 0.8,0 --AF_max_loc 1 --sort=merge --frag=no -F200,1 --secondary=no -a -o /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}_${Pattern}.sam /home/alserm/minimap2-alser/Data/GCF_000001405.40_GRCh38.p14_genomic.fna /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/FASTQ-files/Human-HiFi1.fastq"


Pattern='100'
PatternLength='3'
Name='Human'
sbatch --exclusive -w kratos5 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}-command-out_${Pattern}.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}-time-memory_${Pattern}.txt" /home/alserm/minimap2-alser/Use-cases/kmer-counting/Genome-on-Diet-16Aug2022-ConatinmentIndexing/minimap2_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 -i 0.2 -N 1 -r 0.04,400,800 -n 0.8,0 --AF_max_loc 1 --sort=merge --frag=no -F200,1 --secondary=no -a -o /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}_${Pattern}.sam /home/alserm/minimap2-alser/Data/GCF_000001405.40_GRCh38.p14_genomic.fna /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/FASTQ-files/Human-HiFi1.fastq"


Pattern='110'
PatternLength='3'
Name='Human'
sbatch --exclusive -w kratos5 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}-command-out_${Pattern}.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}-time-memory_${Pattern}.txt" /home/alserm/minimap2-alser/Use-cases/kmer-counting/Genome-on-Diet-16Aug2022-ConatinmentIndexing/minimap2_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 -i 0.2 -N 1 -r 0.04,400,800 -n 0.8,0 --AF_max_loc 1 --sort=merge --frag=no -F200,1 --secondary=no -a -o /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}_${Pattern}.sam /home/alserm/minimap2-alser/Data/GCF_000001405.40_GRCh38.p14_genomic.fna /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/FASTQ-files/Human-HiFi1.fastq"


Pattern='1110'
PatternLength='4'
Name='Human'
sbatch --exclusive -w kratos5 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}-command-out_${Pattern}.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}-time-memory_${Pattern}.txt" /home/alserm/minimap2-alser/Use-cases/kmer-counting/Genome-on-Diet-16Aug2022-ConatinmentIndexing/minimap2_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 -i 0.2 -N 1 -r 0.04,400,800 -n 0.8,0 --AF_max_loc 1 --sort=merge --frag=no -F200,1 --secondary=no -a -o /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}_${Pattern}.sam /home/alserm/minimap2-alser/Data/GCF_000001405.40_GRCh38.p14_genomic.fna /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/FASTQ-files/Human-HiFi1.fastq"


#################################################
#Largest genome
cd /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/
Name='Pinus_taeda'
sbatch --exclusive -w kratos3 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_minimap2_${Name}-command-out.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_minimap2_${Name}-time-memory.txt" /home/alserm/minimap2-alser/minimap2-original-7June2022/minimap2 -t 40 -x map-hifi -I 300G -k 19 -w 16 --idx-no-seq -d /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_minimap2_${Name}.sam /mnt/batty-shared/species/Pinus_taeda_loblolly-pine/GCA_001447015.2_Sugar_pine_JHU_assembly_genomic.fna /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/FASTQ-files/Pinus_taeda-HiFi1.fastq"


Pattern='11'
PatternLength='2'
Name='Pinus_taeda'
sbatch --exclusive -w kratos4 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}-command-out_${Pattern}.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}-time-memory_${Pattern}.txt" /home/alserm/minimap2-alser/Use-cases/kmer-counting/Genome-on-Diet-16Aug2022-ConatinmentIndexing/minimap2_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 -i 0.2 -N 1 -r 0.04,400,800 -n 0.8,0 --AF_max_loc 1 --sort=merge --frag=no -F200,1 --secondary=no -a -o /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}_${Pattern}.sam /mnt/batty-shared/species/Pinus_taeda_loblolly-pine/GCA_001447015.2_Sugar_pine_JHU_assembly_genomic.fna /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/FASTQ-files/Pinus_taeda-HiFi1.fastq"

Pattern='10'
PatternLength='2'
Name='Pinus_taeda'
sbatch --exclusive -w kratos4 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}-command-out_${Pattern}.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}-time-memory_${Pattern}.txt" /home/alserm/minimap2-alser/Use-cases/kmer-counting/Genome-on-Diet-16Aug2022-ConatinmentIndexing/minimap2_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 -i 0.2 -N 1 -r 0.04,400,800 -n 0.8,0 --AF_max_loc 1 --sort=merge --frag=no -F200,1 --secondary=no -a -o /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}_${Pattern}.sam /mnt/batty-shared/species/Pinus_taeda_loblolly-pine/GCA_001447015.2_Sugar_pine_JHU_assembly_genomic.fna"


Pattern='100'
PatternLength='3'
Name='Pinus_taeda'
sbatch --exclusive -w kratos4 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}-command-out_${Pattern}.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}-time-memory_${Pattern}.txt" /home/alserm/minimap2-alser/Use-cases/kmer-counting/Genome-on-Diet-16Aug2022-ConatinmentIndexing/minimap2_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 -i 0.2 -N 1 -r 0.04,400,800 -n 0.8,0 --AF_max_loc 1 --sort=merge --frag=no -F200,1 --secondary=no -a -o /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}_${Pattern}.sam /mnt/batty-shared/species/Pinus_taeda_loblolly-pine/GCA_001447015.2_Sugar_pine_JHU_assembly_genomic.fna /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/FASTQ-files/Pinus_taeda-HiFi1.fastq"


Pattern='110'
PatternLength='3'
Name='Pinus_taeda'
sbatch --exclusive -w kratos4 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}-command-out_${Pattern}.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}-time-memory_${Pattern}.txt" /home/alserm/minimap2-alser/Use-cases/kmer-counting/Genome-on-Diet-16Aug2022-ConatinmentIndexing/minimap2_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 -i 0.2 -N 1 -r 0.04,400,800 -n 0.8,0 --AF_max_loc 1 --sort=merge --frag=no -F200,1 --secondary=no -a -o /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}_${Pattern}.sam /mnt/batty-shared/species/Pinus_taeda_loblolly-pine/GCA_001447015.2_Sugar_pine_JHU_assembly_genomic.fna /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/FASTQ-files/Pinus_taeda-HiFi1.fastq"


Pattern='1110'
PatternLength='4'
Name='Pinus_taeda'
sbatch --exclusive -w kratos4 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}-command-out_${Pattern}.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}-time-memory_${Pattern}.txt" /home/alserm/minimap2-alser/Use-cases/kmer-counting/Genome-on-Diet-16Aug2022-ConatinmentIndexing/minimap2_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 -i 0.2 -N 1 -r 0.04,400,800 -n 0.8,0 --AF_max_loc 1 --sort=merge --frag=no -F200,1 --secondary=no -a -o /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}_${Pattern}.sam /mnt/batty-shared/species/Pinus_taeda_loblolly-pine/GCA_001447015.2_Sugar_pine_JHU_assembly_genomic.fna /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/FASTQ-files/Pinus_taeda-HiFi1.fastq"


#################################################
#RefSeq-100GB
cd /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/

Name='RefSeq-100GB'
sbatch --exclusive -w kratos1 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_minimap2_${Name}-command-out.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_minimap2_${Name}-time-memory.txt" /home/alserm/minimap2-alser/minimap2-original-7June2022/minimap2 -t 40 -x map-hifi -I 300G -k 19 -w 16 --idx-no-seq -d /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_minimap2_${Name}.sam /home/alserm/minimap2-alser/Data/all-organisms.fna /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/FASTQ-files/RefSeq-100GB-HiFi1.fastq"

Pattern='11'
PatternLength='2'
Name='RefSeq-100GB'
sbatch --exclusive -w kratos2 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}-command-out_${Pattern}.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}-time-memory_${Pattern}.txt" /home/alserm/minimap2-alser/Use-cases/kmer-counting/Genome-on-Diet-16Aug2022-ConatinmentIndexing/minimap2_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 -i 0.2 -N 1 -r 0.04,400,800 -n 0.8,0 --AF_max_loc 1 --sort=merge --frag=no -F200,1 --secondary=no -a -o /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}_${Pattern}.sam /home/alserm/minimap2-alser/Data/all-organisms.fna /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/FASTQ-files/RefSeq-100GB-HiFi1.fastq"

Pattern='10'
PatternLength='2'
Name='RefSeq-100GB'
sbatch --exclusive -w kratos3 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}-command-out_${Pattern}.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}-time-memory_${Pattern}.txt" /home/alserm/minimap2-alser/Use-cases/kmer-counting/Genome-on-Diet-16Aug2022-ConatinmentIndexing/minimap2_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 -i 0.2 -N 1 -r 0.04,400,800 -n 0.8,0 --AF_max_loc 1 --sort=merge --frag=no -F200,1 --secondary=no -a -o /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}_${Pattern}.sam /home/alserm/minimap2-alser/Data/all-organisms.fna /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/FASTQ-files/RefSeq-100GB-HiFi1.fastq"


Pattern='100'
PatternLength='3'
Name='RefSeq-100GB'
sbatch --exclusive -w kratos4 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}-command-out_${Pattern}.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}-time-memory_${Pattern}.txt" /home/alserm/minimap2-alser/Use-cases/kmer-counting/Genome-on-Diet-16Aug2022-ConatinmentIndexing/minimap2_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 -i 0.2 -N 1 -r 0.04,400,800 -n 0.8,0 --AF_max_loc 1 --sort=merge --frag=no -F200,1 --secondary=no -a -o /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}_${Pattern}.sam /home/alserm/minimap2-alser/Data/all-organisms.fna /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/FASTQ-files/RefSeq-100GB-HiFi1.fastq"



Pattern='110'
PatternLength='3'
Name='RefSeq-100GB'
sbatch --exclusive -w kratos5 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}-command-out_${Pattern}.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}-time-memory_${Pattern}.txt" /home/alserm/minimap2-alser/Use-cases/kmer-counting/Genome-on-Diet-16Aug2022-ConatinmentIndexing/minimap2_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 -i 0.2 -N 1 -r 0.04,400,800 -n 0.8,0 --AF_max_loc 1 --sort=merge --frag=no -F200,1 --secondary=no -a -o /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}_${Pattern}.sam /home/alserm/minimap2-alser/Data/all-organisms.fna /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/FASTQ-files/RefSeq-100GB-HiFi1.fastq"


Pattern='1110'
PatternLength='4'
Name='RefSeq-100GB'
sbatch --exclusive -w kratos1 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}-command-out_${Pattern}.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}-time-memory_${Pattern}.txt" /home/alserm/minimap2-alser/Use-cases/kmer-counting/Genome-on-Diet-16Aug2022-ConatinmentIndexing/minimap2_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 -i 0.2 -N 1 -r 0.04,400,800 -n 0.8,0 --AF_max_loc 1 --sort=merge --frag=no -F200,1 --secondary=no -a -o /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}_${Pattern}.sam /home/alserm/minimap2-alser/Data/all-organisms.fna /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/FASTQ-files/RefSeq-100GB-HiFi1.fastq"



##################################################
##RefSeq-300GB
#cd /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/
#Name='RefSeq-300GB'
#sbatch --exclusive -w kratos1 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_minimap2_${Name}-command-out.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_minimap2_${Name}-time-memory.txt" /home/alserm/minimap2-alser/minimap2-original-7June2022/minimap2 -t 40 -x map-hifi -I 300G -k 19 -w 16 --idx-no-seq -d /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_minimap2_${Name}.sam /home/alserm/minimap2-alser/Data/all-organisms-300GB.fna"
#
#
#Pattern='11'
#PatternLength='2'
#Name='RefSeq-300GB'
#sbatch --exclusive -w kratos2 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}-command-out_${Pattern}.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}-time-memory_${Pattern}.txt" /home/alserm/minimap2-alser/Use-cases/kmer-counting/Genome-on-Diet-16Aug2022-ConatinmentIndexing/minimap2_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 -i 0.2 -N 1 -r 0.04,400,800 -n 0.8,0 --AF_max_loc 1 --sort=merge --frag=no -F200,1 --secondary=no -a -o /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}_${Pattern}.sam /home/alserm/minimap2-alser/Data/all-organisms-300GB.fna"
#
#Pattern='10'
#PatternLength='2'
#Name='RefSeq-300GB'
#sbatch --exclusive -w kratos3 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}-command-out_${Pattern}.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}-time-memory_${Pattern}.txt" /home/alserm/minimap2-alser/Use-cases/kmer-counting/Genome-on-Diet-16Aug2022-ConatinmentIndexing/minimap2_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 -i 0.2 -N 1 -r 0.04,400,800 -n 0.8,0 --AF_max_loc 1 --sort=merge --frag=no -F200,1 --secondary=no -a -o /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}_${Pattern}.sam /home/alserm/minimap2-alser/Data/all-organisms-300GB.fna"
#
#
#Pattern='100'
#PatternLength='3'
#Name='RefSeq-300GB'
#sbatch --exclusive -w kratos4 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}-command-out_${Pattern}.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}-time-memory_${Pattern}.txt" /home/alserm/minimap2-alser/Use-cases/kmer-counting/Genome-on-Diet-16Aug2022-ConatinmentIndexing/minimap2_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 -i 0.2 -N 1 -r 0.04,400,800 -n 0.8,0 --AF_max_loc 1 --sort=merge --frag=no -F200,1 --secondary=no -a -o /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}_${Pattern}.sam /home/alserm/minimap2-alser/Data/all-organisms-300GB.fna"
#
#
#
#Pattern='110'
#PatternLength='3'
#Name='RefSeq-300GB'
#sbatch --exclusive -w kratos5 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}-command-out_${Pattern}.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}-time-memory_${Pattern}.txt" /home/alserm/minimap2-alser/Use-cases/kmer-counting/Genome-on-Diet-16Aug2022-ConatinmentIndexing/minimap2_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 -i 0.2 -N 1 -r 0.04,400,800 -n 0.8,0 --AF_max_loc 1 --sort=merge --frag=no -F200,1 --secondary=no -a -o /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}_${Pattern}.sam /home/alserm/minimap2-alser/Data/all-organisms-300GB.fna"
#
#
#Pattern='1110'
#PatternLength='4'
#Name='RefSeq-300GB'
#sbatch --exclusive -w kratos6 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}-command-out_${Pattern}.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}-time-memory_${Pattern}.txt" /home/alserm/minimap2-alser/Use-cases/kmer-counting/Genome-on-Diet-16Aug2022-ConatinmentIndexing/minimap2_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 -i 0.2 -N 1 -r 0.04,400,800 -n 0.8,0 --AF_max_loc 1 --sort=merge --frag=no -F200,1 --secondary=no -a -o /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Mapping_${Name}_${Pattern}.sam /home/alserm/minimap2-alser/Data/all-organisms-300GB.fna"



######### Print statistics
cd /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/
grep 'seconds' *time-memory*.txt
grep 'kbytes' *time-memory*.txt
du *.sam
