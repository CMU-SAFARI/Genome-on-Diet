#!/bin/bash
set -e

#################################################
#Human Genome
cd /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/
Name='Human'
sbatch --exclusive -w kratos5 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_minimap2_${Name}-command-out.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_minimap2_${Name}-time-memory.txt" /home/alserm/minimap2-alser/minimap2-original-7June2022/minimap2 -t 40 -x map-hifi -I 300G -k 19 -w 16 --idx-no-seq -d /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_minimap2_${Name}.mmi /home/alserm/minimap2-alser/Data/GCF_000001405.40_GRCh38.p14_genomic.fna"

Pattern='11'
PatternLength='2'
Name='Human'
sbatch --exclusive -w kratos5 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}-command-out_${Pattern}.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}-time-memory_${Pattern}.txt" /home/alserm/minimap2-alser/src2/minimap2_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 --idx-no-seq -d /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}_${Pattern}.mmi /home/alserm/minimap2-alser/Data/GCF_000001405.40_GRCh38.p14_genomic.fna"

Pattern='10'
PatternLength='2'
Name='Human'
sbatch --exclusive -w kratos5 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}-command-out_${Pattern}.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}-time-memory_${Pattern}.txt" /home/alserm/minimap2-alser/src2/minimap2_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 --idx-no-seq -d /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}_${Pattern}.mmi /home/alserm/minimap2-alser/Data/GCF_000001405.40_GRCh38.p14_genomic.fna"


Pattern='100'
PatternLength='3'
Name='Human'
sbatch --exclusive -w kratos5 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}-command-out_${Pattern}.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}-time-memory_${Pattern}.txt" /home/alserm/minimap2-alser/src2/minimap2_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 --idx-no-seq -d /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}_${Pattern}.mmi /home/alserm/minimap2-alser/Data/GCF_000001405.40_GRCh38.p14_genomic.fna"


Pattern='110'
PatternLength='3'
Name='Human'
sbatch --exclusive -w kratos5 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}-command-out_${Pattern}.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}-time-memory_${Pattern}.txt" /home/alserm/minimap2-alser/src2/minimap2_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 --idx-no-seq -d /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}_${Pattern}.mmi /home/alserm/minimap2-alser/Data/GCF_000001405.40_GRCh38.p14_genomic.fna"


Pattern='1110'
PatternLength='4'
Name='Human'
sbatch --exclusive -w kratos5 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}-command-out_${Pattern}.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}-time-memory_${Pattern}.txt" /home/alserm/minimap2-alser/src2/minimap2_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 --idx-no-seq -d /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}_${Pattern}.mmi /home/alserm/minimap2-alser/Data/GCF_000001405.40_GRCh38.p14_genomic.fna"


#################################################
#Largest genome
cd /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/
Name='Pinus_taeda'
sbatch --exclusive -w kratos5 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_minimap2_${Name}-command-out.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_minimap2_${Name}-time-memory.txt" /home/alserm/minimap2-alser/minimap2-original-7June2022/minimap2 -t 40 -x map-hifi -I 300G -k 19 -w 16 --idx-no-seq -d /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_minimap2_${Name}.mmi /mnt/batty-shared/species/Pinus_taeda_loblolly-pine/GCA_001447015.2_Sugar_pine_JHU_assembly_genomic.fna"


Pattern='11'
PatternLength='2'
Name='Pinus_taeda'
sbatch --exclusive -w kratos4 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}-command-out_${Pattern}.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}-time-memory_${Pattern}.txt" /home/alserm/minimap2-alser/src2/minimap2_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 --idx-no-seq -d /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}_${Pattern}.mmi /mnt/batty-shared/species/Pinus_taeda_loblolly-pine/GCA_001447015.2_Sugar_pine_JHU_assembly_genomic.fna"

Pattern='10'
PatternLength='2'
Name='Pinus_taeda'
sbatch --exclusive -w kratos4 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}-command-out_${Pattern}.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}-time-memory_${Pattern}.txt" /home/alserm/minimap2-alser/src2/minimap2_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 --idx-no-seq -d /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}_${Pattern}.mmi /mnt/batty-shared/species/Pinus_taeda_loblolly-pine/GCA_001447015.2_Sugar_pine_JHU_assembly_genomic.fna"


Pattern='100'
PatternLength='3'
Name='Pinus_taeda'
sbatch --exclusive -w kratos4 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}-command-out_${Pattern}.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}-time-memory_${Pattern}.txt" /home/alserm/minimap2-alser/src2/minimap2_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 --idx-no-seq -d /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}_${Pattern}.mmi /mnt/batty-shared/species/Pinus_taeda_loblolly-pine/GCA_001447015.2_Sugar_pine_JHU_assembly_genomic.fna"


Pattern='110'
PatternLength='3'
Name='Pinus_taeda'
sbatch --exclusive -w kratos4 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}-command-out_${Pattern}.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}-time-memory_${Pattern}.txt" /home/alserm/minimap2-alser/src2/minimap2_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 --idx-no-seq -d /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}_${Pattern}.mmi /mnt/batty-shared/species/Pinus_taeda_loblolly-pine/GCA_001447015.2_Sugar_pine_JHU_assembly_genomic.fna"


Pattern='1110'
PatternLength='4'
Name='Pinus_taeda'
sbatch --exclusive -w kratos4 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}-command-out_${Pattern}.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}-time-memory_${Pattern}.txt" /home/alserm/minimap2-alser/src2/minimap2_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 --idx-no-seq -d /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}_${Pattern}.mmi /mnt/batty-shared/species/Pinus_taeda_loblolly-pine/GCA_001447015.2_Sugar_pine_JHU_assembly_genomic.fna"


#################################################
#RefSeq-100GB
cd /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/

Name='RefSeq-100GB'
sbatch --exclusive -w kratos1 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_minimap2_${Name}-command-out.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_minimap2_${Name}-time-memory.txt" /home/alserm/minimap2-alser/minimap2-original-7June2022/minimap2 -t 40 -x map-hifi -I 300G -k 19 -w 16 --idx-no-seq -d /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_minimap2_${Name}.mmi /home/alserm/minimap2-alser/Data/all-organisms.fna"

Pattern='11'
PatternLength='2'
Name='RefSeq-100GB'
sbatch --exclusive -w kratos2 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}-command-out_${Pattern}.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}-time-memory_${Pattern}.txt" /home/alserm/minimap2-alser/src2/minimap2_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 --idx-no-seq -d /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}_${Pattern}.mmi /home/alserm/minimap2-alser/Data/all-organisms.fna"

Pattern='10'
PatternLength='2'
Name='RefSeq-100GB'
sbatch --exclusive -w kratos3 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}-command-out_${Pattern}.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}-time-memory_${Pattern}.txt" /home/alserm/minimap2-alser/src2/minimap2_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 --idx-no-seq -d /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}_${Pattern}.mmi /home/alserm/minimap2-alser/Data/all-organisms.fna"


Pattern='100'
PatternLength='3'
Name='RefSeq-100GB'
sbatch --exclusive -w kratos4 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}-command-out_${Pattern}.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}-time-memory_${Pattern}.txt" /home/alserm/minimap2-alser/src2/minimap2_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 --idx-no-seq -d /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}_${Pattern}.mmi /home/alserm/minimap2-alser/Data/all-organisms.fna"



Pattern='110'
PatternLength='3'
Name='RefSeq-100GB'
sbatch --exclusive -w kratos5 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}-command-out_${Pattern}.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}-time-memory_${Pattern}.txt" /home/alserm/minimap2-alser/src2/minimap2_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 --idx-no-seq -d /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}_${Pattern}.mmi /home/alserm/minimap2-alser/Data/all-organisms.fna"


Pattern='1110'
PatternLength='4'
Name='RefSeq-100GB'
sbatch --exclusive -w kratos1 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}-command-out_${Pattern}.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}-time-memory_${Pattern}.txt" /home/alserm/minimap2-alser/src2/minimap2_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 --idx-no-seq -d /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}_${Pattern}.mmi /home/alserm/minimap2-alser/Data/all-organisms.fna"



##################################################
##RefSeq-300GB
#cd /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/
#Name='RefSeq-300GB'
#sbatch --exclusive -w kratos1 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_minimap2_${Name}-command-out.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_minimap2_${Name}-time-memory.txt" /home/alserm/minimap2-alser/minimap2-original-7June2022/minimap2 -t 40 -x map-hifi -I 300G -k 19 -w 16 --idx-no-seq -d /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_minimap2_${Name}.mmi /home/alserm/minimap2-alser/Data/all-organisms-300GB.fna"
#
#
#Pattern='11'
#PatternLength='2'
#Name='RefSeq-300GB'
#sbatch --exclusive -w kratos2 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}-command-out_${Pattern}.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}-time-memory_${Pattern}.txt" /home/alserm/minimap2-alser/src2/minimap2_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 --idx-no-seq -d /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}_${Pattern}.mmi /home/alserm/minimap2-alser/Data/all-organisms-300GB.fna"
#
#Pattern='10'
#PatternLength='2'
#Name='RefSeq-300GB'
#sbatch --exclusive -w kratos3 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}-command-out_${Pattern}.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}-time-memory_${Pattern}.txt" /home/alserm/minimap2-alser/src2/minimap2_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 --idx-no-seq -d /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}_${Pattern}.mmi /home/alserm/minimap2-alser/Data/all-organisms-300GB.fna"
#
#
#Pattern='100'
#PatternLength='3'
#Name='RefSeq-300GB'
#sbatch --exclusive -w kratos4 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}-command-out_${Pattern}.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}-time-memory_${Pattern}.txt" /home/alserm/minimap2-alser/src2/minimap2_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 --idx-no-seq -d /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}_${Pattern}.mmi /home/alserm/minimap2-alser/Data/all-organisms-300GB.fna"
#
#
#
#Pattern='110'
#PatternLength='3'
#Name='RefSeq-300GB'
#sbatch --exclusive -w kratos5 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}-command-out_${Pattern}.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}-time-memory_${Pattern}.txt" /home/alserm/minimap2-alser/src2/minimap2_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 --idx-no-seq -d /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}_${Pattern}.mmi /home/alserm/minimap2-alser/Data/all-organisms-300GB.fna"
#
#
#Pattern='1110'
#PatternLength='4'
#Name='RefSeq-300GB'
#sbatch --exclusive -w kratos6 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}-command-out_${Pattern}.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}-time-memory_${Pattern}.txt" /home/alserm/minimap2-alser/src2/minimap2_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 --idx-no-seq -d /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/Indexing_${Name}_${Pattern}.mmi /home/alserm/minimap2-alser/Data/all-organisms-300GB.fna"



######### Print statistics
cd /home/alserm/minimap2-alser/Use-cases/kmer-counting/ContainmentIndexing/
grep 'seconds' *time-memory*.txt
grep 'kbytes' *time-memory*.txt
du *.mmi