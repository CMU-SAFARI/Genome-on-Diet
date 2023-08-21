#!/bin/bash
set -e

#################################################
#Human Genome
Name='Human'
./GDiet_avx -t 40 -x map-hifi -I 300G -k 19 -w 16 /Data/Mapping_GDiet_${Name}.sam /Data/GCF_000001405.40_GRCh38.p14_genomic.fna /FASTQ-files/Human-HiFi1.fastq

Pattern='11'
PatternLength='2'
Name='Human'
./GDiet_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 -i 0.2 -N 1 -r 0.04,400,800 -n 0.8,0 --AF_max_loc 1 --sort=merge --frag=no -F200,1 --secondary=no -a -o /Data/Mapping_${Name}_${Pattern}.sam /Data/GCF_000001405.40_GRCh38.p14_genomic.fna /FASTQ-files/Human-HiFi1.fastq

Pattern='10'
PatternLength='2'
Name='Human'
./GDiet_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 -i 0.2 -N 1 -r 0.04,400,800 -n 0.8,0 --AF_max_loc 1 --sort=merge --frag=no -F200,1 --secondary=no -a -o /Data/Mapping_${Name}_${Pattern}.sam /Data/GCF_000001405.40_GRCh38.p14_genomic.fna /FASTQ-files/Human-HiFi1.fastq


Pattern='100'
PatternLength='3'
Name='Human'
./GDiet_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 -i 0.2 -N 1 -r 0.04,400,800 -n 0.8,0 --AF_max_loc 1 --sort=merge --frag=no -F200,1 --secondary=no -a -o /Data/Mapping_${Name}_${Pattern}.sam /Data/GCF_000001405.40_GRCh38.p14_genomic.fna /FASTQ-files/Human-HiFi1.fastq


Pattern='110'
PatternLength='3'
Name='Human'
./GDiet_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 -i 0.2 -N 1 -r 0.04,400,800 -n 0.8,0 --AF_max_loc 1 --sort=merge --frag=no -F200,1 --secondary=no -a -o /Data/Mapping_${Name}_${Pattern}.sam /Data/GCF_000001405.40_GRCh38.p14_genomic.fna /FASTQ-files/Human-HiFi1.fastq


Pattern='1110'
PatternLength='4'
Name='Human'
./GDiet_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 -i 0.2 -N 1 -r 0.04,400,800 -n 0.8,0 --AF_max_loc 1 --sort=merge --frag=no -F200,1 --secondary=no -a -o /Data/Mapping_${Name}_${Pattern}.sam /Data/GCF_000001405.40_GRCh38.p14_genomic.fna /FASTQ-files/Human-HiFi1.fastq


#################################################
#Largest genome
Name='Pinus_taeda'
./GDiet_avx -t 40 -x map-hifi -I 300G -k 19 -w 16 --idx-no-seq -d /Data/Mapping_GDiet_${Name}.sam /Data/Pinus_taeda_loblolly-pine/GCA_001447015.2_Sugar_pine_JHU_assembly_genomic.fna /FASTQ-files/Pinus_taeda-HiFi1.fastq


Pattern='11'
PatternLength='2'
Name='Pinus_taeda'
./GDiet_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 -i 0.2 -N 1 -r 0.04,400,800 -n 0.8,0 --AF_max_loc 1 --sort=merge --frag=no -F200,1 --secondary=no -a -o /Data/Mapping_${Name}_${Pattern}.sam /Data/Pinus_taeda_loblolly-pine/GCA_001447015.2_Sugar_pine_JHU_assembly_genomic.fna /FASTQ-files/Pinus_taeda-HiFi1.fastq

Pattern='10'
PatternLength='2'
Name='Pinus_taeda'
./GDiet_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 -i 0.2 -N 1 -r 0.04,400,800 -n 0.8,0 --AF_max_loc 1 --sort=merge --frag=no -F200,1 --secondary=no -a -o /Data/Mapping_${Name}_${Pattern}.sam /Data/Pinus_taeda_loblolly-pine/GCA_001447015.2_Sugar_pine_JHU_assembly_genomic.fna /FASTQ-files/Pinus_taeda-HiFi1.fastq


Pattern='100'
PatternLength='3'
Name='Pinus_taeda'
./GDiet_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 -i 0.2 -N 1 -r 0.04,400,800 -n 0.8,0 --AF_max_loc 1 --sort=merge --frag=no -F200,1 --secondary=no -a -o /Data/Mapping_${Name}_${Pattern}.sam /Data/Pinus_taeda_loblolly-pine/GCA_001447015.2_Sugar_pine_JHU_assembly_genomic.fna /FASTQ-files/Pinus_taeda-HiFi1.fastq


Pattern='110'
PatternLength='3'
Name='Pinus_taeda'
./GDiet_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 -i 0.2 -N 1 -r 0.04,400,800 -n 0.8,0 --AF_max_loc 1 --sort=merge --frag=no -F200,1 --secondary=no -a -o /Data/Mapping_${Name}_${Pattern}.sam /Data/Pinus_taeda_loblolly-pine/GCA_001447015.2_Sugar_pine_JHU_assembly_genomic.fna /FASTQ-files/Pinus_taeda-HiFi1.fastq


Pattern='1110'
PatternLength='4'
Name='Pinus_taeda'
./GDiet_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 -i 0.2 -N 1 -r 0.04,400,800 -n 0.8,0 --AF_max_loc 1 --sort=merge --frag=no -F200,1 --secondary=no -a -o /Data/Mapping_${Name}_${Pattern}.sam /Data/Pinus_taeda_loblolly-pine/GCA_001447015.2_Sugar_pine_JHU_assembly_genomic.fna /FASTQ-files/Pinus_taeda-HiFi1.fastq


#################################################
#RefSeq-100GB

Name='RefSeq-100GB'
./GDiet_avx -t 40 -x map-hifi -I 300G -k 19 -w 16 --idx-no-seq -d /Data/Mapping_GDiet_${Name}.sam /Data/all-organisms.fna /FASTQ-files/RefSeq-100GB-HiFi1.fastq

Pattern='11'
PatternLength='2'
Name='RefSeq-100GB'
./GDiet_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 -i 0.2 -N 1 -r 0.04,400,800 -n 0.8,0 --AF_max_loc 1 --sort=merge --frag=no -F200,1 --secondary=no -a -o /Data/Mapping_${Name}_${Pattern}.sam /Data/all-organisms.fna /FASTQ-files/RefSeq-100GB-HiFi1.fastq

Pattern='10'
PatternLength='2'
Name='RefSeq-100GB'
./GDiet_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 -i 0.2 -N 1 -r 0.04,400,800 -n 0.8,0 --AF_max_loc 1 --sort=merge --frag=no -F200,1 --secondary=no -a -o /Data/Mapping_${Name}_${Pattern}.sam /Data/all-organisms.fna /FASTQ-files/RefSeq-100GB-HiFi1.fastq


Pattern='100'
PatternLength='3'
Name='RefSeq-100GB'
./GDiet_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 -i 0.2 -N 1 -r 0.04,400,800 -n 0.8,0 --AF_max_loc 1 --sort=merge --frag=no -F200,1 --secondary=no -a -o /Data/Mapping_${Name}_${Pattern}.sam /Data/all-organisms.fna /FASTQ-files/RefSeq-100GB-HiFi1.fastq



Pattern='110'
PatternLength='3'
Name='RefSeq-100GB'
./GDiet_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 -i 0.2 -N 1 -r 0.04,400,800 -n 0.8,0 --AF_max_loc 1 --sort=merge --frag=no -F200,1 --secondary=no -a -o /Data/Mapping_${Name}_${Pattern}.sam /Data/all-organisms.fna /FASTQ-files/RefSeq-100GB-HiFi1.fastq


Pattern='1110'
PatternLength='4'
Name='RefSeq-100GB'
./GDiet_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 -i 0.2 -N 1 -r 0.04,400,800 -n 0.8,0 --AF_max_loc 1 --sort=merge --frag=no -F200,1 --secondary=no -a -o /Data/Mapping_${Name}_${Pattern}.sam /Data/all-organisms.fna /FASTQ-files/RefSeq-100GB-HiFi1.fastq


#################################################
#RefSeq-200GB

Name='RefSeq-200GB'
./GDiet_avx -t 40 -x map-hifi -I 300G -k 19 -w 16 --idx-no-seq -d /Data/Mapping_GDiet_${Name}.sam /Data/all-organisms.fna /FASTQ-files/RefSeq-200GB-HiFi1.fastq

Pattern='11'
PatternLength='2'
Name='RefSeq-200GB'
./GDiet_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 -i 0.2 -N 1 -r 0.04,400,800 -n 0.8,0 --AF_max_loc 1 --sort=merge --frag=no -F200,1 --secondary=no -a -o /Data/Mapping_${Name}_${Pattern}.sam /Data/all-organisms.fna /FASTQ-files/RefSeq-200GB-HiFi1.fastq

Pattern='10'
PatternLength='2'
Name='RefSeq-200GB'
./GDiet_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 -i 0.2 -N 1 -r 0.04,400,800 -n 0.8,0 --AF_max_loc 1 --sort=merge --frag=no -F200,1 --secondary=no -a -o /Data/Mapping_${Name}_${Pattern}.sam /Data/all-organisms.fna /FASTQ-files/RefSeq-200GB-HiFi1.fastq


Pattern='100'
PatternLength='3'
Name='RefSeq-200GB'
./GDiet_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 -i 0.2 -N 1 -r 0.04,400,800 -n 0.8,0 --AF_max_loc 1 --sort=merge --frag=no -F200,1 --secondary=no -a -o /Data/Mapping_${Name}_${Pattern}.sam /Data/all-organisms.fna /FASTQ-files/RefSeq-200GB-HiFi1.fastq



Pattern='110'
PatternLength='3'
Name='RefSeq-200GB'
./GDiet_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 -i 0.2 -N 1 -r 0.04,400,800 -n 0.8,0 --AF_max_loc 1 --sort=merge --frag=no -F200,1 --secondary=no -a -o /Data/Mapping_${Name}_${Pattern}.sam /Data/all-organisms.fna /FASTQ-files/RefSeq-200GB-HiFi1.fastq


Pattern='1110'
PatternLength='4'
Name='RefSeq-200GB'
./GDiet_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 -i 0.2 -N 1 -r 0.04,400,800 -n 0.8,0 --AF_max_loc 1 --sort=merge --frag=no -F200,1 --secondary=no -a -o /Data/Mapping_${Name}_${Pattern}.sam /Data/all-organisms.fna /FASTQ-files/RefSeq-200GB-HiFi1.fastq

