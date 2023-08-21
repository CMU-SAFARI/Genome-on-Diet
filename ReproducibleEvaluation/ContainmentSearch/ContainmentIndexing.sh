#!/bin/bash
set -e

#################################################
#Human Genome
Name='Human'
./GDiet_avx -t 40 -x map-hifi -I 300G -k 19 -w 16 --idx-no-seq -d Indexing_GDiet_${Name}.mmi /Data/GCF_000001405.40_GRCh38.p14_genomic.fna

Pattern='11'
PatternLength='2'
Name='Human'
./GDiet_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 --idx-no-seq -d /Data/Indexing_${Name}_${Pattern}.mmi /Data/GCF_000001405.40_GRCh38.p14_genomic.fna

Pattern='10'
PatternLength='2'
Name='Human'
./GDiet_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 --idx-no-seq -d /Data/Indexing_${Name}_${Pattern}.mmi /Data/GCF_000001405.40_GRCh38.p14_genomic.fna


Pattern='100'
PatternLength='3'
Name='Human'
./GDiet_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 --idx-no-seq -d /Data/Indexing_${Name}_${Pattern}.mmi /Data/GCF_000001405.40_GRCh38.p14_genomic.fna


Pattern='110'
PatternLength='3'
Name='Human'
./GDiet_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 --idx-no-seq -d /Data/Indexing_${Name}_${Pattern}.mmi /Data/GCF_000001405.40_GRCh38.p14_genomic.fna


Pattern='1110'
PatternLength='4'
Name='Human'
./GDiet_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 --idx-no-seq -d /Data/Indexing_${Name}_${Pattern}.mmi /Data/GCF_000001405.40_GRCh38.p14_genomic.fna


#################################################
#Largest genome
cd /Data/
Name='Pinus_taeda'
./GDiet_avx -t 40 -x map-hifi -I 300G -k 19 -w 16 --idx-no-seq -d /Data/Indexing_GDiet_${Name}.mmi /Data/Pinus_taeda_loblolly-pine/GCA_001447015.2_Sugar_pine_JHU_assembly_genomic.fna


Pattern='11'
PatternLength='2'
Name='Pinus_taeda'
./GDiet_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 --idx-no-seq -d /Data/Indexing_${Name}_${Pattern}.mmi /Data/Pinus_taeda_loblolly-pine/GCA_001447015.2_Sugar_pine_JHU_assembly_genomic.fna

Pattern='10'
PatternLength='2'
Name='Pinus_taeda'
./GDiet_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 --idx-no-seq -d /Data/Indexing_${Name}_${Pattern}.mmi /Data/Pinus_taeda_loblolly-pine/GCA_001447015.2_Sugar_pine_JHU_assembly_genomic.fna


Pattern='100'
PatternLength='3'
Name='Pinus_taeda'
./GDiet_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 --idx-no-seq -d /Data/Indexing_${Name}_${Pattern}.mmi /Data/Pinus_taeda_loblolly-pine/GCA_001447015.2_Sugar_pine_JHU_assembly_genomic.fna


Pattern='110'
PatternLength='3'
Name='Pinus_taeda'
./GDiet_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 --idx-no-seq -d /Data/Indexing_${Name}_${Pattern}.mmi /Data/Pinus_taeda_loblolly-pine/GCA_001447015.2_Sugar_pine_JHU_assembly_genomic.fna


Pattern='1110'
PatternLength='4'
Name='Pinus_taeda'
./GDiet_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 --idx-no-seq -d /Data/Indexing_${Name}_${Pattern}.mmi /Data/Pinus_taeda_loblolly-pine/GCA_001447015.2_Sugar_pine_JHU_assembly_genomic.fna


#################################################
#RefSeq-100GB
cd /Data/

Name='RefSeq-100GB'
./GDiet_avx -t 40 -x map-hifi -I 300G -k 19 -w 16 --idx-no-seq -d /Data/Indexing_GDiet_${Name}.mmi /Data/all-organisms.fna

Pattern='11'
PatternLength='2'
Name='RefSeq-100GB'
./GDiet_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 --idx-no-seq -d /Data/Indexing_${Name}_${Pattern}.mmi /Data/all-organisms.fna

Pattern='10'
PatternLength='2'
Name='RefSeq-100GB'
./GDiet_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 --idx-no-seq -d /Data/Indexing_${Name}_${Pattern}.mmi /Data/all-organisms.fna


Pattern='100'
PatternLength='3'
Name='RefSeq-100GB'
./GDiet_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 --idx-no-seq -d /Data/Indexing_${Name}_${Pattern}.mmi /Data/all-organisms.fna



Pattern='110'
PatternLength='3'
Name='RefSeq-100GB'
./GDiet_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 --idx-no-seq -d /Data/Indexing_${Name}_${Pattern}.mmi /Data/all-organisms.fna


Pattern='1110'
PatternLength='4'
Name='RefSeq-100GB'
./GDiet_avx -t 40 -x map-hifi -I 300G -Z ${Pattern} -W ${PatternLength} -k 19 -w 16 --idx-no-seq -d /Data/Indexing_${Name}_${Pattern}.mmi /Data/all-organisms.fna

