

#!/bin/bash
set -e

# 180, 360
#cd /home/alserm/Metalign/Metalign/ERR_test/organism_files_unzip/
#
#files=(*)
#
#for (( i = 0; i < ${#files[@]}; i += 360 ))
#do
#    cat -- "${files[@]:$i:360}" > "../organism_files_combined/XLargeCombinedFile_$i.fna"
#   
#done
#    
#cp /home/alserm/Metalign/Metalign/ERR_test/organism_files_unzip/taxid_2088_genomic.fna /home/alserm/Metalign/Metalign/ERR_test/organism_files_combined/
#cp /home/alserm/Metalign/Metalign/ERR_test/organism_files_unzip/taxid_2371_1_genomic.fna /home/alserm/Metalign/Metalign/ERR_test/organism_files_combined/
#cp /home/alserm/Metalign/Metalign/ERR_test/organism_files_unzip/taxid_328515_1_genomic.fna /home/alserm/Metalign/Metalign/ERR_test/organism_files_combined/
#cp /home/alserm/Metalign/Metalign/ERR_test/organism_files_unzip/taxid_266_3_genomic.fna /home/alserm/Metalign/Metalign/ERR_test/organism_files_combined/
#cp /home/alserm/Metalign/Metalign/ERR_test/organism_files_unzip/taxid_990712_1_genomic.fna /home/alserm/Metalign/Metalign/ERR_test/organism_files_combined/
#cp /home/alserm/Metalign/Metalign/ERR_test/organism_files_unzip/taxid_1004304_genomic.fna /home/alserm/Metalign/Metalign/ERR_test/organism_files_combined/
#cp /home/alserm/Metalign/Metalign/ERR_test/organism_files_unzip/taxid_173053_genomic.fna /home/alserm/Metalign/Metalign/ERR_test/organism_files_combined/
#cp /home/alserm/Metalign/Metalign/ERR_test/organism_files_unzip/taxid_404881_genomic.fna /home/alserm/Metalign/Metalign/ERR_test/organism_files_combined/
#cp /home/alserm/Metalign/Metalign/ERR_test/organism_files_unzip/taxid_266_1_genomic.fna /home/alserm/Metalign/Metalign/ERR_test/organism_files_combined/
#cp /home/alserm/Metalign/Metalign/ERR_test/organism_files_unzip/taxid_266_2_genomic.fna /home/alserm/Metalign/Metalign/ERR_test/organism_files_combined/
#cp /home/alserm/Metalign/Metalign/ERR_test/organism_files_unzip/taxid_1123015_1_genomic.fna /home/alserm/Metalign/Metalign/ERR_test/organism_files_combined/
#cp /home/alserm/Metalign/Metalign/ERR_test/organism_files_unzip/taxid_644383_genomic.fna /home/alserm/Metalign/Metalign/ERR_test/organism_files_combined/
# 
# 
#cd /home/alserm/Metalign/Metalign/ERR_test/
#
#for file in `find organism_files_combined/ -type f -name "*.fna"`
#do
#   # sbatch --wrap="gzip -k -f "../organism_files_combined/CombinedFile_$i.fna""
#   
#done

> /home/alserm/Metalign/Metalign/ERR_test/organism_files_combined_results_CAMI_Low_k28w40/cmashed_db.fna
> /home/alserm/Metalign/Metalign/ERR_test/organism_files_combined_results_CAMI_High_k28w40/cmashed_db.fna

cd /home/alserm/Metalign/Metalign/ERR_test/
i=0
for file in `find organism_files_combined/ -type f -name "XLargeCombinedFile_*.fna"` 
do
	
	## k=28, w =40
    ##### add flag --sam-hit-only to get the exact number of mapped reads per strain/species.
	sbatch --exclusive -w kratos5 -o "/home/alserm/Metalign/Metalign/ERR_test/organism_files_combined_results_CAMI_Low_k28w40/command-GoD-RL_S001_insert_270_hits$i.txt" -c 40 -J bwa_hg38 --wrap="/usr/bin/time -v --output "/home/alserm/Metalign/Metalign/ERR_test/organism_files_combined_results_CAMI_Low_k28w40/time-memory-GoD-RL_S001_insert_270_hits$i.txt" /home/alserm/minimap2-alser/Genome-on-Diet-voting-KMC-Aug2022/src/minimap2_avx -t 40 -ax sr -Z 10 -W 2 -i 2 -k 28 -w 48 -N 100000 -r 0.05,120,200 -n 0.9,0 --AF_max_loc 100000 --secondary=yes -p 100000 -A 1 -I 30G ${file} /home/alserm/minimap2-alser/CAMI-datasets/CAMI_low/RL_S001__insert_270_SingleEnd.fq >> /home/alserm/Metalign/Metalign/ERR_test/organism_files_combined_results_CAMI_Low_k28w40/cmashed_db.fna"
    
    ##### add flag --sam-hit-only to get the exact number of mapped reads per strain/species.
    sbatch --exclusive -w kratos5 -o "/home/alserm/Metalign/Metalign/ERR_test/organism_files_combined_results_CAMI_High_k28w40/command-GoD-RH_S001__insert_270_hits$i.txt" -c 40 -J bwa_hg38 --wrap="/usr/bin/time -v --output "/home/alserm/Metalign/Metalign/ERR_test/organism_files_combined_results_CAMI_High_k28w40/time-memory-GoD-RH_S001__insert_270_hits$i.txt" /home/alserm/minimap2-alser/Genome-on-Diet-voting-KMC-Aug2022/src/minimap2_avx -t 40 -ax sr -Z 10 -W 2 -i 2 -k 28 -w 40 -N 100000 -r 0.05,120,200 -n 0.9,0 --AF_max_loc 100000 --secondary=yes -p 100000 -A 1 -I 30G ${file} /home/alserm/minimap2-alser/CAMI-datasets/CAMI_high/RH_S001__insert_270_SingleEnd.fq >> /home/alserm/Metalign/Metalign/ERR_test/organism_files_combined_results_CAMI_High_k28w40/cmashed_db.fna"
    
    
    #sbatch --exclusive -w kratos5 -o "/home/alserm/Metalign/Metalign/ERR_test/organism_files_combined_results_CAMI_High_k28w40/command-GoD-RH_S001__insert_270_hits_allgenomes.txt" -c 40 -J bwa_hg38 --wrap="/usr/bin/time -v --output "/home/alserm/Metalign/Metalign/ERR_test/organism_files_combined_results_CAMI_High_k28w40/time-memory-GoD-RH_S001__insert_270_hits_allgenomes.txt" /home/alserm/minimap2-alser/Genome-on-Diet-voting-KMC-Aug2022/src/minimap2_avx -t 40 -ax sr -Z 10 -W 2 -i 2 -k 28 -w 40 -N 100000 -r 0.05,120,200 -n 0.9,0 --AF_max_loc 100000 --secondary=yes -p 100000 -A 1 /home/alserm/minimap2-alser/Data/all-organisms.fna /home/alserm/minimap2-alser/CAMI-datasets/CAMI_high/RH_S001__insert_270_SingleEnd.fq >> /home/alserm/Metalign/Metalign/ERR_test/organism_files_combined_results_CAMI_High_k28w40/cmashed_db_allgenomes.fna"
    
    #sbatch --exclusive -w kratos4 -o "/home/alserm/Metalign/Metalign/ERR_test/organism_files_combined_results_CAMI_Low_k28w40/command-GoD-RL_S001_insert_270_hits_allgenomes.txt" -c 40 -J bwa_hg38 --wrap="/usr/bin/time -v --output "/home/alserm/Metalign/Metalign/ERR_test/organism_files_combined_results_CAMI_Low_k28w40/time-memory-GoD-RL_S001_insert_270_hits_allgenomes.txt" /home/alserm/minimap2-alser/Genome-on-Diet-voting-KMC-Aug2022/src/minimap2_avx -t 40 -ax sr -Z 10 -W 2 -i 2 -k 28 -w 40 -N 100000 -r 0.05,120,200 -n 0.9,0 --AF_max_loc 100000 --secondary=yes -p 100000 -A 1 /home/alserm/minimap2-alser/Data/all-organisms.fna /home/alserm/minimap2-alser/CAMI-datasets/CAMI_low/RL_S001__insert_270_SingleEnd.fq >> /home/alserm/Metalign/Metalign/ERR_test/organism_files_combined_results_CAMI_Low_k28w40/cmashed_db_allgenomes.fna"
    
    
    ## k=21, w =11
	##### add flag --sam-hit-only to get the exact number of mapped reads per strain/species.
	#sbatch --exclusive -w kratos3 -o "/home/alserm/Metalign/Metalign/ERR_test/organism_files_combined_results_CAMI_Low_k21w11/command-GoD-RL_S001_insert_270_hits$i.txt" -c 40 -J bwa_hg38 --wrap="/usr/bin/time -v --output "/home/alserm/Metalign/Metalign/ERR_test/organism_files_combined_results_CAMI_Low_k21w11/time-memory-GoD-RL_S001_insert_270_hits$i.txt" /home/alserm/minimap2-alser/Genome-on-Diet-voting-KMC-Aug2022/src/minimap2_avx -t 40 -ax sr -Z 10 -W 2 -i 2 -k 21 -w 11 -N 100000 -r 0.05,120,200 -n 0.9,0 --AF_max_loc 100000 --secondary=yes -p 100000 -A 10 ${file} /home/alserm/minimap2-alser/CAMI-datasets/CAMI_low/RL_S001__insert_270_SingleEnd.fq >> /home/alserm/Metalign/Metalign/ERR_test/organism_files_combined_results_CAMI_Low_k21w11/cmashed_db.fna"
    
    ##### add flag --sam-hit-only to get the exact number of mapped reads per strain/species.
    #sbatch --exclusive -w kratos4 -o "/home/alserm/Metalign/Metalign/ERR_test/organism_files_combined_results_CAMI_High_k21w11/command-GoD-RH_S001__insert_270_hits$i.txt" -c 40 -J bwa_hg38 --wrap="/usr/bin/time -v --output "/home/alserm/Metalign/Metalign/ERR_test/organism_files_combined_results_CAMI_High_k21w11/time-memory-GoD-RH_S001__insert_270_hits$i.txt" /home/alserm/minimap2-alser/Genome-on-Diet-voting-KMC-Aug2022/src/minimap2_avx -t 40 -ax sr -Z 10 -W 2 -i 2 -k 21 -w 11 -N 100000 -r 0.05,120,200 -n 0.9,0 --AF_max_loc 100000 --secondary=yes -p 100000 -A 10 ${file} /home/alserm/minimap2-alser/CAMI-datasets/CAMI_high/RH_S001__insert_270_SingleEnd.fq >> /home/alserm/Metalign/Metalign/ERR_test/organism_files_combined_results_CAMI_High_k21w11/cmashed_db.fna"
    
    ((i=i+1))
done

grep '>' /home/alserm/Metalign/Metalign/tmp-Metalign-CAMI_Low/cmashed_db.fna | awk -F' ' '{print substr($1,1,7)}' | uniq -c

grep '>' /home/alserm/Metalign/Metalign/ERR_test/organism_files_combined_results_CAMI_Low_k28w40/cmashed_db.fna | awk -F' ' '{print substr($1,1,7)}' | uniq -c



grep '>' /home/alserm/Metalign/Metalign/tmp-Metalign-CAMI_High/cmashed_db.fna | awk -F' ' '{print substr($1,1,7)}' | uniq -c

grep '>' /home/alserm/Metalign/Metalign/ERR_test/organism_files_combined_results_CAMI_High_k28w40/cmashed_db.fna | awk -F' ' '{print substr($1,1,7)}' | uniq -c



grep 'seconds' /home/alserm/Metalign/Metalign/ERR_test/organism_files_combined_results_CAMI_High_k28w40/* |awk -F' ' '{sum+=$5;} END{print sum;}'

grep 'seconds' /home/alserm/Metalign/Metalign/ERR_test/organism_files_combined_results_CAMI_Low_k28w40/* |awk -F' ' '{sum+=$5;} END{print sum;}'

grep 'kbytes' /home/alserm/Metalign/Metalign/ERR_test/organism_files_combined_results_CAMI_High_k28w40/* |awk -F' ' 'max<$7 || NR==1{ max=$7; } END{ print max }'

grep 'kbytes' /home/alserm/Metalign/Metalign/ERR_test/organism_files_combined_results_CAMI_Low_k28w40/* |awk -F' ' 'max<$7 || NR==1{ max=$7; } END{ print max }'





# /home/alserm/minimap2-alser/Genome-on-Diet-LatestMinimap2-KMC/src/minimap2 -t 40 -ax sr -Z 10 -W 2 -i 2 -k 28 -w 28 -N 6 --AF_dis 1.2 -n 0.9 --AF_max_loc 5 --secondary=yes -p 100 /home/alserm/Metalign/Metalign/ERR_test/organism_files_combined/taxid_143453_genomic.fna /home/alserm/minimap2-alser/CAMI-datasets/CAMI_low/RL_S001__insert_270.fq



#/home/alserm/minimap2-alser/Genome-on-Diet-LatestMinimap2-KMC/src/minimap2_avx -t 40 -ax sr -Z 10 -W 2 -i 2 -k 21 -w 40 -N 6 --AF_dis 1.2 -n 0.9 --AF_max_loc 5 --secondary=yes -p 100 /home/alserm/Metalign/Metalign/ERR_test/organism_files_combined/taxid_143453_genomic.fna /home/alserm/minimap2-alser/CAMI-datasets/CAMI_high/RH_S001__insert_270.fq

#
#cat /home/alserm/Metalign/Metalign/ERR_test/organism_files_combined_results_CAMI_Low/cmashed_db_Genome-on-Diet_RL_S001_* > tst3.txt
#grep '>' /home/alserm/Metalign/Metalign/ERR_test/organism_files_combined_results_CAMI_Low_k28w40/cmashed_db.fna | cut -c 2- > tst3.txt
#
#sort -k 1,1 tst3.txt > sorted_tst3.txt
#
#sort -k 1,1 /home/alserm/Metalign/Metalign/tmp-CAMI_Low_RL_S001__insert_270/subset_db_info.txt > sorted_subset_db_info.txt 
#
#join sorted_tst3.txt sorted_subset_db_info.txt | awk -F' ' '{print $3}' | uniq > tst4.txt
#
#
#
#cat /home/alserm/Metalign/Metalign/ERR_test/organism_files_combined_results_CAMI_High/cmashed_db_Genome-on-Diet_RL_S001_* > tst3.txt
#
#grep '>' /home/alserm/Metalign/Metalign/ERR_test/organism_files_combined_results_CAMI_High_k28w40/cmashed_db.fna | cut -c 2- > tst3.txt
#
#sort -k 1,1 tst3.txt > sorted_tst3.txt
#
#sort -k 1,1 /home/alserm/Metalign/Metalign/tmp-CAMI_High_RH_S001__insert_270/subset_db_info.txt > sorted_subset_db_info.txt 
#
#join sorted_tst3.txt sorted_subset_db_info.txt | awk -F' ' '{print $3}' | uniq > tst4.txt
##



