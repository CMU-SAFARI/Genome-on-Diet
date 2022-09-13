

sed 's./..g' /home/alserm/minimap2-alser/CAMI-datasets/CAMI_high/RH_S001__insert_270.fq > /home/alserm/minimap2-alser/CAMI-datasets/CAMI_high/RH_S001__insert_270_SingleEnd.fq

sed 's./..g' /home/alserm/minimap2-alser/CAMI-datasets/CAMI_low/RL_S001__insert_270.fq > /home/alserm/minimap2-alser/CAMI-datasets/CAMI_low/RL_S001__insert_270_SingleEnd.fq


#################################################       
####### Original Metalign
## CAMI High k21w11
sbatch --exclusive -w kratos5 -o /home/alserm/Metalign/Metalign/tmp-Metalign-CAMI_High/sbatch-metalign.txt -c 40 -J bwa_hg38 --wrap="/usr/bin/time -v --output /home/alserm/Metalign/Metalign/tmp-Metalign-CAMI_High/time-memory-metalign.txt python3 /home/alserm/minimap2-alser/Use-cases/Metagenomics/metalign/scripts/metalign.py /home/alserm/minimap2-alser/CAMI-datasets/CAMI_high/RH_S001__insert_270_SingleEnd.fq /home/alserm/Metalign/Metalign/ERR_test/ --threads 40 --precise --keep_temp_files --minimap2 /home/alserm/minimap2-alser/Use-cases/Metagenomics/metalign/minimap2-2.17/minimap2 --temp_dir /home/alserm/Metalign/Metalign/tmp-Metalign-CAMI_High/ --output /home/alserm/Metalign/Metalign/tmp-Metalign-CAMI_High/metalign_CAMI_High.tsv"
## CAMI Low k21w11
sbatch --exclusive -w kratos5 -o /home/alserm/Metalign/Metalign/tmp-Metalign-CAMI_Low/sbatch-metalign.txt -c 40 -J bwa_hg38 --wrap="/usr/bin/time -v --output /home/alserm/Metalign/Metalign/tmp-Metalign-CAMI_Low/time-memory-metalign.txt python3 /home/alserm/minimap2-alser/Use-cases/Metagenomics/metalign/scripts/metalign.py /home/alserm/minimap2-alser/CAMI-datasets/CAMI_low/RL_S001__insert_270_SingleEnd.fq /home/alserm/Metalign/Metalign/ERR_test/ --threads 40 --precise --keep_temp_files --minimap2 /home/alserm/minimap2-alser/Use-cases/Metagenomics/metalign/minimap2-2.17/minimap2 --temp_dir /home/alserm/Metalign/Metalign/tmp-Metalign-CAMI_Low/ --output /home/alserm/Metalign/Metalign/tmp-Metalign-CAMI_Low/metalign_CAMI_Low.tsv"



cp -r /home/alserm/Metalign/Metalign/tmp-Metalign-CAMI_High/* /home/alserm/Metalign/Metalign/tmp-Metalign-CAMI_High_k28w40
cp -r /home/alserm/Metalign/Metalign/tmp-Metalign-CAMI_Low/* /home/alserm/Metalign/Metalign/tmp-Metalign-CAMI_Low_k28w40
## CAMI High K28w40
sbatch --exclusive -w kratos5 -o /home/alserm/Metalign/Metalign/tmp-Metalign-CAMI_High_k28w40/sbatch-metalign.txt -c 40 -J bwa_hg38 --wrap="/usr/bin/time -v --output /home/alserm/Metalign/Metalign/tmp-Metalign-CAMI_High_k28w40/time-memory-metalign.txt python3 /home/alserm/minimap2-alser/Use-cases/Metagenomics/metalign-k28w40/scripts/metalign.py /home/alserm/minimap2-alser/CAMI-datasets/CAMI_high/RH_S001__insert_270_SingleEnd.fq /home/alserm/Metalign/Metalign/ERR_test/ --threads 40 --precise --keep_temp_files --minimap2 /home/alserm/minimap2-alser/Use-cases/Metagenomics/metalign/minimap2-2.17/minimap2 --temp_dir /home/alserm/Metalign/Metalign/tmp-Metalign-CAMI_High_k28w40/ --output /home/alserm/Metalign/Metalign/tmp-Metalign-CAMI_High_k28w40/metalign_CAMI_High.tsv"
## CAMI Low K28w40
sbatch --exclusive -w kratos5 -o /home/alserm/Metalign/Metalign/tmp-Metalign-CAMI_Low_k28w40/sbatch-metalign.txt -c 40 -J bwa_hg38 --wrap="/usr/bin/time -v --output /home/alserm/Metalign/Metalign/tmp-Metalign-CAMI_Low_k28w40/time-memory-metalign.txt python3 /home/alserm/minimap2-alser/Use-cases/Metagenomics/metalign-k28w40/scripts/metalign.py /home/alserm/minimap2-alser/CAMI-datasets/CAMI_low/RL_S001__insert_270_SingleEnd.fq /home/alserm/Metalign/Metalign/ERR_test/ --threads 40 --precise --keep_temp_files --minimap2 /home/alserm/minimap2-alser/Use-cases/Metagenomics/metalign/minimap2-2.17/minimap2 --temp_dir /home/alserm/Metalign/Metalign/tmp-Metalign-CAMI_Low_k28w40/ --output /home/alserm/Metalign/Metalign/tmp-Metalign-CAMI_Low_k28w40/metalign_CAMI_Low.tsv"





#################################################       
####### map_and_profile_only Genome-on-Diet
cp -r /home/alserm/Metalign/Metalign/tmp-Metalign-CAMI_High/* /home/alserm/Metalign/Metalign/tmp-GoD-CAMI_High_k21w11
cp -r /home/alserm/Metalign/Metalign/tmp-Metalign-CAMI_Low/* /home/alserm/Metalign/Metalign/tmp-GoD-CAMI_Low_k21w11
## CAMI High k21w11
sbatch --exclusive -w kratos5 -o /home/alserm/Metalign/Metalign/tmp-GoD-CAMI_High_k21w11/sbatch-GoD.txt -c 40 -J bwa_hg38 --wrap="/usr/bin/time -v --output /home/alserm/Metalign/Metalign/tmp-GoD-CAMI_High_k21w11/time-memory-GoD.txt python3 /home/alserm/minimap2-alser/Use-cases/Metagenomics/GoD-map_n_profile_only_k21w11_High/scripts/metalign.py /home/alserm/minimap2-alser/CAMI-datasets/CAMI_high/RH_S001__insert_270_SingleEnd.fq /home/alserm/Metalign/Metalign/ERR_test/ --threads 40 --precise --keep_temp_files --minimap2 /home/alserm/minimap2-alser/Use-cases/Metagenomics/GoD-map_n_profile_only_k21w11/Genome-on-Diet-1Sept2022/minimap2_avx --temp_dir /home/alserm/Metalign/Metalign/tmp-GoD-CAMI_High_k21w11/ --output /home/alserm/Metalign/Metalign/tmp-GoD-CAMI_High_k21w11/GoD_CAMI_High.tsv"
## CAMI Low k21w11
sbatch --exclusive -w kratos5 -o /home/alserm/Metalign/Metalign/tmp-GoD-CAMI_Low_k21w11/sbatch-GoD.txt -c 40 -J bwa_hg38 --wrap="/usr/bin/time -v --output /home/alserm/Metalign/Metalign/tmp-GoD-CAMI_Low_k21w11/time-memory-GoD.txt python3 /home/alserm/minimap2-alser/Use-cases/Metagenomics/GoD-map_n_profile_only_k21w11/scripts/metalign.py /home/alserm/minimap2-alser/CAMI-datasets/CAMI_low/RL_S001__insert_270_SingleEnd.fq /home/alserm/Metalign/Metalign/ERR_test/ --threads 40 --precise --keep_temp_files --minimap2 /home/alserm/minimap2-alser/Use-cases/Metagenomics/GoD-map_n_profile_only_k21w11/Genome-on-Diet-1Sept2022/minimap2_avx --temp_dir /home/alserm/Metalign/Metalign/tmp-GoD-CAMI_Low_k21w11/ --output /home/alserm/Metalign/Metalign/tmp-GoD-CAMI_Low_k21w11/GoD_CAMI_Low.tsv"



cp -r /home/alserm/Metalign/Metalign/tmp-Metalign-CAMI_High/* /home/alserm/Metalign/Metalign/tmp-GoD-CAMI_High_k28w40
cp -r /home/alserm/Metalign/Metalign/tmp-Metalign-CAMI_Low/* /home/alserm/Metalign/Metalign/tmp-GoD-CAMI_Low_k28w40
## CAMI High k28w40
sbatch --exclusive -w kratos5 -o /home/alserm/Metalign/Metalign/tmp-GoD-CAMI_High_k28w40/sbatch-GoD.txt -c 40 -J bwa_hg38 --wrap="/usr/bin/time -v --output /home/alserm/Metalign/Metalign/tmp-GoD-CAMI_High_k28w40/time-memory-GoD.txt python3 /home/alserm/minimap2-alser/Use-cases/Metagenomics/GoD-map_n_profile_only_k28w40/scripts/metalign.py /home/alserm/minimap2-alser/CAMI-datasets/CAMI_high/RH_S001__insert_270_SingleEnd.fq /home/alserm/Metalign/Metalign/ERR_test/ --threads 40 --precise --keep_temp_files --minimap2 /home/alserm/minimap2-alser/Use-cases/Metagenomics/GoD-map_n_profile_only_k28w40/Genome-on-Diet-1Sept2022/minimap2_avx --temp_dir /home/alserm/Metalign/Metalign/tmp-GoD-CAMI_High_k28w40/ --output /home/alserm/Metalign/Metalign/tmp-GoD-CAMI_High_k28w40/GoD_CAMI_High.tsv"
## CAMI Low k28w40
sbatch --exclusive -w kratos5 -o /home/alserm/Metalign/Metalign/tmp-GoD-CAMI_Low_k28w40/sbatch-GoD.txt -c 40 -J bwa_hg38 --wrap="/usr/bin/time -v --output /home/alserm/Metalign/Metalign/tmp-GoD-CAMI_Low_k28w40/time-memory-GoD.txt python3 /home/alserm/minimap2-alser/Use-cases/Metagenomics/GoD-map_n_profile_only_k28w40/scripts/metalign.py /home/alserm/minimap2-alser/CAMI-datasets/CAMI_low/RL_S001__insert_270_SingleEnd.fq /home/alserm/Metalign/Metalign/ERR_test/ --threads 40 --precise --keep_temp_files --minimap2 /home/alserm/minimap2-alser/Use-cases/Metagenomics/GoD-map_n_profile_only_k28w40/Genome-on-Diet-1Sept2022/minimap2_avx --temp_dir /home/alserm/Metalign/Metalign/tmp-GoD-CAMI_Low_k28w40/ --output /home/alserm/Metalign/Metalign/tmp-GoD-CAMI_Low_k28w40/GoD_CAMI_Low.tsv"





#################################################       
####### Genome-on-Diet + Genome-on-Diet
cp -r /home/alserm/Metalign/Metalign/tmp-Metalign-CAMI_High/* /home/alserm/Metalign/Metalign/tmp-GoD-GoD-CAMI_High_k21w11
cp -r /home/alserm/Metalign/Metalign/tmp-Metalign-CAMI_Low/* /home/alserm/Metalign/Metalign/tmp-GoD-GoD-CAMI_Low_k21w11
cp /home/alserm/Metalign/Metalign/ERR_test/organism_files_combined_results_CAMI_Low_k28w40/cmashed_db.fna /home/alserm/Metalign/Metalign/tmp-GoD-GoD-CAMI_Low_k21w11/
grep 'NZ_PPGA0' /home/alserm/Metalign/Metalign/ERR_test/db_info.txt >> /home/alserm/Metalign/Metalign/tmp-GoD-GoD-CAMI_Low_k21w11/subset_db_info.txt

## CAMI High k21w11
sbatch --exclusive -w kratos5 -o /home/alserm/Metalign/Metalign/tmp-GoD-GoD-CAMI_High_k21w11/sbatch-GoD.txt -c 40 -J bwa_hg38 --wrap="/usr/bin/time -v --output /home/alserm/Metalign/Metalign/tmp-GoD-GoD-CAMI_High_k21w11/time-memory-GoD.txt python3 /home/alserm/minimap2-alser/Use-cases/Metagenomics/GoD-map_n_profile_only_k21w11_High/scripts/metalign.py /home/alserm/minimap2-alser/CAMI-datasets/CAMI_high/RH_S001__insert_270_SingleEnd.fq /home/alserm/Metalign/Metalign/ERR_test/ --threads 40 --precise --keep_temp_files --minimap2 /home/alserm/minimap2-alser/Use-cases/Metagenomics/GoD-map_n_profile_only_k21w11/Genome-on-Diet-1Sept2022/minimap2_avx --temp_dir /home/alserm/Metalign/Metalign/tmp-GoD-GoD-CAMI_High_k21w11/ --output /home/alserm/Metalign/Metalign/tmp-GoD-GoD-CAMI_High_k21w11/GoD_CAMI_High.tsv"
## CAMI Low k21w11
sbatch --exclusive -w kratos5 -o /home/alserm/Metalign/Metalign/tmp-GoD-GoD-CAMI_Low_k21w11/sbatch-GoD.txt -c 40 -J bwa_hg38 --wrap="/usr/bin/time -v --output /home/alserm/Metalign/Metalign/tmp-GoD-GoD-CAMI_Low_k21w11/time-memory-GoD.txt python3 /home/alserm/minimap2-alser/Use-cases/Metagenomics/GoD-map_n_profile_only_k21w11/scripts/metalign.py /home/alserm/minimap2-alser/CAMI-datasets/CAMI_low/RL_S001__insert_270_SingleEnd.fq /home/alserm/Metalign/Metalign/ERR_test/ --threads 40 --precise --keep_temp_files --minimap2 /home/alserm/minimap2-alser/Use-cases/Metagenomics/GoD-map_n_profile_only_k21w11/Genome-on-Diet-1Sept2022/minimap2_avx --temp_dir /home/alserm/Metalign/Metalign/tmp-GoD-GoD-CAMI_Low_k21w11/ --output /home/alserm/Metalign/Metalign/tmp-GoD-GoD-CAMI_Low_k21w11/GoD_CAMI_Low.tsv"


#################################################       
####### Genome-on-Diet + minimap2
cp -r /home/alserm/Metalign/Metalign/tmp-Metalign-CAMI_High/* /home/alserm/Metalign/Metalign/tmp-GoD-minimap2-CAMI_High_k21w11
cp -r /home/alserm/Metalign/Metalign/tmp-Metalign-CAMI_Low/* /home/alserm/Metalign/Metalign/tmp-GoD-minimap2-CAMI_Low_k21w11
cp /home/alserm/Metalign/Metalign/ERR_test/organism_files_combined_results_CAMI_Low_k28w40/cmashed_db.fna /home/alserm/Metalign/Metalign/tmp-GoD-minimap2-CAMI_Low_k21w11/
grep 'NZ_PPGA0' /home/alserm/Metalign/Metalign/ERR_test/db_info.txt >> /home/alserm/Metalign/Metalign/tmp-GoD-minimap2-CAMI_Low_k21w11/subset_db_info.txt

## CAMI High k21w11
sbatch --exclusive -w kratos5 -o /home/alserm/Metalign/Metalign/tmp-GoD-minimap2-CAMI_High_k21w11/sbatch-GoD.txt -c 40 -J bwa_hg38 --wrap="/usr/bin/time -v --output /home/alserm/Metalign/Metalign/tmp-GoD-minimap2-CAMI_High_k21w11/time-memory-GoD.txt python3 /home/alserm/minimap2-alser/Use-cases/Metagenomics/metalign-map-and-profile/scripts/metalign.py /home/alserm/minimap2-alser/CAMI-datasets/CAMI_high/RH_S001__insert_270_SingleEnd.fq /home/alserm/Metalign/Metalign/ERR_test/ --threads 40 --precise --keep_temp_files --minimap2 /home/alserm/minimap2-alser/Use-cases/Metagenomics/metalign/minimap2-2.17/minimap2 --temp_dir /home/alserm/Metalign/Metalign/tmp-GoD-minimap2-CAMI_High_k21w11/ --output /home/alserm/Metalign/Metalign/tmp-GoD-minimap2-CAMI_High_k21w11/GoD_CAMI_High.tsv"
## CAMI Low k21w11
sbatch --exclusive -w kratos5 -o /home/alserm/Metalign/Metalign/tmp-GoD-minimap2-CAMI_Low_k21w11/sbatch-GoD.txt -c 40 -J bwa_hg38 --wrap="/usr/bin/time -v --output /home/alserm/Metalign/Metalign/tmp-GoD-minimap2-CAMI_Low_k21w11/time-memory-GoD.txt python3 /home/alserm/minimap2-alser/Use-cases/Metagenomics/metalign-map-and-profile/scripts/metalign.py /home/alserm/minimap2-alser/CAMI-datasets/CAMI_low/RL_S001__insert_270_SingleEnd.fq /home/alserm/Metalign/Metalign/ERR_test/ --threads 40 --precise --keep_temp_files --minimap2 /home/alserm/minimap2-alser/Use-cases/Metagenomics/metalign/minimap2-2.17/minimap2 --temp_dir /home/alserm/Metalign/Metalign/tmp-GoD-minimap2-CAMI_Low_k21w11/ --output /home/alserm/Metalign/Metalign/tmp-GoD-minimap2-CAMI_Low_k21w11/GoD_CAMI_Low.tsv"

#################################################       
####### OPAL

conda activate OPAL  

cd /home/alserm/OPAL-master

# KMC+CMash as containment index
# High
python3 /home/alserm/OPAL-master/opal.py -g /home/alserm/Metalign/Metalign/tmp-Metalign-CAMI_High/metalign_CAMI_High.tsv /home/alserm/Metalign/Metalign/tmp-Metalign-CAMI_High/metalign_CAMI_High.tsv /home/alserm/Metalign/Metalign/tmp-GoD-CAMI_High_k21w11/GoD_CAMI_High.tsv -l "Metalign-k21w11, Genome-on-Diet-k21w11" -o /home/alserm/OPAL-master/CAMI_High_k21w11
  
# Low
python3 /home/alserm/OPAL-master/opal.py -g /home/alserm/Metalign/Metalign/tmp-Metalign-CAMI_Low/metalign_CAMI_Low.tsv /home/alserm/Metalign/Metalign/tmp-Metalign-CAMI_Low/metalign_CAMI_Low.tsv /home/alserm/Metalign/Metalign/tmp-GoD-CAMI_Low_k21w11/GoD_CAMI_Low.tsv -l "Metalign-k21w11, Genome-on-Diet-k21w11" -o /home/alserm/OPAL-master/CAMI_Low_k21w11


# Genome-on-Diet as containment index
# High
python3 /home/alserm/OPAL-master/opal.py -g /home/alserm/Metalign/Metalign/tmp-GoD-minimap2-CAMI_High_k21w11/GoD_CAMI_High.tsv /home/alserm/Metalign/Metalign/tmp-GoD-GoD-CAMI_High_k21w11/GoD_CAMI_High.tsv /home/alserm/Metalign/Metalign/tmp-GoD-minimap2-CAMI_High_k21w11/GoD_CAMI_High.tsv -l "GoD-GoD-k21w11, GoD-miniam2-k21w11" -o /home/alserm/OPAL-master/CAMI_High-GoD_k21w11
# Low
python3 /home/alserm/OPAL-master/opal.py -g /home/alserm/Metalign/Metalign/tmp-GoD-minimap2-CAMI_Low_k21w11/GoD_CAMI_Low.tsv /home/alserm/Metalign/Metalign/tmp-GoD-GoD-CAMI_Low_k21w11/GoD_CAMI_Low.tsv /home/alserm/Metalign/Metalign/tmp-GoD-minimap2-CAMI_Low_k21w11/GoD_CAMI_Low.tsv -l "GoD-GoD-k21w11, GoD-miniam2-k21w11" -o /home/alserm/OPAL-master/CAMI_Low-GoD_k21w11

