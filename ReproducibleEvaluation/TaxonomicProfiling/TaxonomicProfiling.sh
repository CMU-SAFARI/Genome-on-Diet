

#sed 's./..g' /CAMI-datasets/CAMI_high/RH_S001__insert_270.fq > /CAMI-datasets/CAMI_high/RH_S001__insert_270_SingleEnd.fq

#sed 's./..g' /CAMI-datasets/CAMI_low/RL_S001__insert_270.fq > /CAMI-datasets/CAMI_low/RL_S001__insert_270_SingleEnd.fq


#################################################       
####### Original Metalign
## CAMI High k21w11
python3 /Metalign/scripts/metalign.py /CAMI-datasets/CAMI_high/RH_S001__insert_270_SingleEnd.fq /Metalign/ERR_test/ --threads 40 --precise --keep_temp_files --minimap2 /Metalign/minimap2-2.17/minimap2 --temp_dir /Metalign/tmp-Metalign-CAMI_High/ --output /Metalign/tmp-Metalign-CAMI_High/metalign_CAMI_High.tsv

## CAMI Low k21w11
python3 /Metalign/scripts/metalign.py /CAMI-datasets/CAMI_low/RL_S001__insert_270_SingleEnd.fq /Metalign/ERR_test/ --threads 40 --precise --keep_temp_files --minimap2 /Metalign/minimap2-2.17/minimap2 --temp_dir /Metalign/tmp-Metalign-CAMI_Low/ --output /Metalign/tmp-Metalign-CAMI_Low/metalign_CAMI_Low.tsv



cp -r /Metalign/tmp-Metalign-CAMI_High/* /Metalign/tmp-Metalign-CAMI_High_k28w40
cp -r /Metalign/tmp-Metalign-CAMI_Low/* /Metalign/tmp-Metalign-CAMI_Low_k28w40
## CAMI High K28w40
python3 /Metalign-k28w40/scripts/metalign.py /CAMI-datasets/CAMI_high/RH_S001__insert_270_SingleEnd.fq /Metalign/ERR_test/ --threads 40 --precise --keep_temp_files --minimap2 /Metalign/minimap2-2.17/minimap2 --temp_dir /Metalign/tmp-Metalign-CAMI_High_k28w40/ --output /Metalign/tmp-Metalign-CAMI_High_k28w40/metalign_CAMI_High.tsv

## CAMI Low K28w40
python3 /Metalign-k28w40/scripts/metalign.py /CAMI-datasets/CAMI_low/RL_S001__insert_270_SingleEnd.fq /Metalign/ERR_test/ --threads 40 --precise --keep_temp_files --minimap2 /Metalign/minimap2-2.17/minimap2 --temp_dir /Metalign/tmp-Metalign-CAMI_Low_k28w40/ --output /Metalign/tmp-Metalign-CAMI_Low_k28w40/metalign_CAMI_Low.tsv





#################################################       
####### map_and_profile_only Genome-on-Diet
cp -r /Metalign/tmp-Metalign-CAMI_High/* /Metalign/tmp-GoD-CAMI_High_k21w11
cp -r /Metalign/tmp-Metalign-CAMI_Low/* /Metalign/tmp-GoD-CAMI_Low_k21w11
## CAMI High k21w11
python3 /GoD-map_n_profile_only_k21w11_High/scripts/metalign.py /CAMI-datasets/CAMI_high/RH_S001__insert_270_SingleEnd.fq /Metalign/ERR_test/ --threads 40 --precise --keep_temp_files --minimap2 /GoD-map_n_profile_only_k21w11/GDiet_avx --temp_dir /Metalign/tmp-GoD-CAMI_High_k21w11/ --output /Metalign/tmp-GoD-CAMI_High_k21w11/GoD_CAMI_High.tsv

## CAMI Low k21w11
python3 /GoD-map_n_profile_only_k21w11/scripts/metalign.py /CAMI-datasets/CAMI_low/RL_S001__insert_270_SingleEnd.fq /Metalign/ERR_test/ --threads 40 --precise --keep_temp_files --minimap2 /GoD-map_n_profile_only_k21w11/GDiet_avx --temp_dir /Metalign/tmp-GoD-CAMI_Low_k21w11/ --output /Metalign/tmp-GoD-CAMI_Low_k21w11/GoD_CAMI_Low.tsv



cp -r /Metalign/tmp-Metalign-CAMI_High/* /Metalign/tmp-GoD-CAMI_High_k28w40
cp -r /Metalign/tmp-Metalign-CAMI_Low/* /Metalign/tmp-GoD-CAMI_Low_k28w40
## CAMI High k28w40
python3 /GoD-map_n_profile_only_k28w40/scripts/metalign.py /CAMI-datasets/CAMI_high/RH_S001__insert_270_SingleEnd.fq /Metalign/ERR_test/ --threads 40 --precise --keep_temp_files --minimap2 /GoD-map_n_profile_only_k28w40/GDiet_avx --temp_dir /Metalign/tmp-GoD-CAMI_High_k28w40/ --output /Metalign/tmp-GoD-CAMI_High_k28w40/GoD_CAMI_High.tsv

## CAMI Low k28w40
python3 /GoD-map_n_profile_only_k28w40/scripts/metalign.py /CAMI-datasets/CAMI_low/RL_S001__insert_270_SingleEnd.fq /Metalign/ERR_test/ --threads 40 --precise --keep_temp_files --minimap2 /GoD-map_n_profile_only_k28w40/GDiet_avx --temp_dir /Metalign/tmp-GoD-CAMI_Low_k28w40/ --output /Metalign/tmp-GoD-CAMI_Low_k28w40/GoD_CAMI_Low.tsv





#################################################       
####### Genome-on-Diet + Genome-on-Diet
cp -r /Metalign/tmp-Metalign-CAMI_High/* /Metalign/tmp-GoD-GoD-CAMI_High_k21w11
cp -r /Metalign/tmp-Metalign-CAMI_Low/* /Metalign/tmp-GoD-GoD-CAMI_Low_k21w11
cp /Metalign/ERR_test/organism_files_combined_results_CAMI_Low_k28w40/cmashed_db.fna /Metalign/tmp-GoD-GoD-CAMI_Low_k21w11/
grep 'NZ_PPGA0' /Metalign/ERR_test/db_info.txt >> /Metalign/tmp-GoD-GoD-CAMI_Low_k21w11/subset_db_info.txt

## CAMI High k21w11
python3 /GoD-map_n_profile_only_k21w11_High/scripts/metalign.py /CAMI-datasets/CAMI_high/RH_S001__insert_270_SingleEnd.fq /Metalign/ERR_test/ --threads 40 --precise --keep_temp_files --minimap2 /GoD-map_n_profile_only_k21w11/GDiet_avx --temp_dir /Metalign/tmp-GoD-GoD-CAMI_High_k21w11/ --output /Metalign/tmp-GoD-GoD-CAMI_High_k21w11/GoD_CAMI_High.tsv

## CAMI Low k21w11
python3 /GoD-map_n_profile_only_k21w11/scripts/metalign.py /CAMI-datasets/CAMI_low/RL_S001__insert_270_SingleEnd.fq /Metalign/ERR_test/ --threads 40 --precise --keep_temp_files --minimap2 /GoD-map_n_profile_only_k21w11/GDiet_avx --temp_dir /Metalign/tmp-GoD-GoD-CAMI_Low_k21w11/ --output /Metalign/tmp-GoD-GoD-CAMI_Low_k21w11/GoD_CAMI_Low.tsv


#################################################       
####### Genome-on-Diet + minimap2
cp -r /Metalign/tmp-Metalign-CAMI_High/* /Metalign/tmp-GoD-minimap2-CAMI_High_k21w11
cp -r /Metalign/tmp-Metalign-CAMI_Low/* /Metalign/tmp-GoD-minimap2-CAMI_Low_k21w11
cp /Metalign/ERR_test/organism_files_combined_results_CAMI_Low_k28w40/cmashed_db.fna /Metalign/tmp-GoD-minimap2-CAMI_Low_k21w11/
grep 'NZ_PPGA0' /Metalign/ERR_test/db_info.txt >> /Metalign/tmp-GoD-minimap2-CAMI_Low_k21w11/subset_db_info.txt

## CAMI High k21w11
python3 /Metalign-map-and-profile/scripts/metalign.py /CAMI-datasets/CAMI_high/RH_S001__insert_270_SingleEnd.fq /Metalign/ERR_test/ --threads 40 --precise --keep_temp_files --minimap2 /Metalign/minimap2-2.17/minimap2 --temp_dir /Metalign/tmp-GoD-minimap2-CAMI_High_k21w11/ --output /Metalign/tmp-GoD-minimap2-CAMI_High_k21w11/GoD_CAMI_High.tsv

## CAMI Low k21w11
python3 /Metalign-map-and-profile/scripts/metalign.py /CAMI-datasets/CAMI_low/RL_S001__insert_270_SingleEnd.fq /Metalign/ERR_test/ --threads 40 --precise --keep_temp_files --minimap2 /Metalign/minimap2-2.17/minimap2 --temp_dir /Metalign/tmp-GoD-minimap2-CAMI_Low_k21w11/ --output /Metalign/tmp-GoD-minimap2-CAMI_Low_k21w11/GoD_CAMI_Low.tsv

#################################################       
####### OPAL

conda activate OPAL  

cd /home/alserm/OPAL-master

# KMC+CMash as containment index
# High
python3 /home/alserm/OPAL-master/opal.py -g /Metalign/tmp-Metalign-CAMI_High/metalign_CAMI_High.tsv /Metalign/tmp-Metalign-CAMI_High/metalign_CAMI_High.tsv /Metalign/tmp-GoD-CAMI_High_k21w11/GoD_CAMI_High.tsv -l "Metalign-k21w11, Genome-on-Diet-k21w11" -o /home/alserm/OPAL-master/CAMI_High_k21w11
  
# Low
python3 /home/alserm/OPAL-master/opal.py -g /Metalign/tmp-Metalign-CAMI_Low/metalign_CAMI_Low.tsv /Metalign/tmp-Metalign-CAMI_Low/metalign_CAMI_Low.tsv /Metalign/tmp-GoD-CAMI_Low_k21w11/GoD_CAMI_Low.tsv -l "Metalign-k21w11, Genome-on-Diet-k21w11" -o /home/alserm/OPAL-master/CAMI_Low_k21w11


# Genome-on-Diet as containment index
# High
python3 /home/alserm/OPAL-master/opal.py -g /Metalign/tmp-GoD-minimap2-CAMI_High_k21w11/GoD_CAMI_High.tsv /Metalign/tmp-GoD-GoD-CAMI_High_k21w11/GoD_CAMI_High.tsv /Metalign/tmp-GoD-minimap2-CAMI_High_k21w11/GoD_CAMI_High.tsv -l "GoD-GoD-k21w11, GoD-miniam2-k21w11" -o /home/alserm/OPAL-master/CAMI_High-GoD_k21w11
# Low
python3 /home/alserm/OPAL-master/opal.py -g /Metalign/tmp-GoD-minimap2-CAMI_Low_k21w11/GoD_CAMI_Low.tsv /Metalign/tmp-GoD-GoD-CAMI_Low_k21w11/GoD_CAMI_Low.tsv /Metalign/tmp-GoD-minimap2-CAMI_Low_k21w11/GoD_CAMI_Low.tsv -l "GoD-GoD-k21w11, GoD-miniam2-k21w11" -o /home/alserm/OPAL-master/CAMI_Low-GoD_k21w11

