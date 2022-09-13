
#############################################
############### Performing variant calling using sniffles
#############################################

conda create -n SNIFFLES sniffles=2.0
conda activate SNIFFLES
conda install -c bioconda samtools
conda install -c conda-forge ncurses
conda install -c bioconda bcftools
conda install -c bioconda rtg-tools

conda install -c bioconda gatk4




cd /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Illumina


sbatch --exclusive -w kratos5 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Genome-on-Diet-GRCh38-Illumina-stats-command-out_k21w11_VC.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Genome-on-Diet-GRCh38-Illumina-stats-time-memory_k21w11_VC.txt"  /home/alserm/minimap2-alser/src2/minimap2_avx -t 40 -ax sr -Z 10 -W 2 -i 2 -k 21 -w 11 -N 1 -r 0.05,120,200 -n 0.9,0 --AF_max_loc 10 --secondary=yes -a -o "/home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Illumina/Genome-on-Diet-GRCh38-Illumina-stats_k21w11_VC.sam" /home/alserm/minimap2-alser/Data/GCF_000001405.40_GRCh38.p14_genomic.fna /home/alserm/Molecules2Variations/ERR240727_1.fastq"

# Convert to bam file
sbatch --exclusive -w kratos5 -c 40 -J bwa_hg38 --wrap="samtools view -b /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Illumina/Genome-on-Diet-GRCh38-Illumina-stats_k21w11_VC.sam  > /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Illumina/Genome-on-Diet-GRCh38-Illumina-stats_k21w11.bam"

# Sort the bam file
sbatch --exclusive -w kratos5 -c 40 -J bwa_hg38 --wrap="samtools sort -o /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Illumina/Genome-on-Diet-GRCh38-Illumina-stats_k21w11_sorted.bam /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Illumina/Genome-on-Diet-GRCh38-Illumina-stats_k21w11.bam"

# create an index file
sbatch --exclusive -w kratos5 -c 40 -J bwa_hg38 --wrap="samtools index /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Illumina/Genome-on-Diet-GRCh38-Illumina-stats_k21w11_sorted.bam"


#gatk MarkDuplicates -I /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Illumina/Genome-on-Diet-GRCh38-Illumina-stats_k21w11_sorted.bam -O marked_duplicates.bam -M marked_dup_metrics.txt

sbatch --exclusive -w kratos5 -c 40 -J bwa_hg38 --wrap="samtools faidx /home/alserm/minimap2-alser/Data/GCF_000001405.40_GRCh38.p14_genomic.fna" 

sbatch --exclusive -w kratos5 -c 40 -J bwa_hg38 --wrap="gatk CreateSequenceDictionary -R /home/alserm/minimap2-alser/Data/GCF_000001405.40_GRCh38.p14_genomic.fna"
 
#gatk ValidateSamFile -I /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Illumina/Genome-on-Diet-GRCh38-Illumina-stats_k21w11_sorted.bam -MODE SUMMARY

sbatch --exclusive -w kratos5 -c 40 -J bwa_hg38 --wrap="gatk AddOrReplaceReadGroups -I /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Illumina/Genome-on-Diet-GRCh38-Illumina-stats_k21w11_sorted.bam -O /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Illumina/Genome-on-Diet-GRCh38-Illumina-stats_k21w11_sorted_fixed.bam --SORT_ORDER coordinate --RGID foo --RGPU unit1 --RGLB bar --RGPL illumina --RGSM Sample1 --CREATE_INDEX True"


sbatch --exclusive -w kratos5 -c 40 -J bwa_hg38 --wrap="gatk HaplotypeCaller -R /home/alserm/minimap2-alser/Data/GCF_000001405.40_GRCh38.p14_genomic.fna -I /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Illumina/Genome-on-Diet-GRCh38-Illumina-stats_k21w11_sorted_fixed.bam -O /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Illumina/variants-GRCh38-Genome-on-Diet-Illumina.vcf"

sbatch --exclusive -w kratos5 -c 40 -J bwa_hg38 --wrap="gatk CountVariants -V /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Illumina/variants-GRCh38-Genome-on-Diet-Illumina.vcf > /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Illumina/variants-GRCh38-Genome-on-Diet-Illumina.txt"





sbatch --exclusive -w kratos2 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/minimap2-GRCh38-Illumina-stats-command-out_k21w11.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/minimap2-GRCh38-Illumina-stats-time-memory_k21w11.txt" /home/alserm/minimap2-alser/minimap2-original-7June2022/minimap2 -t 40 -ax sr -N 1 -a -k 21 -w 11 --secondary=yes -o "/home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Illumina/minimap2-GRCh38-Illumina-stats_k21w11.sam" /home/alserm/minimap2-alser/Data/GCF_000001405.40_GRCh38.p14_genomic.fna /home/alserm/Molecules2Variations/ERR240727_1.fastq"

# Convert to bam file
sbatch --exclusive -w kratos2 -c 40 -J bwa_hg38 --wrap="samtools view -b /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Illumina/minimap2-GRCh38-Illumina-stats_k21w11.sam  > /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Illumina/minimap2-GRCh38-Illumina-stats_k21w11.bam"

# Sort the bam file
sbatch --exclusive -w kratos2 -c 40 -J bwa_hg38 --wrap="samtools sort -o /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Illumina/minimap2-GRCh38-Illumina-stats_k21w11_sorted.bam /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Illumina/minimap2-GRCh38-Illumina-stats_k21w11.bam"

# create an index file
sbatch --exclusive -w kratos2 -c 40 -J bwa_hg38 --wrap="samtools index /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Illumina/minimap2-GRCh38-Illumina-stats_k21w11_sorted.bam"


#gatk MarkDuplicates -I /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Illumina/minimap2-GRCh38-Illumina-stats_k21w11_sorted.bam -O marked_duplicates.bam -M marked_dup_metrics.txt

sbatch --exclusive -w kratos2 -c 40 -J bwa_hg38 --wrap="samtools faidx /home/alserm/minimap2-alser/Data/GCF_000001405.40_GRCh38.p14_genomic.fna" 

sbatch --exclusive -w kratos2 -c 40 -J bwa_hg38 --wrap="gatk CreateSequenceDictionary -R /home/alserm/minimap2-alser/Data/GCF_000001405.40_GRCh38.p14_genomic.fna" 

#gatk ValidateSamFile -I /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Illumina/minimap2-GRCh38-Illumina-stats_k21w11_sorted.bam -MODE SUMMARY

sbatch --exclusive -w kratos2 -c 40 -J bwa_hg38 --wrap="gatk AddOrReplaceReadGroups -I /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Illumina/minimap2-GRCh38-Illumina-stats_k21w11_sorted.bam -O /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Illumina/minimap2-GRCh38-Illumina-stats_k21w11_sorted_fixed.bam --SORT_ORDER coordinate --RGID foo --RGPU unit1 --RGLB bar --RGPL illumina --RGSM Sample1 --CREATE_INDEX True"


sbatch --exclusive -w kratos2 -c 40 -J bwa_hg38 --wrap="gatk HaplotypeCaller -R /home/alserm/minimap2-alser/Data/GCF_000001405.40_GRCh38.p14_genomic.fna -I /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Illumina/minimap2-GRCh38-Illumina-stats_k21w11_sorted_fixed.bam -O /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Illumina/variants-GRCh38-minimap2-Illumina.vcf"


sbatch --exclusive -w kratos2 -c 40 -J bwa_hg38 --wrap="gatk CountVariants -V /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Illumina/variants-GRCh38-minimap2-Illumina.vcf > /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Illumina/variants-GRCh38-minimap2-Illumina.txt"













cd /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/HiFi

sbatch --exclusive -w kratos5 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Genome-on-Diet-GRCh38-HiFi-stats-command-out_k19w19_VC.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Genome-on-Diet-GRCh38-HiFi-stats-time-memory_k19w19_VC.txt"  /home/alserm/minimap2-alser/src2/minimap2_avx -t 40 --MD -ax map-hifi -Z 10 -W 2 -i 0.2 -k 19 -w 19 -N 1 -r 0.04,400,800 -n 0.8,0.005 --AF_max_loc 10 --sort=merge --frag=no -F200,1 --secondary=yes -a -o "/home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/HiFi/Genome-on-Diet-GRCh38-HiFi-stats_k19w19_VC.sam" /home/alserm/minimap2-alser/Data/GCF_000001405.40_GRCh38.p14_genomic.fna /home/alserm/Molecules2Variations/PBmixSequel729_1_A01_PBTH_30hours_19kbV2PD_70pM_HumanHG003.fastq"



# Convert to bam file
sbatch --exclusive -w kratos5 -c 40 -J bwa_hg38 --wrap="samtools view --threads 40 -bS /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/HiFi/Genome-on-Diet-GRCh38-HiFi-stats_k19w19_VC.sam  > /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/HiFi/Genome-on-Diet-GRCh38-HiFi-stats_k19w19.bam"

# Sort the bam file
sbatch --exclusive -w kratos5 -c 40 -J bwa_hg38 --wrap="samtools sort -o /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/HiFi/Genome-on-Diet-GRCh38-HiFi-stats_k19w19_sorted.bam /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/HiFi/Genome-on-Diet-GRCh38-HiFi-stats_k19w19.bam"

# create an index file
sbatch --exclusive -w kratos5 -c 40 -J bwa_hg38 --wrap="samtools index /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/HiFi/Genome-on-Diet-GRCh38-HiFi-stats_k19w19_sorted.bam"


sbatch --exclusive -w kratos5 -c 40 -J bwa_hg38 --wrap="sniffles --input /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/HiFi/Genome-on-Diet-GRCh38-HiFi-stats_k19w19_sorted.bam --allow-overwrite --vcf /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/HiFi/variants-GRCh38-Genome-on-Diet-HiFi.vcf --threads 40"




sbatch --exclusive -w kratos1 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/minimap2-GRCh38-HiFi-stats-command-out_k19w19_VC.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/minimap2-GRCh38-HiFi-stats-time-memory_k19w19_VC.txt" /home/alserm/minimap2-alser/minimap2-original-7June2022/minimap2 -t 40 --MD -ax map-hifi -N 1 -a -k 19 -w 19 --secondary=yes -o "/home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/HiFi/minimap2-GRCh38-HiFi-stats_k19w19_VC.sam" /home/alserm/minimap2-alser/Data/GCF_000001405.40_GRCh38.p14_genomic.fna /home/alserm/Molecules2Variations/PBmixSequel729_1_A01_PBTH_30hours_19kbV2PD_70pM_HumanHG003.fastq"



# Convert to bam file
sbatch --exclusive -w kratos1 -c 40 -J bwa_hg38 --wrap="samtools view --threads 40 -bS /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/HiFi/minimap2-GRCh38-HiFi-stats_k19w19_VC.sam  > /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/HiFi/minimap2-GRCh38-HiFi-stats_k19w19.bam"

# Sort the bam file
sbatch --exclusive -w kratos1 -c 40 -J bwa_hg38 --wrap="samtools sort -o /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/HiFi/minimap2-GRCh38-HiFi-stats_k19w19_sorted.bam /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/HiFi/minimap2-GRCh38-HiFi-stats_k19w19.bam"

# create an index file
sbatch --exclusive -w kratos1 -c 40 -J bwa_hg38 --wrap="samtools index /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/HiFi/minimap2-GRCh38-HiFi-stats_k19w19_sorted.bam"


sbatch --exclusive -w kratos1 -c 40 -J bwa_hg38 --wrap="sniffles --input /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/HiFi/minimap2-GRCh38-HiFi-stats_k19w19_sorted.bam --allow-overwrite --vcf /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/HiFi/variants-GRCh38-minimap2-HiFi.vcf --threads 40"



sbatch --exclusive -w kratos1 -c 40 -J bwa_hg38 --wrap="rtg vcfstats /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/HiFi/variants-GRCh38-minimap2-HiFi.vcf /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/HiFi/variants-GRCh38-Genome-on-Diet-HiFi.vcf"













cd /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/ONT

sbatch --exclusive -w kratos5 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Genome-on-Diet-GRCh38-ONT-stats-command-out_k15w10_VC.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Genome-on-Diet-GRCh38-ONT-stats-time-memory_k15w10_VC.txt" /home/alserm/minimap2-alser/src2/minimap2_avx -t 40 -ax map-ont --MD -Z 10 -W 2 -i 0.2 -k 15 -w 10 -N 1 -r 0.04,400,800 -n 0.2,0.0005 --AF_max_loc 10 --sort=merge --frag=no -F200,1 --secondary=yes -a -o "/home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/ONT/Genome-on-Diet-GRCh38-ONT-stats_k15w10_VC.sam" /home/alserm/minimap2-alser/Data/GCF_000001405.40_GRCh38.p14_genomic.fna /home/alserm/Molecules2Variations/GM24149_1_300filtered_2Mreads.fastq"


cd /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/ONT
# Convert to bam file
sbatch --exclusive -w kratos5 -c 40 -J bwa_hg38 --wrap="samtools view --threads 40 -bS /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/ONT/Genome-on-Diet-GRCh38-ONT-stats_k15w10_VC.sam  > /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/ONT/Genome-on-Diet-GRCh38-ONT-stats_k15w10_VC.bam"

# Sort the bam file
sbatch --exclusive -w kratos5 -c 40 -J bwa_hg38 --wrap="samtools sort -o /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/ONT/Genome-on-Diet-GRCh38-ONT-stats_k15w10_sorted_VC.bam /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/ONT/Genome-on-Diet-GRCh38-ONT-stats_k15w10_VC.bam"

# create an index file
sbatch --exclusive -w kratos5 -c 40 -J bwa_hg38 --wrap="samtools index /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/ONT/Genome-on-Diet-GRCh38-ONT-stats_k15w10_sorted_VC.bam"


sbatch --exclusive -w kratos5 -c 40 -J bwa_hg38 --wrap="sniffles --input /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/ONT/Genome-on-Diet-GRCh38-ONT-stats_k15w10_sorted_VC.bam --allow-overwrite --vcf /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/ONT/variants-GRCh38-Genome-on-Diet-ONT.vcf --threads 40"




sbatch --exclusive -w kratos2 -c 40 -J bwa_hg38 -o "/home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/minimap2-GRCh38-ONT-stats-command-out_k15w10_VC.txt" --wrap="/usr/bin/time -v --output "/home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/minimap2-GRCh38-ONT-stats-time-memory_k15w10_VC.txt" /home/alserm/minimap2-alser/minimap2-original-7June2022/minimap2 -t 40 -ax map-ont --MD -N 1 -a -k 15 -w 10 -o "/home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/ONT/minimap2-GRCh38-ONT-stats_k15w10_VC.sam" /home/alserm/minimap2-alser/Data/GCF_000001405.40_GRCh38.p14_genomic.fna /home/alserm/Molecules2Variations/GM24149_1_300filtered_2Mreads.fastq"



cd /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/ONT

# Convert to bam file
sbatch --exclusive -w kratos2 -c 40 -J bwa_hg38 --wrap="samtools view --threads 40 -bS /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/ONT/minimap2-GRCh38-ONT-stats_k15w10_VC.sam  > /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/ONT/minimap2-GRCh38-ONT-stats_k15w10_VC.bam"

# Sort the bam file
sbatch --exclusive -w kratos2 -c 40 -J bwa_hg38 --wrap="samtools sort -o /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/ONT/minimap2-GRCh38-ONT-stats_k15w10_sorted_VC.bam /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/ONT/minimap2-GRCh38-ONT-stats_k15w10_VC.bam"

# create an index file
sbatch --exclusive -w kratos2 -c 40 -J bwa_hg38 --wrap="samtools index /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/ONT/minimap2-GRCh38-ONT-stats_k15w10_sorted_VC.bam"


sbatch --exclusive -w kratos2 -c 40 -J bwa_hg38 --wrap="sniffles --input /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/ONT/minimap2-GRCh38-ONT-stats_k15w10_sorted_VC.bam --allow-overwrite --vcf /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/ONT/variants-GRCh38-minimap2-ONT.vcf --threads 40"




sbatch --exclusive -w kratos2 -c 40 -J bwa_hg38 --wrap="rtg vcfstats /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/ONT/variants-GRCh38-minimap2-ONT.vcf /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/ONT/variants-GRCh38-Genome-on-Diet-ONT.vcf" 













############# Find common variants
#get the ground-truth variants from annotated genome: ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/
#https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg19&c=chrX&g=giab

#echo "/home/alserm/minimap2-alser/Use-cases/Read-mapping/HG003_GRCh38_1_22_v4.2.1_all.vcf" > lsvar.txt
# Add header of VCF and Filter by DP https://gatk.broadinstitute.org/hc/en-us/articles/360035531692-VCF-Variant-Call-Format
awk -F':' '(substr($1,1,1)=="#" || $7>=3){print $0}' /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Illumina/variants-GRCh38-Genome-on-Diet-Illumina.vcf > /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Illumina/variants-GRCh38-Genome-on-Diet-Illumina_filtered.vcf

awk -F':' '(substr($1,1,1)=="#" || $7>=3){print $0}' /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Illumina/variants-GRCh38-minimap2-Illumina.vcf > /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Illumina/variants-GRCh38-minimap2-Illumina_filtered.vcf

bgzip --force /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Illumina/variants-GRCh38-minimap2-Illumina_filtered.vcf
tabix --force -p vcf /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Illumina/variants-GRCh38-minimap2-Illumina_filtered.vcf.gz

bgzip --force /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Illumina/variants-GRCh38-Genome-on-Diet-Illumina_filtered.vcf
tabix --force -p vcf /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Illumina/variants-GRCh38-Genome-on-Diet-Illumina_filtered.vcf.gz

vcf-compare /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Illumina/variants-GRCh38-minimap2-Illumina_filtered.vcf.gz /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Illumina/variants-GRCh38-Genome-on-Diet-Illumina_filtered.vcf.gz



#bgzip -d /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Illumina/variants-GRCh38-minimap2-Illumina_filtered.vcf.gz
#bgzip -d /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Illumina/variants-GRCh38-Genome-on-Diet-Illumina_filtered.vcf.gz
#
#
#ls /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/Illumina/variants-GRCh38*_filtered.vcf > lsvar.txt
#
#SURVIVOR merge lsvar.txt 1000 1 1 0 0 1 /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/sample_merged.vcf 
#
#SURVIVOR genComp /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/sample_merged.vcf 0 /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/sample_merged.mat.txt
#
#perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/' /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/sample_merged.vcf | sed -e 's/\(.\)/\1 /g' | grep '1 1' |wc -l
#perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/' /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/sample_merged.vcf | sed -e 's/\(.\)/\1 /g' | grep '1 0' |wc -l
#perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/' /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/sample_merged.vcf | sed -e 's/\(.\)/\1 /g' | grep '0 1' |wc -l
#
#
#
#SURVIVOR merge lsvar.txt 1000 1 1 0 0 51 /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/sample_merged.vcf 
#
#SURVIVOR genComp /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/sample_merged.vcf 0 /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/sample_merged.mat.txt
#
#perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/' /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/sample_merged.vcf | sed -e 's/\(.\)/\1 /g' | grep '1 1' |wc -l
#perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/' /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/sample_merged.vcf | sed -e 's/\(.\)/\1 /g' | grep '1 0' |wc -l
#perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/' /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/sample_merged.vcf | sed -e 's/\(.\)/\1 /g' | grep '0 1' |wc -l






# Add header of VCF and Filter by Support number of reads
awk -F'\t' '(substr($1,1,1)=="#"){print $0}' /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/HiFi/variants-GRCh38-minimap2-HiFi.vcf > /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/HiFi/variants-GRCh38-minimap2-HiFi_filtered.vcf

grep --invert-match "SUPPORT=[0-2];" /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/HiFi/variants-GRCh38-minimap2-HiFi.vcf >> /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/HiFi/variants-GRCh38-minimap2-HiFi_filtered.vcf


awk -F'\t' '(substr($1,1,1)=="#"){print $0}' /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/HiFi/variants-GRCh38-Genome-on-Diet-HiFi.vcf > /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/HiFi/variants-GRCh38-Genome-on-Diet-HiFi_filtered.vcf

grep --invert-match "SUPPORT=[0-2];" /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/HiFi/variants-GRCh38-Genome-on-Diet-HiFi.vcf >> /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/HiFi/variants-GRCh38-Genome-on-Diet-HiFi_filtered.vcf


ls /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/HiFi/variants-GRCh38*_filtered.vcf > lsvar.txt

SURVIVOR merge lsvar.txt 1000 1 1 0 0 1 /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/sample_merged.vcf 

SURVIVOR genComp /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/sample_merged.vcf 0 /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/sample_merged.mat.txt

perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/' /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/sample_merged.vcf | sed -e 's/\(.\)/\1 /g' | grep '1 1' | wc -l > HiFi-variant.txt
perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/' /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/sample_merged.vcf | sed -e 's/\(.\)/\1 /g' | grep '1 0' | wc -l >> HiFi-variant.txt
perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/' /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/sample_merged.vcf | sed -e 's/\(.\)/\1 /g' | grep '0 1' | wc -l >> HiFi-variant.txt



SURVIVOR merge lsvar.txt 1000 1 1 0 0 51 /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/sample_merged.vcf 

SURVIVOR genComp /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/sample_merged.vcf 0 /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/sample_merged.mat.txt

perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/' /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/sample_merged.vcf | sed -e 's/\(.\)/\1 /g' | grep '1 1' | wc -l >> HiFi-variant.txt
perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/' /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/sample_merged.vcf | sed -e 's/\(.\)/\1 /g' | grep '1 0' | wc -l >> HiFi-variant.txt
perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/' /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/sample_merged.vcf | sed -e 's/\(.\)/\1 /g' | grep '0 1' | wc -l >> HiFi-variant.txt















# Add header of VCF and Filter by Support number of reads
awk -F'\t' '(substr($1,1,1)=="#"){print $0}' /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/ONT/variants-GRCh38-minimap2-ONT.vcf > /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/ONT/variants-GRCh38-minimap2-ONT_filtered.vcf

grep --invert-match "SUPPORT=[0-2];" /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/ONT/variants-GRCh38-minimap2-ONT.vcf >> /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/ONT/variants-GRCh38-minimap2-ONT_filtered.vcf


awk -F'\t' '(substr($1,1,1)=="#"){print $0}' /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/ONT/variants-GRCh38-Genome-on-Diet-ONT.vcf > /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/ONT/variants-GRCh38-Genome-on-Diet-ONT_filtered.vcf

grep --invert-match "SUPPORT=[0-2];" /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/ONT/variants-GRCh38-Genome-on-Diet-ONT.vcf >> /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/ONT/variants-GRCh38-Genome-on-Diet-ONT_filtered.vcf


ls /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/ONT/variants-GRCh38*_filtered.vcf > lsvar.txt

SURVIVOR merge lsvar.txt 1000 1 1 0 0 1 /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/sample_merged.vcf 

SURVIVOR genComp /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/sample_merged.vcf 0 /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/sample_merged.mat.txt

perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/' /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/sample_merged.vcf | sed -e 's/\(.\)/\1 /g' | grep '1 1' | wc -l > ONT-variant.txt
perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/' /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/sample_merged.vcf | sed -e 's/\(.\)/\1 /g' | grep '1 0' | wc -l >> ONT-variant.txt
perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/' /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/sample_merged.vcf | sed -e 's/\(.\)/\1 /g' | grep '0 1' | wc -l >> ONT-variant.txt



SURVIVOR merge lsvar.txt 1000 1 1 0 0 51 /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/sample_merged.vcf 

SURVIVOR genComp /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/sample_merged.vcf 0 /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/sample_merged.mat.txt

perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/' /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/sample_merged.vcf | sed -e 's/\(.\)/\1 /g' | grep '1 1' | wc -l >> ONT-variant.txt
perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/' /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/sample_merged.vcf | sed -e 's/\(.\)/\1 /g' | grep '1 0' | wc -l >> ONT-variant.txt
perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/' /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/sample_merged.vcf | sed -e 's/\(.\)/\1 /g' | grep '0 1' | wc -l >> ONT-variant.txt













# https://www.melbournebioinformatics.org.au/tutorials/tutorials/longread_sv_calling/longread_sv_calling/

#awk '/SVTYPE=/ {split($8,infoArr,";"); print substr(infoArr[9],8), $3, $1, substr(infoArr[3],6), $2, substr(infoArr[4],5), substr(infoArr[11],7)}' /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/ONT/variants-GRCh38-minimap2-ONT.vcf


#awk '/SVTYPE=/ {split($8,infoArr,";"); print substr(infoArr[9],8), $3, $1, substr(infoArr[3],6), $2, substr(infoArr[4],5), substr(infoArr[11],7)}' /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/ONT/variants-GRCh38-Genome-on-Diet-ONT.vcf


