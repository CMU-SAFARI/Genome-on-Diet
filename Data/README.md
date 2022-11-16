# Sequencing Reads:
The short reads (Illumina), the accurate long reads (HiFi), and the ultra-long reads (ONT) are obtained from the NIST's Genome-in-a-Bottle (GIAB)  project: 
```
https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/
```

## PacBio HiFi reads:
```
https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb_20kb_chemistry2/reads/m64011_190830_220126.fastq.gz
```


## ONT ultra-long reads:
We consider only the first 2 million reads whose length is greater than or equal 1000 bp (using NanoFilt --length 1000)
```
https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/Ultralong_OxfordNanopore/guppy-V3.4.5/HG002_ONT-UL_GIAB_20200204.fastq.gz
```

## Illumina 250bp reads:
```
https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/reads/D1_S1_L001_R1_001.fastq.gz
 
https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/reads/D1_S1_L001_R1_002.fastq.gz

https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/reads/D1_S1_L001_R1_003.fastq.gz
```
---

# CAMI Metagenomic Reads:
https://data.cami-challenge.org/
## CAMI Low Complexity
```
RL_S001__insert_270.fq
```
## CAMI High Complexity
```
RH_S001__insert_270.fq
```

---

# Reference Genomes:

## The complete Human Genome GRCh38 was used for variant calling

```
https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz 
```

## The complete Human Genome GRCh38.p14 (GCF_000001405.40, release date 3 February 2022) 

```
https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.40
```

## The largest sequenced reference genome, pinus taeda (also known as loblolly pine, GCA_000404065.3, release date 9 January 2017)

```
https://www.ncbi.nlm.nih.gov/assembly/GCA_000404065.3/
```

---
 
# Metagenomes:
We obtain our metagenomes for building the reference database from RefSeq (https://www.ncbi.nlm.nih.gov/refseq) database.
The list of TAXIDs for the metagenomes we choose is listed as follows:
## RefSeq1
```
RefSeq1_taxids.txt
```
## RefSeq2
```
RefSeq2_taxids.txt
```


