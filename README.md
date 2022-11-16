# Genome-on-Diet: Taming Large-Scale Genomic Analyses via Sparsified Genomics 

Searching for similar genomic sequences is an essential and fundamental step in biomedical research and an overwhelming majority of genomic analyses. State-of-the-art computational methods performing such comparisons fail to cope with the exponential growth of genomic sequencing data. We are witnessing many tens of petabases (e.g., 17.5x10<sup>15</sup> bases in GenBank alone) of publicly available sequencing data. Global efforts are underway to identify the earth’s virome in preparation for the next pandemic and improve genomic diversity through building population-specific reference genomes. The genomic data that needs to be quickly analyzed will quickly grow with such efforts and ever-more powerful sequencing technologies (e.g., the portable Nanopore technology or the production-scale Illumina technology). 

We introduce the concept of sparsified genomics where we systematically exclude a large number of bases from genomic sequences and enable the processing of the sparsified, shorter genomic sequence while maintaining the same or better accuracy than that of processing non-sparsified sequences. Sparsified genomics provides significant benefits to most genomic analyses and has broad applicability. Sparsifying genomic sequences greatly accelerates state-of-the-art read mappers by 2.98-12.11x using real Illumina, HiFi, and ONT reads, while providing a higher number of mapped reads with the highest mapping quality and more detected indels and complex structural variations. Sparsifying genomic sequences makes containment search through very large genomes and very large databases 72.7-75.88x faster and 723.3x more space-efficient than searching through non-sparsified genomic sequences. Sparsifying genomic sequences enables robust microbiome discovery by providing 54.15-61.88x faster and 720x more space-efficient taxonomic profiling of metagenomic samples.

Described by Alser et al. (preliminary version at https://arxiv.org/abs/........).


## <a name="started"></a>Getting Started
```sh
git clone https://github.com/CMU-SAFARI/Genome-on-Diet
cd Genome-on-Diet/src && make

# Illumina sequences
./GoD_avx -t 40 -ax sr -Z 10 -W 2 -i 2 -k 21 -w 11 -N 1 -r 0.05,120,200 -n 0.9,0 --AF_max_loc 10 --secondary=yes ../Data/GRCh38.p14_NC_000001.11.fna ../Data/ERR240727_1.fastq

# HiFi sequences
./GoD_avx -t 40 -ax map-hifi -Z 10 -W 2 -i 0.2 -k 19 -w 19 -N 1 -r 0.04,400,800 -n 0.8,0.005  --AF_max_loc 10 --sort=merge --secondary=yes ../Data/GRCh38.p14_NC_000001.11.fna ../Data/PBmixSequel729_1_A01_PBTH_30hours_19kbV2PD_70pM_HumanHG003.fastq

# ONT sequences
./GoD_avx -t 40 -ax map-ont -Z 10 -W 2 -i 0.2 -k 15 -w 10 -N 1 -r 0.04,400,800 -n 0.2,0.0005 --AF_max_loc 10 --sort=merge --secondary=yes ../Data/GRCh38.p14_NC_000001.11.fna ../Data/GM24149_1_300filtered_2Mreads.fastq
```

## Table of Contents
- [Getting Started](#started)
- [Key Idea](#idea)
- [Benefits of Genome-on-Diet](#results)
- [Using Genome-on-Diet](#usage)
- [Directory Structure](#directory)
- [Getting help](#contact)
- [Citing Genome-on-Diet](#cite)

##  <a name="idea"></a>The key Idea 
We now need more than ever to catalyze and greatly accelerate genomic analyses by addressing the four critical limitations. Our goal is to enable ultra-fast and highly-efficient indexing and seeding steps in various genomic analyses so that pre-building genome indexes for each genome assembly is no longer a requirement for quickly and directly running large-scale genomic analyses using large genomes and various versions of genome assembly.

We introduce the new concept of sparsified genomics where we systematically exclude a large number of bases from genomic sequences and enable the processing of the sparsified, shorter genomic sequence while maintaining the same or better accuracy than that of processing non-sparsified sequences. We exploit redundancy in genomic sequences to eliminate some regions in the genomic sequences to reduce the input workload of each step of genomic analysis and accelerate overall performance (check te figure below). To demonstrate the benefits of sparsified genomics in real genomic applications, we introduce Genome-on-Diet, the first highly parallel, memory-frugal, and accurate framework for sparsifying genomic sequences.

Genome-on-Diet is based on four key ideas. 
1. Using a repeating pattern sequence to decide on which base in the genomic sequence should be excluded and which one should be included. The pattern sequence is a user-defined shortest repeating substring that represents included and excluded bases by 1’s and 0’s, respectively. Genome-on-Diet provides lossless compression to genomic sequences as it only manipulates a copy of the genomic sequence that is used for the initial steps of an analysis, such as indexing and seeding. The original genomic sequence is still maintained for performing accuracy-critical steps of an analysis, such as base-level sequence alignment where all bases must be accounted for.
2. Inferring at which location of the query sequence the pattern should be applied in order to correctly match the included bases of the query sequence with the included bases of the target sequences. Applying the pattern always starting from the first base of the read sequence can lead to poor results due to a possible lack of seed matches.
3. Providing a highly parallel and highly optimized implementation using modern single-instruction multiple-data (SIMD) instructions employed in modern microprocessors.
4. Introducing four key optimization strategies for making Genome-on-Diet highly parallel, efficient, and accurate. All parameters and optimizations can be conveniently configured, enabled, and disabled using input parameter values entered in the command line of Genome-on-Diet.


![alt text](https://github.com/CMU-SAFARI/Genome-on-Diet/blob/main/Figures/Genome-on-Diet-main.png?raw=true)


##  <a name="results"></a>Benefits of Genome-on-Diet 
Sparsifying genomic sequences makes large-scale analyses feasible and efficient. The use of ‘1110’, ‘110’, ‘10’, and ‘100’ patterns accelerates the analyses we evaluated by 1.3x, 1.4x, 1.9x, and 2.7x, respectively, showing a reduction in the execution time by almost 1/4, 1/3, 1/2, and 2/3, respectively. This demonstrates that the execution time scales linearly with the number of zeros determined in the pattern sequence. The use of ‘1110’, ‘110’, ‘10’, and ‘100’ patterns also directly reduces the size of the index by 1/4, 1/3, 1/2, and 2/3, respectively.

Sparsifying genomic sequences greatly accelerates state-of-the-art read mappers by 2.98-12.11x using real Illumina, HiFi, and ONT reads, while providing a higher number of mapped reads with the highest mapping quality and more detected indels and complex structural variations. Sparsifying genomic sequences makes containment search through very large genomes and very large databases 72.7-75.88x faster and 723.3x more space-efficient than searching through non-sparsified genomic sequences. Sparsifying genomic sequences enables robust microbiome discovery by providing 54.15-61.88x faster and 720x more space-efficient taxonomic profiling of metagenomic samples.

##  <a name="usage"></a>Using Genome-on-Diet:
The concept of sparsified genomics can be implemented and leveraged using different implementations and algorithms.
We introduce one way to implement sparsified genomics using highly-efficient and well-optimized framework, called Genome-on-Diet. The Genome-on-Diet framework is a five-step procedure: compressed indexing, pattern alignment, compressed seeding, location voting, and sequence alignment. These steps can be used individually or collectively together depending on the target genomic analysis.

##  <a name="directory"></a>Directory Structure:
```
Genome-on-Diet-master
├───1. Data
└───2. src
├───3. ReproducibleEvaluation
    └───4. ReadMapping
    └───5. ContainmentSearch
    └───6. TaxonomicProfiling        
├───7. RawResults
    └───8. ReadMapping
    └───9. ContainmentSearch
    └───10. TaxonomicProfiling 
├───11. Figures
```            
1. In the "Data" directory, you will find a description of the sequencing reads and metagenomic reads that we used in our evaluation. You will also find details on how to obtain the reference genomes and RefSeq database that we used in our evaluation. This enables you to use the scripts we have in the "ReproducibleEvaluation" directory to reproduce the exact same experimental evaluations.
2. In the "src" directory, you will find the source code of the Genome-on-Diet implementation. We use the source code of minimap2 as a baseline, where we applied our changes directly to minimap2. This enables researchers that are familiar with the source code of minimap2 to be also familiar with the source code of Genome-on-Diet with minimal efforts.
3. In the "ReproducibleEvaluation" directory, you will find all shell scripts and commands we used to run the experimental evaluations we presented in the paper. You will find three subdirectories, (4), (5), and (6), covering the three applications that we evaluated.
7. In the "RawResults" directory, you will find all raw results and exact values that we presented in the paper for the scripts we provided in (3). You will find three subdirectories, (8), (9), and (10), covering the three applications that we evaluated.
11. In the "Figures" directory, you will find the highh-quality figures that we included in the paper.


##  <a name="contact"></a>Getting Help
If you have any suggestion for improvement, new applications, or collaboration, please contact alserm at ethz dot ch
If you encounter bugs or have further questions or requests, you can raise an issue at the [issue page][issue].

## <a name="cite"></a>Citing Genome-on-Diet

If you use Genome-on-Diet in your work, please cite:

> Mohammed Alser,  Julien Eudine, and Onur Mutlu. 
> "Taming Large-Scale Genomic Analyses via Sparsified Genomics" 
> (2022). [link](https://doi.org/.....)

Below is bibtex format for citation.

```bibtex
@article{Genome-on-Diet2022,
    author = {Alser, Mohammed and Eudine, Julien and Mutlu, Onur},
    title = "{Taming Large-Scale Genomic Analyses via Sparsified Genomics}",
    year = {2020},
}
```


[issue]: https://github.com/CMU-SAFARI/Genome-on-Diet/issues
