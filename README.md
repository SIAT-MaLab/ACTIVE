# ACTIVE--Activity Comparison Test for Induced Viral Elements
A pipeline for detecting active prophage or other viral regions in bacterial genomes using sequencing depth analysis. It tests whether predicted prophage regions exhibit significantly elevated sequencing read depth compared to the host genomic background — a hallmark of active phage replication or lytic induction.

## Overview
Prophage activity detection is based on a core rationale: Active prophage replication increases their genome copy number relative to their hosts, which manifests as a discrepancy in sequencing depth distributions between the prophage and its host background. ACTIVE quantifies this difference using robust statistical tests and effect size estimation.
![image](https://github.com/SIAT-MaLab/ACTIVE/blob/main/figures1a.png)

The pipeline proceeds through four main stages:

Host CDS prediction — Prodigal predicts coding sequences on all host genomes (fast mode only).

Marker gene extraction — FetchMGs identifies universal single-copy marker genes (MGS) as a proxy for host background depth (fast mode only).

Read mapping — Illumina paired-end reads are aligned with BBMap; Oxford Nanopore reads (if provided) are aligned with Minimap2 and merged into a single BAM.

Statistical analysis — Per-region depth distributions are compared between the prophage region and the host background using normality-aware hypothesis tests (t-test / Welch's t-test / Brunner-Munzel test), with effect sizes (Cohen's D or Cliff's Δ) reported for significant hits.


## Pipeline Modes
| Mode | Background reference | Suitable for |
|------|------|------|
| fast (default) | Single-copy marker gene (MGS) regions identified by FetchMGs | Most use cases, especailly metagenomic generated MAGs |
| precise (`-m`) | All non-prophage genomic positions | High-confidence analysis; slower in statstical test |

In **fast** mode, host background depth is sampled from the genomic coordinates of universal single-copy marker genes, which are constitutively expressed and expected to reflect baseline chromosomal copy number. This provides a biologically meaningful and computationally efficient reference.

In **precise** mode, the entire non-prophage genome is used as background, offering a more conservative and comprehensive comparison at the cost of substantially increased computation time.

## Dependencies
The following tools must be available in environment
| Tool | Purpose |
|------|------|
| Prodigal | CDS prediction (fast mode) |
| FetchMGs=2.1.0 | Single-copy marker gene extraction (fast mode) |
| BBMap | Short-read alignment |
| Minimap2 | Long-read alignment (optional, ONT) |
| samtools=1.21 | BAM processing and depth extraction |
| bedtools | Genomic interval operations |
| Python3 | MGS coordinate extraction (stdlib only; no extra packages required) |
| R=4.3.3 | Statistical analysis |

**Required R packages:**
```r
install.packages(c("tseries", "brunnermunzel", "effsize", "doParallel", "foreach"))
```

## Installation
No installation is required. Download the single script and make it executable:
```bash
bash active.sh -h
```
## Input
### Host genomes (`-g`)

A **directory** containing one or more host genome FASTA files. Each file must have the `.fasta` extension. The filename stem (without extension) is used as the genome identifier throughout the pipeline.
```
genomes/
    strainA.fasta
    strainB.fasta
    strainC.fasta
```
> **Contig naming convention:** The pipeline infers genome identity from the first underscore-delimited field of each contig name. Contigs must follow the pattern <genome_id>_<number> (e.g., EcoliK12_1, EcoliK12_2). **Avoid underscores in the genome identifier itself**, as the pipeline will split on the first _ and treat only the preceding text as the genome prefix. This convention is especially important when multiple genomes are concatenated into a single FASTA reference.

> **Concatenated reference genomes:** Multiple genomes may be concatenated into a single `.fasta` reference file(especially in MAGs from same sample), provided all contigs follow the naming convention above. An **index table (`-i`)** must then be used to direct sample reads to the correct reference file. Note that concatenating closely related genomes can introduce ambiguous multi-mapping reads, which may affect depth estimates even with `ambiguous=random` (BBMap).

### Sequencing data (`-s`)
A **directory** containing sequencing reads. Naming conventions must be followed exactly:
| Data type | Required filename pattern |
|------|------|
| Illumina R1 | `<sample>_1_clean.fastq` or `<sample>_1_clean.fastq.gz` |
| Illumina R2 | `<sample>_2_clean.fastq` or `<sample>_2_clean.fastq.gz` |
| Nanopore (optional) | `<sample>_ont.fastq` or `<sample>_ont.fastq.gz` |

Illumima data is mapped by BBMap (`minid=0.95`, `ambiguous=random`)

ONT data is available, mapped by minimap2

When both Illumina and ONT data are present for the same sample, BAMs are merged with `samtools merge` before depth extraction.


### Phage region coordinates (`-p`)
A tab-separated plain-text file **(no header)** listing predicted prophage regions: contig name; start; end
```
strainA_1      45231    78940
strainA_3      120000   155000
strainB_2      88200    119400
```

This file can be generated from any prophage prediction tool (e.g., PHASTER, VirSorter2, geNomad, DeepVirFinder) by extracting the relevant coordinate columns.

### Index table (`-i`, optional)
A tab-separated file **with a header** row mapping sequencing sample to reference genome. Required when sample identifiers do not match genome filenames, or when one sample should be mapped against multiple reference genomes.
```
sample    ref
sampleA   strainA
sampleB   strainB
sampleC   strainC
```

Thus, sampleA_1_clean.fastq/sampleA_2_clean.fastq/sampleA_ont.fastq will be mapped to reference StrainA.fasta

## Basic Usage
### Fast mode
```bash
bash active.sh \
    -g /path/to/genomes/ \
    -s /path/to/reads/ \
    -p prophage_regions.tsv \
    -o output_dir/
```
### Percies mode
```bash
bash active.sh \
    -g /path/to/genomes/ \
    -s /path/to/reads/ \
    -p prophage_regions.tsv \
    -o output_dir/ \
    -m
```
### With index table and custom parameters
```bash
bash active.sh \
    -g /path/to/genomes/ \
    -s /path/to/reads/ \
    -p prophage_regions.tsv \
    -o output_dir/ \
    -i sample_genome_index.tsv \
    -t 32 \
    -r 16 \
    --score 350 \
    --cov 0.80 \
    --iter 1000
```

### Remove intermediate file
```bash
# Delete all intermediate files; keep only final results
bash active.sh \
    -g /path/to/genomes/ \
    -s /path/to/reads/ \
    -p prophage_regions.tsv \
    -o output_dir/ \
    --rmall
```

## Parameters
### Basic required
```
  -g   Directory containing host genome(s); files must use .fasta extension
  -s   Directory containing sequencing data samples
         Short reads : *_1_clean.fastq[.gz]  and  *_2_clean.fastq[.gz]
         Long reads  : *_ont.fastq[.gz]  (optional)
  -p   Phage coordinate file (tab-separated: contig  start  end)
  -o   Output directory
  -t   Threads for mapping steps (default: 20)
  -r   Threads for R statistical analysis (default: 20)
```
### Mode
```
  -m   Run in precise mode (default: fast mode)
```
### Optional
```
  -i   Index table mapping sample names to reference genome names
         (tab-separated, first row = header)
  --score
       FetchMGs bit-score threshold; fast mode only (default: 300)
  --cov
       Minimum coverage fraction required for a phage region to be tested
       (default: 0.75; passed to the R statistics step)
  --iter
       Number of permutation iterations for Cliff's delta in precise mode
       (default: 500)
  --rmbam
       Delete BAM files after depth extraction to save disk space
       (default: keep)
  --rmdepth
       Delete per-region depth .txt files after R analysis
       (default: keep)
  --rmprodigal
       Delete Prodigal output directory after pipeline completes
       (fast mode only; default: keep)
  --rmmgs
       Delete all FetchMGs-related files and directories after pipeline
       completes (fast mode only; default: keep)
  --clean
       Shorthand for --rmbam --rmdepth: delete all BAM and depth files
  --rmall
       Delete ALL intermediate files and directories after the pipeline
       completes. Only the final activity_results/<mode>/ folder is kept.
       Equivalent to --rmbam --rmdepth --rmprodigal --rmmgs, and also
       removes the entire mapping_results/ tree (except the CSV summaries
       which have already been moved to activity_results/).
       (default: keep everything)
  -h, --help
       Show this help message and exit
```

## Output
```
output_dir/
├── activity_results/
│   └── <mode>/                                   # Final results (always retained)
│       └── <sample>_phage_depth_stats_summary_final.csv
├── mapping_results/
│   └── <sample>/
│       ├── <sample>_<ref>.bam                    # Aligned reads
│       ├── <sample>_<ref>.bam.bai
│       ├── <sample>_phage_basic_stats_summary.csv     # Descriptive statistics (fast mode)
│       ├── <sample>_basic_stats_summary.csv           # Descriptive statistics (precise mode)
│       ├── global/                               # Contig/prefix lists and region subsets
│       └── depth_to_stat/
│           └── <genome_prefix>/
│               ├── phage_depth_*.txt             # Per-position depth in prophage regions
│               ├── *_mgs_depth.txt               # Host background depth (fast mode)
│               └── *_host_nonphage_depth.txt     # Host background depth (precise mode)
├── prodigal_results/                             # Prodigal GFF/FAA (fast mode)
├── fetchmgs_results/                             # FetchMGs scores (fast mode)
├── mgs_ext_results/                              # MGS genomic coordinates (fast mode)
├── mgs_ordinate.tsv                              # Simplified MGS coordinate table
└── target_cogs_list.txt                          # COG IDs passing score threshold
```

## Output example
### Main output: Stastical test results
`<sample>_phage_depth_stats_summary_final.csv`. 

Key columns:

| Column | Description |
|------|------|
| Folder | Sample name |
| Genome| Host genome perfix |
| Phage_File | Prophage region |
| Coverage | Fraction of prophage region positions with depth > 0 |
| Test_Type| Statistical test applied (Student's t / Welch's t / Brunner-Munzel) |
| p_value | One-sided p-value (phage depth > host background) |
| Significant | `Yes` / `No` / `Coverage failed` (p_value= 0.01) |
| BM_Prob_Superiority| Brunner-Munzel probability of superiority (when applicable) |
| Q95 / Q99 | Whole genome 95th/99th depth percentiles (precise mode) |
| Phage_OverQ95 / Phage_OverQ99 | Sites in prophage region exceeding Q95/Q99 (precise mode) |
| RR_Q95 / RR_Q99 | Risk ratio relative to host background at Q95/Q99 (precise mode) |
| Delta_Mean| Mean Cliff's Δ from permutation analysis (precise mode, significant only) |
| Delta_CI_2.5 / Delta_CI_97.5 | 95% confidence interval of Cliff's Δ (precise mode, significant only) |
| Effect_Size_Type/ Effect_Size_Value | Effect size metric and value (fast mode) |

### Other output: Descriptive statistics
 `*_phage_basic_stats_summary.csv` / `*_basic_stats_summary.csv`

 Key columns:

| Column | Description |
|------|------|
| Folder | Sample name |
| Genome| Region type: `phage`, `host_mgs` (fast mode), or `host_non_phage` (precise mode) |
| File| Source depth file name |
| Mean| Mean depth of given region |
| Meadian| Median depth of given region |
| Variance| Variance of sequencing depth |
| COVERAGE| Fraction of positions with depth > 0 (only prophage region) |

A companion file recording depth summary statistics for every region processed, regardless of whether it passed the coverage threshold or hypothesis testing. This file is produced unconditionally alongside the test results. **Notion**: --rmall will remove this results

## Statistical Methods
Test selection is based on normality assessment using the Jarque-Bera test (alpha = 0.05) applied independently to the phage and host depth distributions:

Both normal: Levene's F-test for variance equality -> Student's t-test (equal variance) or Welch's t-test (unequal variance)
Either non-normal: Brunner-Munzel test (one-sided, testing phage > host)

All tests are one-sided with alpha = 0.01.
For significant Brunner-Munzel results (precise mode), effect size is estimated using permutation-based Cliff's delta:

Host depth vectors are subsampled (up to 10,000 positions per iteration) and Cliff's delta is computed against the full phage depth vector
The process is repeated --iter times (default: 500) in parallel
The mean and 95% bootstrap CI of delta are reported

For significant t-test results (fast mode), Hedges' g (bias-corrected Cohen's D) is reported.
A prophage region is excluded from hypothesis testing if fewer than --cov (default: 75%) of its positions have non-zero depth, as sparse coverage may indicate assembly artifacts or insufficient read recruitment rather than biological inactivity. Excluded regions are still recorded in the descriptive statistics file with their coverage value and depth summary statistics.
![image](https://github.com/SIAT-MaLab/ACTIVE/blob/main/figures1b.png)

## Citation
If you use active in your research, please cite:
> Y.H et al. A large-scale atlas of active gut prophages reveals significantly active prophage lineages in the human gut

## License
This project is licensed under the MIT License.
