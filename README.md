# Autumn–winter progression in Arabidopsis

Analysis scripts accompanying Ormancey et al. (2026), *Autumn–winter progression triggers distinct transcriptional programs in Arabidopsis thaliana*. [Genome Biology](https://doi.org/10.1186/s13059-026-04289-3).

The scripts cover RNA-seq processing, differential expression, expression-pattern classification, temperature/time modelling, ATAC-seq processing, motif enrichment and TF networks. Run commands from this directory. Outputs go into `results/`.

`scripts/` contains the individual analysis steps. Input tables, cluster assignments and motif sequences are in `data/reference/`. Sample information is in `metadata/`.

## Data and setup

Raw sequencing data are available under [GSE283018](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE283018), comprising [RNA-seq: GSE283016](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE283016) and [ATAC-seq: GSE283017](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE283017).

To start with the processed data, download `GeneRPKM.tsv`, `GeneCounts.tsv` and `atac_peak_rpkm.tsv` from the supplementary files on those pages, remove any accession prefix from the filenames, and place them in `data/processed/`. Alternatively:

```bash
bash scripts/download_processed.sh
```

Install the R packages with:

```bash
Rscript scripts/install_packages.R
```

The raw-read scripts additionally require cutadapt, STAR, StringTie and its `prepDE.py3` script for RNA-seq; fastp, Bowtie 2, samtools, bedtools, MACS3 and UCSC `bedGraphToBigWig` for ATAC-seq. Motif discovery/enrichment uses MEME Suite (`streme`, `sea`, `centrimo`, `tomtom`). Python 3 is used to combine ATAC counts.

## Follow the analyses

These commands use the processed expression matrices and the inputs in `data/reference/`:

```bash
Rscript scripts/03_expression_patterns.R
Rscript scripts/04_rnaseq_figures.R
Rscript scripts/05_temperature_time.R

Rscript scripts/07_atac_figures.R
Rscript scripts/09_tf_networks.R
Rscript scripts/10_motif_figures.R
Rscript scripts/11_gene_ontology.R

Rscript scripts/12_tf_enrichment.R
Rscript scripts/13_rnaseq_supplement.R
Rscript scripts/14_go_figures.R
```

The pattern script reads correlations and model statistics from `data/reference/rna_pattern_statistics.tsv.gz`; add `--refit` to calculate them from counts. Heatmaps use the cluster assignments in `data/reference/`. Motif and GO figures read the enrichment tables there; `11_gene_ontology.R --refresh` reruns the GO queries.

Script 04 also produces cluster-average profiles, the VIN3 comparison with Liu et al., VIN3-correlated gene profiles and the U1 gene profiles. Script 12 plots TF motif enrichment against expression changes. Scripts 13 and 14 produce the RNA-seq correlation and differential-expression plots, condition-specific GO plots and ATAC-cluster GO plots.

To repeat normalisation, differential expression or motif discovery:

```bash
Rscript scripts/01_rnaseq_normalization.R
Rscript scripts/02_rnaseq_differential_expression.R
Rscript scripts/06_atac_normalization.R

bash scripts/08_motif_tools.sh
```

The normalisation scripts write matrices into `results/`. To use these in subsequent analyses, set `RNA_RPKM=results/rnaseq/gene_rpkm_matrix_ruvseq.tsv` or `ATAC_RPKM=results/atac/atac_peak_counts.tsv` before running the relevant R script. The motif shell script writes its results into `results/motifs/`; the figure and network scripts read motif inputs from `data/reference/` by default.

Known TF motif enrichment uses the JASPAR2024 Arabidopsis motifs. To use enrichment results from script 08 in script 12, set `TF_ENRICHMENT_DIR=results/motifs/known_tf`.

## Figure outputs

The scripts produce individual plots for the panels listed below. Output paths are relative to `results/`. Entries marked “RNA only” cover the RNA-seq profiles within a panel; the Figure 1F entry provides the cluster-average expression curves.

| Figure panel | Script | Output files |
|---|---|---|
| 1D | [04_rnaseq_figures.R](scripts/04_rnaseq_figures.R) | `rnaseq/figures/Up_regulated_genes.pdf`, `Down_regulated_genes.pdf` |
| 1E | [04_rnaseq_figures.R](scripts/04_rnaseq_figures.R), [11_gene_ontology.R](scripts/11_gene_ontology.R) | `rnaseq/figures/RNAseq_clusters.pdf`; `gene_ontology/RNAseq_clusters_GO.pdf` includes GO word clouds |
| 1F — expression curves | [04_rnaseq_figures.R](scripts/04_rnaseq_figures.R) | `rnaseq/figures/cluster_mean_profiles.pdf` |
| 2A | [04_rnaseq_figures.R](scripts/04_rnaseq_figures.R) | `rnaseq/figures/EMB2746.pdf`, `KIN1.pdf`, `RD29B.pdf` |
| 2D | [04_rnaseq_figures.R](scripts/04_rnaseq_figures.R) | `rnaseq/figures/Liu2022_Inc_cold_*.pdf` |
| 3A, 3B | [04_rnaseq_figures.R](scripts/04_rnaseq_figures.R) | `rnaseq/figures/VIN3.pdf`, `Liu2022_VIN3.pdf` |
| 3C | [05_temperature_time.R](scripts/05_temperature_time.R) | `temperature_time/AllGenes.pdf` |
| 4A | [04_rnaseq_figures.R](scripts/04_rnaseq_figures.R) | `rnaseq/figures/CBF1.pdf`, `CBF2.pdf`, `CBF3.pdf` |
| 5A | [07_atac_figures.R](scripts/07_atac_figures.R) | `atac/figures/ATACseq_clusters.pdf` |
| 5B | [10_motif_figures.R](scripts/10_motif_figures.R) | `motifs/figures/motif_enrichment.pdf`, `motif_logos.pdf` |
| 5C | [12_tf_enrichment.R](scripts/12_tf_enrichment.R) | `motifs/tf_enrichment/TF_enrichment.pdf` |
| 5G | [09_tf_networks.R](scripts/09_tf_networks.R) | `networks/COOLING_network.pdf`, `CHILLING_network.pdf` |
| 6A | [04_rnaseq_figures.R](scripts/04_rnaseq_figures.R) | `rnaseq/figures/CCA1.pdf`, `ELF4.pdf`, `COR28.pdf`, `HY5.pdf`, `HYH.pdf`, `COR27.pdf` |
| 6B | [14_go_figures.R](scripts/14_go_figures.R) | `gene_ontology/figures/AVG_GO.pdf`, `FLU_GO.pdf` |
| 7A | [07_atac_figures.R](scripts/07_atac_figures.R) | `atac/figures/top500_up.pdf`, `top500_down.pdf` |
| S2A–C | [13_rnaseq_supplement.R](scripts/13_rnaseq_supplement.R) | `rnaseq/supplement/sample_correlations.pdf`, `initial_timepoint_comparisons.pdf`, `AVG_FLU_weekly_comparisons.pdf` |
| S3A | [04_rnaseq_figures.R](scripts/04_rnaseq_figures.R) | `rnaseq/figures/Inc_cold_*.pdf`, `Dec_cold_*.pdf` |
| S3B | [04_rnaseq_figures.R](scripts/04_rnaseq_figures.R) | `rnaseq/figures/VIN3_correlated_*.pdf` |
| S3C | [05_temperature_time.R](scripts/05_temperature_time.R) | `temperature_time/AllGenes_rev.pdf`, `early.pdf`, `linear.pdf`, `late.pdf` |
| S4 — RNA only | [04_rnaseq_figures.R](scripts/04_rnaseq_figures.R) | `rnaseq/figures/HCAR.pdf`, `SCPL49.pdf`, `RD29A.pdf`, `TPR.pdf`, `EARLI1.pdf`, `AZI1.pdf` |
| S5 | [04_rnaseq_figures.R](scripts/04_rnaseq_figures.R) | `rnaseq/figures/U1_gene_profiles.pdf` |
| S7A | [14_go_figures.R](scripts/14_go_figures.R) | `gene_ontology/figures/C1_GO.pdf`, `C2_GO.pdf` |
| S7B — RNA only | [04_rnaseq_figures.R](scripts/04_rnaseq_figures.R) | `rnaseq/figures/KIN1.pdf` |
| S11A | [04_rnaseq_figures.R](scripts/04_rnaseq_figures.R) | `rnaseq/figures/CO.pdf`, `FKF1.pdf`, `PRR5.pdf`, `PRR9.pdf` |

Where several filenames share a directory, it is shown on the first filename. The scripts also write associated data tables alongside many of the plots.

## Start from raw reads

Download paired FASTQ files through the SRA links on the GEO pages. Match libraries using `metadata/rnaseq_samples.tsv` and `metadata/atac_samples.tsv`, and name the paired files according to their `read1` and `read2` columns:

```text
data/raw/rnaseq/JQ_2WS1_1.fq.gz
data/raw/rnaseq/JQ_2WS1_2.fq.gz

data/raw/atac/ATAC_cold_MO_19_1.fq.gz
data/raw/atac/ATAC_cold_MO_19_2.fq.gz
```

Each sample table records submitted filenames, condition, replicate and week. The `sample` column uses weeks since the start of cold; `analysis` gives the matrix column names used by the R scripts. There are 39 RNA-seq libraries and 16 ATAC-seq libraries.

Place the TAIR10 genome FASTA and Araport11 GFF3 in `data/genome/`, uncompressed. The FASTA chromosome names should match `data/reference/chrom.sizes` (`1`–`5`, `Mt`, `Pt`). Then set the file paths and run:

```bash
export GENOME_FASTA="$PWD/data/genome/TAIR10.fa"
export ARAPORT_GFF="$PWD/data/genome/Araport11.gff3"
export PREPDE="/path/to/stringtie/prepDE.py3"
export THREADS=6

bash scripts/00_rnaseq_reads.sh
bash scripts/00_atac_reads.sh

Rscript scripts/01_rnaseq_normalization.R results/rnaseq/gene_count_matrix.csv
Rscript scripts/02_rnaseq_differential_expression.R results/rnaseq/gene_count_matrix.csv
Rscript scripts/06_atac_normalization.R results/atac/counts.tsv
```

ATAC counts are collected over `data/reference/atac_peaks.bed`. Peak calls are also written to `results/atac/peaks/`.

The Liu et al. comparison uses the FPKM table from [Liu et al. (2022)](https://doi.org/10.3389/fpls.2022.983460), associated with PRJNA854395.

The MIT licence covers the analysis code. Reference data retain their source terms and attribution.

## Software and packages

Install the tools needed for the scripts you plan to run.

| Software | Used for |
|---|---|
| R, including `Rscript` | Statistical analyses and figures |
| Bash, `awk`, `gzip`/`gunzip`, `curl` and standard Unix utilities | Shell scripts and data downloads |
| Python 3 | StringTie count extraction and combining ATAC counts; no additional Python packages are needed |
| cutadapt | RNA-seq adapter trimming |
| STAR | RNA-seq alignment |
| StringTie and `prepDE.py3` | Transcript quantification and RNA-seq count matrices |
| fastp | ATAC-seq adapter trimming |
| Bowtie 2, including `bowtie2-build` | ATAC-seq alignment and genome indexing |
| samtools | BAM processing and duplicate removal |
| bedtools | Interval operations, coverage and peak counts |
| MACS3 | ATAC-seq peak calling |
| UCSC `bedGraphToBigWig` | BigWig track generation |
| MEME Suite: `streme`, `sea`, `centrimo`, `tomtom` | Motif discovery, enrichment and comparison |

**R packages from CRAN:** BiocManager, readr, dplyr, ggplot2, ggrepel, cowplot, igraph, matrixStats, gprofiler2, circlize and viridisLite.

**R packages from Bioconductor:** edgeR, RUVSeq, GenomicFeatures, rtracklayer, Biostrings, ComplexHeatmap, simplifyEnrichment, universalmotif and variancePartition.

The `grid` package is included with R. Install the R packages and their dependencies with `Rscript scripts/install_packages.R`.
