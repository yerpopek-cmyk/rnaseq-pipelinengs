# RNA-Seq Data Analysis Pipeline

## Overview
This repository contains an end-to-end automated pipeline for RNA-Seq data analysis. It covers everything from raw read quality control to advanced analyses such as variant calling, gene fusion detection, and eQTL mapping.

The pipeline is implemented as a Bash script utilizing GNU Parallel for efficient multi-threading across samples.

## Features & Pipeline Steps
The pipeline is divided into 7 main stages:
1. **Raw Read QC & Alignment**: Runs FastQC on raw reads, then aligns them to the reference genome using **STAR** (2-pass mode).
2. **Expression Quantification**: Quantifies gene expression using **HTSeq-count** and **StringTie** (both reference-guided and de novo transcript assembly).
3. **BAM Preparation**: Prepares alignments for GATK4 variant calling by adding Read Groups, marking duplicates, splitting 'N' CIGAR reads (introns), and applying Base Quality Score Recalibration (BQSR).
4. **Variant Calling**: Calls variants using **GATK4 HaplotypeCaller** in GVCF mode, performs joint genotyping, and applies RNA-specific hard filtration.
5. **Post-Alignment QC**: Runs FastQC on recalibrated BAMs, **Qualimap** for RNA-Seq specific QC, and aggregates all metrics using **MultiQC**.
6. **Annotation & Fusion Detection**: Annotates variants with **SnpEff** and detects gene fusions using **Arriba**.
7. **ASE & eQTL**: Prepares per-sample VCFs, calculates Allele-Specific Expression (ASE) using **GATK ASEReadCounter**, and performs cis-eQTL mapping using **TensorQTL** and **plink2**.

## Prerequisites
The pipeline relies on several Conda environments to manage dependencies. Ensure the following environments are created and contain the respective tools:
- `STAR`: fastqc, parallel, star, samtools, htseq, stringtie
- `gatk4`: gatk4
- `QC_fastq`: fastqc, qualimap, bcftools, multiqc
- `snpEff`: snpeff, arriba, bcftools
- `pandas-plink`: plink2, tensorqtl, pandas

## Directory Structure
The script expects the following inputs in the working directory before execution:
- `fastq/` directory containing raw `.fastq.gz` read files.
- `STAR/` directory containing the STAR genome index, reference GTF, reference genome FASTA, and `expression_final.bed.gz` for TensorQTL.
- `reference/` directory containing known variants for BQSR (`common_all_20180418.vcf.gz`).

## How to Run
To run the pipeline, simply execute the main shell script. Ensure your Conda initialization is set up properly for non-interactive shells or run it interactively.

```bash
bash rnaseq_pipeline_v2.sh
```

## Outputs
All pipeline results will be saved into the `RNA_Seq/result/` folder. The outputs include:
- `fastqc/`, `qualimap/`: Quality control reports (viewable collectively via `multiqc_report.html`).
- `bam/`, `final_bam/`: Sorted and preprocessed alignments.
- `gvcf/`, `RNA.filt.vcf.gz`: Variant calling and filtration results.
- `fusions/`: Gene fusion predictions.
- `tensorqtl/`: eQTL mapping results.
- `*.counts`, `stringtie/`: Gene expression quantification results.

Only generalized outputs and code are tracked in this repository, intermediate heavy files (BAMs, raw FASTQs) are excluded via `.gitignore`.
