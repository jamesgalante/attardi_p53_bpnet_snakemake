# BPNet and ChromBPNet Analysis of p53 Binding in Mouse Embryonic Fibroblasts

## Overview

This pipeline trains BPNet on p53 ChIP-seq data and ChromBPNet on ATAC-seq data from mouse embryonic fibroblasts (MEFs) to investigate p53 transcription factor binding patterns and chromatin accessibility.

**Biological Context:** Previous ChromBPNet analyses on pancreatic pseudo-bulked datasets did not identify p53 as a discovered motif through TF-MoDISco, despite expected binding at known p53-regulated genes (CDKN1A, BAX). Examination of sequence attributions at canonical p53 binding sites showed minimal sequence importance scores across all conditions tested. This unexpected finding prompted an investigation using orthogonal datasets from MEFs, where p53 binding is well-characterized, to determine whether this represents a biological phenomenon or a technical limitation of the attribution methods.

**Results:** Following BPNet training on ChIP-seq and ChromBPNet training on ATAC-seq from MEFs, quality control metrics from the ChromBPNet evaluation reports showed good model performance. However, TF-MoDISco analysis on the initial subsampled dataset again failed to identify the p53 motif. Visualization of ChromBPNet sequence attributions at ChIP-validated p53 binding sites in IGV revealed the same pattern observed in pancreatic datasets: p53-bound sequences showed no notable sequence attribution scores, suggesting this may reflect fundamental aspects of p53 binding dynamics, however technical limitations due to the models themselves or composition of the training data have not been ruled out.

The final outputs from this pipeline enable systematic comparison of model predictions, contribution scores, and motif discovery results across ChIP-seq and ATAC-seq modalities.

## Repository Structure

```
attardi_p53_bpnet_snakemake/
├── config/
│   ├── config.yaml               # Pipeline configuration parameters
│   └── config_info.txt           # Configuration documentation/reasoning
├── workflow/
│   ├── Snakefile                 # Main workflow definition
│   ├── rules/                    # Snakemake rule definitions
│   │   ├── 0_references.smk      # Reference genome and indices
│   │   ├── 1_alignment.smk       # Read alignment
│   │   ├── 2_processing.smk      # Signal track generation
│   │   ├── 3_peaks.smk           # Peak calling and IDR
│   │   ├── 4_inputs.smk          # BPNet input preparation
│   │   ├── 5_training.smk        # BPNet model training
│   │   ├── 6_denoising_and_attributions.smk  # Predictions and SHAP
│   │   ├── 7_modisco.smk         # TF-MoDISco motif discovery
│   │   ├── 8_process_atac.smk    # ATAC-seq processing
│   │   ├── 9_run_chrombpnet.smk  # ChromBPNet training
│   │   └── 10_post_chrombpnet.smk # ChromBPNet contributions
│   ├── scripts/                  # Analysis scripts
│   │   ├── plot_atac_p53_overlap.R
│   │   └── plot_p53_strand_bias.R
│   └── envs/                     # Conda environment definitions
│       ├── bpnet.yaml
│       ├── chrombpnet.yaml
│       ├── macs3.yaml
│       ├── preprocessing.yaml
│       ├── tfmodisco.yaml
│       └── plotting.yaml
├── resources/                    # Downloaded reference files
├── results/                      # Pipeline outputs
└── README.md
```

## Requirements

### Software Dependencies

- **Snakemake:** 9.13.7
- **Python:** 3.13.9
- **Miniforge3:** For environment management
- **Snakemake executor** 1.9.2

### Environment Management

Dependencies are managed automatically through conda environments defined in `workflow/envs/`. Snakemake will create and activate the appropriate environments for each step.

## Configuration

### Critical Configuration Note

⚠️ **ChromBPNet Bias Model Requirement**

The pipeline requires a pretrained ChromBPNet bias model for mouse ATAC-seq. Currently configured to use a local filepath on Stanford's Oak storage:

```yaml
other_MEF_tracks:
  reference:
    pretrained_bias_url: "/oak/stanford/groups/akundaje/zhangby/mouse-atlas/bias_model/ATAC/ENCSR012YAB/fold_0/threshold_0.9/models/ENCSR012YAB_bias.h5"
```

**TODO:** Update to a publicly accessible URL. For now, ensure you have access to this filepath or train your own bias model following ChromBPNet documentation.

### Main Configuration File

**Reference Genome:**
- `reference.name`: Genome assembly (default: mm10)
- `reference.fasta_url`: FASTA download URL
- `reference.chrom_sizes_url`: Chromosome sizes
- `reference.blacklist_url`: Genomic blacklist regions
- `reference.motif_url`: Motif database (HOCOMOCO)

**Samples:**
```yaml
samples:
  p53_WT:
    experiment: "SRR832846"  # ChIP-seq
    control: "SRR832847"     # Input control
  p53_KO:
    experiment: "SRR832848"
    control: "SRR832849"
```

**ATAC-seq Replicates:**
```yaml
other_MEF_tracks:
  reference:
    atac: ["SRR18371806", "SRR18371805", "SRR18371804"]  # ATAC-seq replicates
```

See `config/config.yaml` for complete parameter documentation.

## Running the Pipeline

### HPC Execution (Recommended)

This pipeline is designed to run on Stanford's Sherlock HPC cluster with SLURM scheduling:

```bash
# Clone repository
git clone https://github.com/jamesgalante/attardi_p53_bpnet_snakemake.git
cd attardi_p53_bpnet_snakemake

# Run complete pipeline
snakemake --profile slurm_profile all # See below for slurm_profile details
```

### SLURM Profile Configuration

Create a Snakemake profile at `~/.config/snakemake/slurm_profile/config.yaml`:

```yaml
# Set the Executor for Cluster Submission
executor: "slurm"

# Global Defaults
jobs: 500
retries: 0

# Environment Management
use-conda: True
conda-frontend: mamba
conda-prefix: path_to_snakemake_conda_envs
notemp: False

# 4. Default Resources (Passed to the SLURM Executor)
default-resources:
    - slurm_account=your_account
    - slurm_partition=your_partition,owners,normal
    - runtime="13h"
    - slurm_extra="--nice"
```

## Pipeline Stages

### 1. Reference Preparation (`0_references.smk`)
- Downloads mm10 reference genome and indices
- Builds Bowtie and Bowtie2 indices
- Downloads motif databases and blacklist regions

### 2. ChIP-seq Alignment (`1_alignment.smk`)
- Downloads p53 ChIP-seq and input control samples
- Aligns reads using Bowtie
- Filters alignments (MAPQ ≥30, removes duplicates)

### 3. Signal Processing (`2_processing.smk`)
- Generates strand-specific bigWig tracks
- Creates 5' cut-site and full-coverage tracks
- Normalizes signal for visualization

### 4. Peak Calling (`3_peaks.smk`)
- Calls peaks using MACS2 on full samples
- Creates pseudoreplicates for each sample
- Runs IDR to identify reproducible peaks
- Filters peaks using IDR threshold

### 5. BPNet Input Preparation (`4_inputs.smk`)
- Removes outlier peaks
- Computes GC content reference
- Generates GC-matched negative regions
- Creates training/validation/test datasets
- Computes loss weights for profile and counts

### 6. BPNet Training (`5_training.smk`)
- Trains BPNet models with bias correction
- Uses cross-validation on chromosome splits
- Generates profile and counts predictions

### 7. Predictions and Attribution (`6_denoising_and_attributions.smk`)
- Generates predictions on test chromosomes
- Creates genome-wide denoised tracks
- Computes SHAP attribution scores
- Produces bigWig files for visualization

### 8. TF-MoDISco Analysis (`7_modisco.smk`)
- Runs TF-MoDISco on SHAP scores
- Discovers sequence motifs
- Generates HTML reports with motif annotations
- Compares to known motif databases

### 9. ATAC-seq Processing (`8_process_atac.smk`)
- Downloads and aligns paired-end ATAC-seq reads
- Removes duplicates and filters for quality
- Pools biological replicates
- Calls accessible chromatin peaks
- Filters blacklisted regions

### 10. ChromBPNet Training (`9_run_chrombpnet.smk`)
- Prepares chromosome splits for ChromBPNet
- Generates GC-matched non-peak regions
- Trains ChromBPNet model with pretrained bias
- Produces evaluation metrics and reports

### 11. ChromBPNet Analysis (`10_post_chrombpnet.smk`)
- Combines ATAC peaks with p53 ChIP peaks
- Computes contribution scores across peaks
- Generates profile and counts predictions
- Creates bigWig tracks for genome browser visualization

## TODO

- [ ] Update `pretrained_bias_url` to public URL when mouse ATAC-seq bias model is released
- [ ] Add optional plotting rules for ATAC/p53 overlap analysis
