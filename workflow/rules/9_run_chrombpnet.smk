
# Subset chromosome sizes for ChromBPNet
rule subset_chrom_sizes:
  input:
    chrom_sizes = "resources/reference/chrom.sizes"
  output:
    subset = "results/atac/chrombpnet_prep/chrom.subset.sizes"
  shell:
    """
    grep -v "chrM" {input.chrom_sizes} > {output.subset}
    """

# Generate chromosome splits for train/val/test
rule generate_splits:
  input:
    subset_chrom_sizes = "results/atac/chrombpnet_prep/chrom.subset.sizes"
  output:
    fold0_json = "results/atac/chrombpnet_prep/splits/fold_0.json"
  conda: "../envs/chrombpnet.yaml"
  shell:
    """
    OUTPUT_PREFIX=$(echo {output.fold0_json} | sed 's/\\.json$//')
    
    chrombpnet prep splits \
      -c {input.subset_chrom_sizes} \
      -tcr chr1 chr8 chr9 \
      -vcr chr2 chr10 \
      -op $OUTPUT_PREFIX
    """

# Generate GC-matched non-peak background regions
rule generate_nonpeaks:
  input:
    fasta = "resources/reference/genome.fa",
    peaks = "results/atac/final/final_peaks.bed",
    chrom_sizes = "results/atac/chrombpnet_prep/chrom.subset.sizes",
    folds = "results/atac/chrombpnet_prep/splits/fold_0.json",
    blacklist = "resources/reference/blacklist.bed"
  output:
    nonpeaks = "results/atac/chrombpnet_prep/nonpeaks/output_negatives.bed"
  conda: "../envs/chrombpnet.yaml"
  resources:
    mem = "16G",
    time = "2:00:00"
  shell:
    """
    OUTPUT_PREFIX=$(echo {output.nonpeaks} | sed 's/_negatives\\.bed$//')
    
    chrombpnet prep nonpeaks \
      -g {input.fasta} \
      -p {input.peaks} \
      -c {input.chrom_sizes} \
      -fl {input.folds} \
      -br {input.blacklist} \
      -o $OUTPUT_PREFIX
    """

# Download pretrained bias model
rule download_pretrained_bias_model:
  output:
    pretrained_bias_model = "resources/chrombpnet/bias_models/ENCSR868FGK_bias_fold_0.h5"
  params:
    url = config['other_MEF_tracks']['reference']['pretrained_bias_url']
  shell:
    """
    wget -q -O {output.pretrained_bias_model} {params.url}
    """
    
# Train ChromBPNet model
rule chrombpnet_pipeline:
  input:
    merged_bam = "results/atac/merged/pooled_reps.bam",
    fasta = "resources/reference/genome.fa",
    chrom_sizes = "results/atac/chrombpnet_prep/chrom.subset.sizes",
    peaks = "results/atac/final/final_peaks.bed",
    nonpeaks = "results/atac/chrombpnet_prep/nonpeaks/output_negatives.bed",
    folds = "results/atac/chrombpnet_prep/splits/fold_0.json",
    bias_model = "resources/chrombpnet/bias_models/ENCSR868FGK_bias_fold_0.h5"
  output:
    html_report = "results/atac/chrombpnet_model/evaluation/overall_report.html",
    pdf_report = "results/atac/chrombpnet_model/evaluation/overall_report.pdf",
    tf_model = "results/atac/chrombpnet_model/models/chrombpnet_nobias.h5",
    full_model = "results/atac/chrombpnet_model/models/chrombpnet.h5"
  conda: "../envs/chrombpnet.yaml"
  threads: 4
  resources:
   slurm_partition = "akundaje",
   gpu = 1,
   tasks_per_gpu = 0,
   runtime = "24h",
   mem = "32G"
  shell:
    """
    # Flush Python output immediately
    export PYTHONUNBUFFERED=1
    
    # Clean up previous run
    # rm -rf results/atac/chrombpnet_model
    
    # Load CUDA/cuDNN modules (adjust versions based on your cluster)
    ml system
    ml cudnn/8.1.1.33
    ml cuda/11.2.0
    ml cairo/1.14.10
    ml pango/1.40.10
    
    # Train ChromBPNet
    chrombpnet pipeline \
        -ibam {input.merged_bam} \
        -d ATAC \
        -g {input.fasta} \
        -c {input.chrom_sizes} \
        -p {input.peaks} \
        -n {input.nonpeaks} \
        -fl {input.folds} \
        -b {input.bias_model} \
        -o results/atac/chrombpnet_model
    """
