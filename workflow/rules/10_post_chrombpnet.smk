rule concatenate_p53_peaks_w_atac:
  input:
    atac_regions = "results/atac/final/final_peaks.bed",
    p53_peaks = "results/p53_WT/data/peaks_inliers.bed"
  output:
    atac_w_p53 = "results/atac/chrombpnet_contribs/atac_w_p53.bed"
  shell:
    """
    ml load biology bedtools
    
    # Concatenate both peak sets
    cat {input.atac_regions} {input.p53_peaks} | \
      bedtools sort -i stdin > {output.atac_w_p53}
    """

rule chrombpnet_contributions:
  input:
    model = "results/atac/chrombpnet_model/models/chrombpnet_nobias.h5",
    atac_w_p53 = "results/atac/chrombpnet_contribs/atac_w_p53.bed",
    fasta = "resources/reference/genome.fa",
    chrom_sizes = "results/atac/chrombpnet_prep/chrom.subset.sizes"
  output:
    profile_h5 = "results/atac/chrombpnet_contribs/contribs.profile_scores.h5",
    profile_bw = "results/atac/chrombpnet_contribs/contribs.profile_scores.bw",
    counts_h5 = "results/atac/chrombpnet_contribs/contribs.counts_scores.h5",
    counts_bw = "results/atac/chrombpnet_contribs/contribs.counts_scores.bw",
    interpreted_regions = "results/atac/chrombpnet_contribs/contribs.interpreted_regions.bed"
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
    ml system
    ml cudnn/8.1.1.33
    ml cuda/11.2.0
    ml cairo/1.14.10
    ml pango/1.40.10
    
    # Derive output prefix from the profile_h5 output file
    OUTPUT_PREFIX=$(echo {output.profile_h5} | sed 's/\\.profile_scores\\.h5$//')
    
    mkdir -p $(dirname $OUTPUT_PREFIX)
    
    chrombpnet contribs_bw \
      -m {input.model} \
      -r {input.atac_w_p53} \
      -g {input.fasta} \
      -c {input.chrom_sizes} \
      -op $OUTPUT_PREFIX \
      -pc profile counts
    """

rule chrombpnet_predictions:
  input:
    bias_model = "resources/chrombpnet/bias_models/pretrained_bias_model.h5",
    chrombpnet_model = "results/atac/chrombpnet_model/models/chrombpnet.h5",
    chrombpnet_nobias = "results/atac/chrombpnet_model/models/chrombpnet_nobias.h5",
    regions = "results/atac/chrombpnet_contribs/atac_w_p53.bed",
    fasta = "resources/reference/genome.fa",
    chrom_sizes = "results/atac/chrombpnet_prep/chrom.subset.sizes"
  output:
    bias_bw = "results/atac/chrombpnet_preds/preds_bias.bw",
    bias_bed = "results/atac/chrombpnet_preds/preds_bias_preds.bed",
    chrombpnet_bw = "results/atac/chrombpnet_preds/preds_chrombpnet.bw",
    chrombpnet_bed = "results/atac/chrombpnet_preds/preds_chrombpnet_preds.bed",
    nobias_bw = "results/atac/chrombpnet_preds/preds_chrombpnet_nobias.bw",
    nobias_bed = "results/atac/chrombpnet_preds/preds_chrombpnet_nobias_preds.bed"
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
    ml system
    ml cudnn/8.1.1.33
    ml cuda/11.2.0
    ml cairo/1.14.10
    ml pango/1.40.10
    
    # Derive output prefix from bias_bw output file
    OUTPUT_PREFIX=$(echo {output.bias_bw} | sed 's/_bias\\.bw$//')
    
    mkdir -p $(dirname $OUTPUT_PREFIX)
    
    chrombpnet pred_bw \
      -bm {input.bias_model} \
      -cm {input.chrombpnet_model} \
      -cmb {input.chrombpnet_nobias} \
      -r {input.regions} \
      -g {input.fasta} \
      -c {input.chrom_sizes} \
      -op $OUTPUT_PREFIX
    """
