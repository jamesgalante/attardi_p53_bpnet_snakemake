# NOTE: this is in MEFs - a differentiated, stable cell type
# It is not reflective of the dynamics of p53 binding

# I don't feel like preprocessing the ATAC data yet
# Might have to do if this turns out to be interesting
# STEP 1
rule download_atac_seq_data:
  output: "resources/Other_MEF_Tracks/atac.bw"
  params: 
    url = config['other_MEF_tracks']['reference']['atac']
  shell: 'wget -O {output} {params.url}'

# Actually... since they only provide the bigwig, I might have to reanalyze the data
# Yeah - have to download each bam from SRA
# Then have to process thsoe bam files and peak call with macs2

# plot statistics of atac seq data and p53 data
# STEP 2
rule plot_preliminary_statistics:
  input:
    atac = "atac_bed_file.bed",
    p53 = "p53_bed_file.bed"
  output:
    atac_statistics = "atac_stats.pdf",
    p53_statistics = "p53_stats.pdf"
  conda: "../envs/overlapping_analysis.yaml"
  resources:
    mem = "8G",
    time = "1:00:00"
  script: "../scripts/plot_preliminary_statistics.R"

# overlap atac seq and p53 chip and plot
# STEP 3
rule plot_overlap:
  input:
    atac = "atac_bed_file.bed",
    p53 = "p53_bed_file.bed"
  output:
    atac_statistics = "atac_stats.pdf",
    p53_statistics = "p53_stats.pdf",
    open_p53 = "open_p53.bed",
    closed_p53 = "closed_p53.bed",
    open_no_p53 = "open_no_p53.bed"
  conda: "../envs/overlapping_analysis.yaml"
  resources:
    mem = "8G",
    time = "1:00:00"
  script: "../scripts/plot_overlap.R"
