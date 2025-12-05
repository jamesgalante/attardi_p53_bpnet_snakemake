# NOTE: this is in MEFs - a differentiated, stable cell type
# It is not reflective of the dynamics of p53 binding

rule download_atac_fastq:
  output:
    fastq1 = "results/atac/raw_fastq/{srr}_1.fastq",
    fastq2 = "results/atac/raw_fastq/{srr}_2.fastq"
  threads: 4
  resources:
    mem = "16G",
    time = "2:00:00"
  shell:
    """
    ml load biology sra-tools/3.0.7
    mkdir -p $(dirname {output.fastq1})
    fasterq-dump --threads {threads} --outdir $(dirname {output.fastq1}) {wildcards.srr}
    """

rule build_bowtie2_index:
  input: fasta = "resources/reference/genome.fa"
  output: marker = "resources/reference/index_bowtie2/build.done"
  threads: 16
  resources:
    mem = "32G",
    time = "4:00:00"
  shell:
    """
    ml load biology bowtie2/2.5.4

    OUT_DIR=$(dirname {output.marker})
    mkdir -p $OUT_DIR

    bowtie2-build \
      --threads {threads} \
      {input.fasta} \
      $OUT_DIR/genome

    touch {output.marker}
    """

# https://www.htslib.org/doc/samtools-markdup.html
rule align_atac_reads:
  input:
    fastq1 = "results/atac/raw_fastq/{srr}_1.fastq",
    fastq2 = "results/atac/raw_fastq/{srr}_2.fastq",
    idx_marker = "resources/reference/index_bowtie2/build.done"
  output:
    bam = "results/atac/aligned/{srr}.sorted.raw.bam",
    flagstat_raw = "results/atac/aligned/{srr}.flagstat_raw.txt"
  threads: 16
  resources:
    mem = "32G",
    time = "4:00:00"
  shell:
    """
    ml load biology bowtie2/2.5.4 samtools/1.16.1

    IDX_PREFIX=$(dirname {input.idx_marker})/genome

    bowtie2 \
      -p {threads} \
      -x $IDX_PREFIX \
      -1 {input.fastq1} \
      -2 {input.fastq2} \
    | samtools fixmate -m - - \
    | samtools sort -@ {threads} -o {output.bam} -

    samtools flagstat {output.bam} > {output.flagstat_raw}
    """

rule filter_dedup_atac:
  input:
    bam = "results/atac/aligned/{srr}.sorted.raw.bam"
  output:
    final_bam = "results/atac/aligned/{srr}.nodup.bam",
    final_bai = "results/atac/aligned/{srr}.nodup.bam.bai",
    flagstat = "results/atac/aligned/{srr}.nodup.flagstat.txt"
  threads: 4
  resources:
    mem = "8G",
    time = "2:00:00"
  shell:
    """
    ml load biology samtools/1.16.1

    samtools markdup -r {input.bam} - \
    | samtools view -b -F 1804 -f 2 -q 30 -o {output.final_bam} -

    samtools index {output.final_bam}
    samtools flagstat {output.final_bam} > {output.flagstat}
    """

rule pool_replicates:
  input: 
    bams = lambda w: expand("results/atac/aligned/{srr}.nodup.bam", srr=config["other_MEF_tracks"]["reference"]["atac"])
  output: 
    pooled_bam = "results/atac/merged/pooled_reps.bam"
  threads: 8
  resources:
    mem = "16G",
    time = "2:00:00"
  shell:
    """
    ml load biology samtools/1.16.1

    samtools merge \
      -@ {threads} \
      {output.pooled_bam} \
      {input.bams}

    samtools index {output.pooled_bam}
    """
    
rule generate_atac_bigwig:
  input:
    bam = "results/atac/merged/pooled_reps.bam",
    sizes = "resources/reference/chrom.sizes"
  output:
    bw = "results/atac/bigwig/merged.bw",
    bg = temp("results/atac/bigwig/merged.bedGraph")
  resources:
    mem = "16G",
    time = "4:00:00"
  shell:
    """
    ml load biology bedtools/2.30.0 samtools/1.16.1 ucsc-utils
    
    mkdir -p $(dirname {output.bw})
    
    samtools view -b {input.bam} $(cut -f 1 {input.sizes}) | \
      bedtools genomecov -ibam stdin -bg | \
      sort -k1,1 -k2,2n > {output.bg}
    
    bedGraphToBigWig {output.bg} {input.sizes} {output.bw}
    """

rule call_peaks_pooled:
  input: 
    bam = "results/atac/merged/pooled_reps.bam"
  output: 
    peaks = "results/atac/peaks/pooled_peaks.narrowPeak"
  params:
    g_size = config["bpnet_params"]["macs_g"]
  conda: "../envs/macs3.yaml"
  shell:
    """
    macs3 callpeak \
      -t {input.bam} \
      -f BAMPE \
      -n pooled \
      -g {params.g_size} \
      -p 0.01 \
      --keep-dup all \
      --outdir $(dirname {output.peaks})
    """

rule filter_blacklist_pooled:
  input:
    peaks = "results/atac/peaks/pooled_peaks.narrowPeak",
    blacklist = "resources/reference/blacklist.bed",
    sizes = "resources/reference/chrom.sizes"
  output: 
    final_peaks = "results/atac/final/final_peaks.unfiltered.bed",
    extended_blacklist = temp("results/atac/final/extended_blacklist.bed")
  shell:
    """
    ml load biology bedtools
    
    # Extend blacklist by 1057 bp on each side
    bedtools slop -i {input.blacklist} -g {input.sizes} -b 1057 > {output.extended_blacklist}
    
    # Filter peaks
    bedtools intersect -v -a {input.peaks} -b {output.extended_blacklist} > {output.final_peaks}
    """

# rule plot_overlap:
#   input:
#     atac = "results/atac/final/final_peaks.bed",
#     p53_full = "results/p53_WT/data/macs2/full_peaks.narrowPeak",
#     p53_idr = "results/p53_WT/data/idr_peaks.bed"
#   output:
#     atac_statistics = "atac_stats.pdf",
#     p53_statistics = "p53_stats.pdf",
#     open_p53 = "open_p53.bed",
#     closed_p53 = "closed_p53.bed",
#     open_no_p53 = "open_no_p53.bed"
#   log: "p53_strand_bias.log"
#   conda: "../envs/overlapping_analysis.yaml"
#   resources:
#     mem = "8G",
#     time = "1:00:00"
#   script: "../scripts/plot_atac_p53_overlap.R"
# 
# rule plot_p53_strand_bias:
#   input:
#     p53_full = "results/p53_WT/data/macs2/full_peaks.narrowPeak",
#     p53_idr = "results/p53_WT/data/idr_peaks.bed",
#     p53_plus = "results/p53_WT/data/experiment_plus.bw.bw",
#     p53_minus = "results/p53_WT/data/experiment_minus.bw",
#     atac = "results/atac/final/final_peaks.bed",
#     control_plus = "results/p53_WT/data/control_plus.bw",
#     control_minus = "results/p53_WT/data/control_minus.bw"
#   output:
#     plot_full = "results/plots/p53_strand_bias_full.png",
#     plot_idr = "results/plots/p53_strand_bias_idr.png",
#     signal_data_full = "results/stats/p53_strand_bias_full_signal.csv",
#     signal_data_idr = "results/stats/p53_strand_bias_idr_signal.csv",
#     summary_stats = "results/stats/p53_strand_bias_summary.csv"
#   log: "p53_strand_bias.log"
#   conda: "../envs/overlapping_analysis.yaml"
#   resources:
#     mem = "8G",
#     time = "1:00:00"
#   script: "../scripts/plot_p53_strand_bias.R"
