rule download_fastq:
  output:
    fastq = "results/raw_fastq/{srr}.fastq"
  threads: 4
  resources:
    mem = "16G",
    time = "2:00:00"
  shell:
    """
    ml load biology sra-tools/3.0.7
    mkdir -p $(dirname {output.fastq})
    fasterq-dump --threads {threads} --outdir $(dirname {output.fastq}) {wildcards.srr}
    """

rule align_reads:
  input:
    fastq = "results/raw_fastq/{srr}.fastq",
    idx_marker = "resources/reference/index/build.done"
  output:
    bam = "results/aligned/{srr}.sorted.raw.bam",
    bai = "results/aligned/{srr}.sorted.raw.bam.bai",
    flagstat_raw = "results/aligned/{srr}.flagstat_raw.txt" 
  threads: 16
  resources:
    mem = "16G",
    time = "2:00:00"
  shell:
    """
    ml load biology bowtie/1.2.2 samtools/1.16.1

    IDX_PREFIX=$(dirname {input.idx_marker})/genome

    bowtie --threads {threads} -S $IDX_PREFIX -q {input.fastq} | \\
    samtools sort -@ {threads} -o {output.bam} -
    
    samtools index {output.bam}

    samtools flagstat {output.bam} > {output.flagstat_raw}
    """

rule post_alignment_qc:
  input:
    bam = "results/aligned/{srr}.sorted.raw.bam"
  output:
    filtered_bam = "results/aligned/{srr}.sorted.bam",
    filtered_bai = "results/aligned/{srr}.sorted.bam.bai",
    flagstat_final = "results/aligned/{srr}.flagstat_final.txt",
  threads: 10
  resources:
    mem = "16G",
    time = "2:00:00"
  shell:
    """
    ml load biology samtools/1.16.1

    samtools view -F 772 -q 30 -b {input.bam} -o {output.filtered_bam}

    samtools index {output.filtered_bam}

    samtools flagstat {output.filtered_bam} > {output.flagstat_final}
    """
