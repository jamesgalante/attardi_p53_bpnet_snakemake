# 1. Call Peaks on Full Sample
rule call_peaks_macs2_full:
  input:
    exp = lambda w: f"results/aligned/{config['samples'][w.sample]['experiment']}.sorted.bam",
    ctl = lambda w: f"results/aligned/{config['samples'][w.sample]['control']}.sorted.bam"
  output:
    peaks = "results/{sample}/data/macs2/full_peaks.narrowPeak"
  params:
    out_dir = "results/{sample}/data/macs2",
    g_size = config["bpnet_params"]["macs_g"]
  resources:
    mem = "16G",
    time = "2:00:00"
  conda: "../envs/preprocessing.yaml"
  shell:
    """
    macs2 callpeak \
      -t {input.exp} -c {input.ctl} \
      -f BAM -n full --outdir {params.out_dir} \
      -q 0.05 -g {params.g_size}
    """

# 2. Create Pseudoreplicates
rule create_pseudoreps:
  input:
    bam = lambda w: f"results/aligned/{config['samples'][w.sample][w.type]}.sorted.bam"
  output:
    pr1 = temp("results/{sample}/data/idr/{type}_pr00.bam"),
    pr2 = temp("results/{sample}/data/idr/{type}_pr01.bam")
  params:
    prefix = "results/{sample}/data/idr/{type}_pr"
  resources:
    mem = "16G",
    time = "2:00:00"
  threads: 4
  shell:
    """
    ml load biology samtools/1.16.1

    samtools view -H {input.bam} > {params.prefix}.header.sam
    
    nlines=$(samtools view -c {input.bam})
    nlines=$(( (nlines + 1) / 2 ))
    
    samtools view {input.bam} | shuf | split -d -l ${{nlines}} - {params.prefix}
    
    cat {params.prefix}.header.sam {params.prefix}00 | samtools view -bS - > {output.pr1}
    cat {params.prefix}.header.sam {params.prefix}01 | samtools view -bS - > {output.pr2}
    
    rm {params.prefix}00 {params.prefix}01 {params.prefix}.header.sam
    """

# 3. Call Peaks on Pseudoreplicates
rule call_peaks_macs2_pseudoreps:
  input:
    exp = "results/{sample}/data/idr/experiment_pr{rep}.bam",
    ctl = "results/{sample}/data/idr/control_pr{rep}.bam"
  output:
    peaks = temp("results/{sample}/data/idr/peaks_pr{rep}_peaks.narrowPeak")
  params:
    out_dir = lambda w, output: os.path.dirname(output.peaks),
    name = lambda w: f"peaks_pr{w.rep}",
    g_size = config["bpnet_params"]["macs_g"]
  resources:
    mem = "16G",
    time = "2:00:00"
  conda: "../envs/preprocessing.yaml"
  shell:
    """
    macs2 callpeak \
      -t {input.exp} -c {input.ctl} \
      -f BAM -n {params.name} --outdir {params.out_dir} \
      -p 1e-3 -g {params.g_size}
    """

# 4. Run IDR and Filter
rule run_idr:
  input:
    p1 = "results/{sample}/data/idr/peaks_pr00_peaks.narrowPeak",
    p2 = "results/{sample}/data/idr/peaks_pr01_peaks.narrowPeak"
  output:
    final_peaks = "results/{sample}/data/idr_peaks.bed",
    plot = "results/{sample}/data/idr/idr_plot.png"
  params:
    idr_out = "results/{sample}/data/idr/idr_results.txt"
  resources:
    mem = "16G",
    time = "2:00:00"
  conda: "../envs/preprocessing.yaml"
  shell:
    """
    sort -k8,8nr {input.p1} > {input.p1}.sorted
    sort -k8,8nr {input.p2} > {input.p2}.sorted
    
    idr --samples {input.p1}.sorted {input.p2}.sorted \
      --input-file-type narrowPeak --rank p.value \
      --output-file {params.idr_out} --plot
    mv {params.idr_out}.png {output.plot}

    awk '$5 >= 540' {params.idr_out} | cut -f 1-10 > {output.final_peaks}
    """
