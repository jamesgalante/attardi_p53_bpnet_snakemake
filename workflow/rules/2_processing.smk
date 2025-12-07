rule generate_bigwigs:
  input:
    exp_bam = lambda w: f"results/aligned/{config['samples'][w.sample]['experiment']}.sorted.bam",
    ctl_bam = lambda w: f"results/aligned/{config['samples'][w.sample]['control']}.sorted.bam",
    sizes = "resources/reference/chrom.sizes"
  output:
    exp_plus = "results/{sample}/data/experiment_plus.bw",
    exp_minus = "results/{sample}/data/experiment_minus.bw",
    ctl_plus = "results/{sample}/data/control_plus.bw",
    ctl_minus = "results/{sample}/data/control_minus.bw",
    ep_bg = temp("results/{sample}/data/experiment_plus.bedGraph"),
    em_bg = temp("results/{sample}/data/experiment_minus.bedGraph"),
    cp_bg = temp("results/{sample}/data/control_plus.bedGraph"),
    cm_bg = temp("results/{sample}/data/control_minus.bedGraph")
  resources:
    mem = "16G",
    time = "4:00:00"
  shell:
    """
    ml load biology bedtools/2.30.0 samtools/1.16.1 ucsc-utils

    samtools view -b {input.exp_bam} $(cut -f 1 {input.sizes}) | \
      bedtools genomecov -5 -bg -strand + -ibam stdin | \
      sort -k1,1 -k2,2n > {output.ep_bg}

    samtools view -b {input.exp_bam} $(cut -f 1 {input.sizes}) | \
      bedtools genomecov -5 -bg -strand - -ibam stdin | \
      sort -k1,1 -k2,2n > {output.em_bg}

    samtools view -b {input.ctl_bam} $(cut -f 1 {input.sizes}) | \
      bedtools genomecov -5 -bg -strand + -ibam stdin | \
      sort -k1,1 -k2,2n > {output.cp_bg}

    samtools view -b {input.ctl_bam} $(cut -f 1 {input.sizes}) | \
      bedtools genomecov -5 -bg -strand - -ibam stdin | \
      sort -k1,1 -k2,2n > {output.cm_bg}

    bedGraphToBigWig {output.ep_bg} {input.sizes} {output.exp_plus}
    bedGraphToBigWig {output.em_bg} {input.sizes} {output.exp_minus}
    bedGraphToBigWig {output.cp_bg} {input.sizes} {output.ctl_plus}
    bedGraphToBigWig {output.cm_bg} {input.sizes} {output.ctl_minus}
    """

rule generate_full_bigwigs:
  input:
    exp_bam = lambda w: f"results/aligned/{config['samples'][w.sample]['experiment']}.sorted.bam",
    ctl_bam = lambda w: f"results/aligned/{config['samples'][w.sample]['control']}.sorted.bam",
    sizes = "resources/reference/chrom.sizes"
  output:
    exp_plus = "results/{sample}/data/experiment_plus_full.bw",
    exp_minus = "results/{sample}/data/experiment_minus_full.bw",
    ctl_plus = "results/{sample}/data/control_plus_full.bw",
    ctl_minus = "results/{sample}/data/control_minus_full.bw",
    ep_bg = temp("results/{sample}/data/experiment_plus_full.bedGraph"),
    em_bg = temp("results/{sample}/data/experiment_minus_full.bedGraph"),
    cp_bg = temp("results/{sample}/data/control_plus_full.bedGraph"),
    cm_bg = temp("results/{sample}/data/control_minus_full.bedGraph")
  resources:
    mem = "16G",
    time = "4:00:00"
  shell:
    """
    ml load biology bedtools/2.30.0 samtools/1.16.1 ucsc-utils
    samtools view -b {input.exp_bam} $(cut -f 1 {input.sizes}) | \
      bedtools genomecov -bg -strand + -ibam stdin | \
      sort -k1,1 -k2,2n > {output.ep_bg}
    samtools view -b {input.exp_bam} $(cut -f 1 {input.sizes}) | \
      bedtools genomecov -bg -strand - -ibam stdin | \
      sort -k1,1 -k2,2n > {output.em_bg}
    samtools view -b {input.ctl_bam} $(cut -f 1 {input.sizes}) | \
      bedtools genomecov -bg -strand + -ibam stdin | \
      sort -k1,1 -k2,2n > {output.cp_bg}
    samtools view -b {input.ctl_bam} $(cut -f 1 {input.sizes}) | \
      bedtools genomecov -bg -strand - -ibam stdin | \
      sort -k1,1 -k2,2n > {output.cm_bg}
    bedGraphToBigWig {output.ep_bg} {input.sizes} {output.exp_plus}
    bedGraphToBigWig {output.em_bg} {input.sizes} {output.exp_minus}
    bedGraphToBigWig {output.cp_bg} {input.sizes} {output.ctl_plus}
    bedGraphToBigWig {output.cm_bg} {input.sizes} {output.ctl_minus}
    """
