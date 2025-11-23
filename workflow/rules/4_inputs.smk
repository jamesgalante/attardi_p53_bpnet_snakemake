rule prepare_outlier_json:
  input:
    plus = "results/{sample}/data/experiment_plus.bw",
    minus = "results/{sample}/data/experiment_minus.bw",
    ctl_plus = "results/{sample}/data/control_plus.bw",
    ctl_minus = "results/{sample}/data/control_minus.bw",
    peaks = "results/{sample}/data/peaks.bed"
  output:
    json = "results/{sample}/configs/input_outliers.json"
  run:
    import json
    data = {
      "0": {
        "signal": {"source": [input.plus, input.minus]},
        "loci": {"source": [input.peaks]},
        "bias": {
          "source": [input.ctl_plus, input.ctl_minus],
          "smoothing": [None, None]
        }
      }
    }
    with open(output.json, 'w') as f:
      json.dump(data, f, indent=4)

rule remove_outliers:
  input:
    json = "results/{sample}/configs/input_outliers.json",
    sizes = "resources/reference/chrom.sizes",
    chroms = "resources/reference/chroms.txt",
    blacklist = "resources/reference/blacklist.bed"
  output:
    peaks_inliers = "results/{sample}/data/peaks_inliers.bed"
  params:
    quantile = config["bpnet_params"]["outlier_quantile"],
    scale = config["bpnet_params"]["outlier_scale_factor"],
    seq_len = config["bpnet_params"]["output_len"]
  conda: "../envs/bpnet.yaml"
  shell:
    """
    bpnet-outliers \
      --input-data {input.json} \
      --quantile {params.quantile} \
      --quantile-value-scale-factor {params.scale} \
      --task 0 \
      --chrom-sizes {input.sizes} \
      --chroms $(paste -s -d ' ' {input.chroms}) \
      --sequence-len {params.seq_len} \
      --blacklist {input.blacklist} \
      --global-sample-weight 1.0 \
      --output-bed {output.peaks_inliers}
    """

rule gc_reference:
  input:
    fasta = "resources/reference/genome.fa",
    sizes = "resources/reference/chrom.sizes"
  output:
    gc_bed = "resources/reference/genomewide_gc.bed"
  params:
    input_len = config["bpnet_params"]["input_len"],
    out_prefix = "resources/reference/genomewide_gc"
  conda: "../envs/bpnet.yaml"
  shell:
    """
    bpnet-gc-reference \
      --ref_fasta {input.fasta} --chrom_sizes {input.sizes} \
      --output_prefix {params.out_prefix} \
      --inputlen {params.input_len} --stride 1000
    
    mv {params.out_prefix}*.bed {output.gc_bed}
    """

rule gc_background:
  input:
    fasta = "resources/reference/genome.fa",
    peaks = "results/{sample}/data/peaks_inliers.bed",
    ref_gc = "resources/reference/genomewide_gc.bed"
  output:
    negatives = "results/{sample}/data/gc_negatives.bed"
  params:
    flank = config["bpnet_params"]["flank_size"],
    ratio = config["bpnet_params"]["neg_to_pos_ratio"],
    out_dir = "results/{sample}/data",
    out_prefix = "results/{sample}/data/gc_negatives"
  conda: "../envs/bpnet.yaml"
  shell:
    """
    bpnet-gc-background \
      --ref_fasta {input.fasta} --peaks_bed {input.peaks} \
      --out_dir {params.out_dir} \
      --ref_gc_bed {input.ref_gc} \
      --output_prefix {params.out_prefix} \
      --flank_size {params.flank} --neg_to_pos_ratio_train {params.ratio}
      
    mv {params.out_prefix}.bed {output.negatives}
    """
