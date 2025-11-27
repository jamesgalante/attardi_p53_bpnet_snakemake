rule prepare_outlier_json:
  input:
    plus = "results/{sample}/data/experiment_plus.bw",
    minus = "results/{sample}/data/experiment_minus.bw",
    ctl_plus = "results/{sample}/data/control_plus.bw",
    ctl_minus = "results/{sample}/data/control_minus.bw",
    peaks = "results/{sample}/data/idr_peaks.bed"
  output:
    json = "results/{sample}/configs/input_outliers.json"
  run:
    import json
    import os
    os.makedirs(os.path.dirname(output.json), exist_ok=True)
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
    input_len = config["bpnet_params"]["input_len"]
  conda: "../envs/bpnet.yaml"
  shell:
    """
    OUT="{output.gc_bed}"
    PREFIX=${{OUT%.bed}}

    bpnet-gc-reference \
      --ref_fasta {input.fasta} \
      --chrom_sizes {input.sizes} \
      --output_prefix $PREFIX \
      --inputlen {params.input_len} \
      --stride 1000
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
    ratio = config["bpnet_params"]["neg_to_pos_ratio"]
  conda: "../envs/bpnet.yaml"
  shell:
    """
    NEG="{output.negatives}"
    OUT_DIR=$(dirname "$NEG")
    
    bpnet-gc-background \
      --ref_fasta {input.fasta} \
      --peaks_bed {input.peaks} \
      --out_dir $OUT_DIR \
      --ref_gc_bed {input.ref_gc} \
      --output_prefix gc_negatives \
      --flank_size {params.flank} \
      --neg_to_pos_ratio_train {params.ratio}
    """
