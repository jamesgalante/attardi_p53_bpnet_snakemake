rule predict_bpnet_denoise:
  input:
    training_done = "results/{sample}/models/training.done",
    input_json = "results/{sample}/configs/input_data.json",
    genome = "resources/reference/genome.fa",
    chroms = "resources/reference/chroms.txt",
    sizes = "resources/reference/chrom.sizes"
  output:
    preds_dir = directory("results/{sample}/predictions_denoise"),
    flag = touch("results/{sample}/predictions_denoise/denoise.done")
  params:
    model_path = "results/{sample}/models/model_split000",
    input_len = config["bpnet_params"]["input_len"],
    output_len = config["bpnet_params"]["output_len"]
  conda: "../envs/bpnet.yaml"
  resources:
    slurm_partition = "akundaje",
    gpu = 1,
    tasks_per_gpu = 0,
    runtime = "6h"
  threads: 1
  shell:
    """
    mkdir -p {output.preds_dir}

    bpnet-predict \
      --model {params.model_path} \
      --chrom-sizes {input.sizes} \
      --chroms $(paste -s -d ' ' {input.chroms}) \
      --test-indices-file None \
      --reference-genome {input.genome} \
      --output-dir {output.preds_dir} \
      --input-data {input.input_json} \
      --sequence-generator-name BPNet \
      --input-seq-len {params.input_len} \
      --output-len {params.output_len} \
      --output-window-size {params.output_len} \
      --batch-size 64 \
      --reverse-complement-average \
      --threads {threads} \
      --generate-predicted-profile-bigWigs
    """

rule compute_shap:
  input:
    training_done = "results/{sample}/models/training.done",
    input_json = "results/{sample}/configs/input_data.json",
    peaks = "results/{sample}/data/peaks_inliers.bed",
    genome = "resources/reference/genome.fa",
    sizes = "resources/reference/chrom.sizes"
  output:
    shap_dir = directory("results/{sample}/shap"),
    flag = touch("results/{sample}/shap/shap.done"),
    counts_scores = "results/{sample}/shap/counts_scores.h5",
    profile_scores = "results/{sample}/shap/profile_scores.h5"
  params:
    model_path = "results/{sample}/models/model_split000",
    input_len = config["bpnet_params"]["input_len"],
    output_len = config["bpnet_params"]["output_len"]
  conda: "../envs/bpnet.yaml"
  resources:
    slurm_partition = "akundaje",
    gpu = 1,
    tasks_per_gpu = 0,
    runtime = "12h"
  threads: 1
  shell:
    """
    mkdir -p {output.shap_dir}

    bpnet-shap \
      --reference-genome {input.genome} \
      --model {params.model_path} \
      --bed-file {input.peaks} \
      --output-dir {output.shap_dir} \
      --input-seq-len {params.input_len} \
      --control-len {params.output_len} \
      --task-id 0 \
      --input-data {input.input_json} \
      --chrom-sizes {input.sizes} \
      --generate-shap-bigWigs
    """
