rule generate_training_configs:
  input:
    plus = "results/{sample}/data/experiment_plus.bw",
    minus = "results/{sample}/data/experiment_minus.bw",
    ctl_plus = "results/{sample}/data/control_plus.bw",
    ctl_minus = "results/{sample}/data/control_minus.bw",
    peaks = "results/{sample}/data/peaks_inliers.bed",
    negatives = "results/{sample}/data/gc_negatives.bed"
  output:
    input_json = "results/{sample}/configs/input_data.json",
    splits_json = "results/{sample}/configs/splits.json",
    params_json = "results/{sample}/configs/bpnet_params.json"
  run:
    import json
    p = config["bpnet_params"]
    
    ratio_decimal = 1.0 / float(p["neg_to_pos_ratio"])
    input_data = {
      "0": {
        "signal": {"source": [input.plus, input.minus]},
        "loci": {"source": [input.peaks]},
        "background_loci": {"source": [input.negatives], "ratio": [ratio_decimal]},
        "bias": {
          "source": [input.ctl_plus, input.ctl_minus],
          "smoothing": [None, None]
        }
      }
    }
    
    splits_data = {"0": p["splits"]}
    
    model_params = {
      "input_len": p["input_len"],
      "output_profile_len": p["output_len"],
      "motif_module_params": {"filters": [p["filters"]], "kernel_sizes": [p["conv_kernel_size"]], "padding": "valid"},
      "syntax_module_params": {"num_dilation_layers": p["num_dilation_layers"], "filters": p["filters"], "kernel_size": 3, "padding": "valid", "pre_activation_residual_unit": True},
      "profile_head_params": {"filters": 1, "kernel_size": p["profile_kernel_size"], "padding": "valid"},
      "counts_head_params": {"units": [1], "dropouts": [0.0], "activations": ["linear"]},
      "loss_weights": p["loss_weights"],
      "counts_loss": "MSE"
    }

    with open(output.input_json, 'w') as f: json.dump(input_data, f, indent=4)
    with open(output.splits_json, 'w') as f: json.dump(splits_data, f, indent=4)
    with open(output.params_json, 'w') as f: json.dump(model_params, f, indent=4)

rule train_bpnet:
  input:
    input_json = "results/{sample}/configs/input_data.json",
    splits_json = "results/{sample}/configs/splits.json",
    params_json = "results/{sample}/configs/bpnet_params.json",
    peaks = "results/{sample}/data/peaks_inliers.bed"
  output:
    model_dir = directory("results/{sample}/bpnet_model/models"),
    flag = touch("results/{sample}/bpnet_model/training.done")
  params:
    genome = "resources/reference/genome.fa",
    sizes = "resources/reference/chrom.sizes",
    chroms = lambda w: " ".join(config["bpnet_params"]["splits"]["train"] + config["bpnet_params"]["splits"]["val"] + config["bpnet_params"]["splits"]["test"]),
    input_len = config["bpnet_params"]["input_len"],
    output_len = config["bpnet_params"]["output_len"]
  conda: "../envs/bpnet.yaml"
  resources: gpu=1
  shell:
    """
    mkdir -p {output.model_dir}
    bpnet-train \
      --input-data {input.input_json} --output-dir {output.model_dir} \
      --reference-genome {params.genome} --chrom-sizes {params.sizes} \
      --chroms {params.chroms} --splits {input.splits_json} \
      --model-arch-name BPNet --model-arch-params-json {input.params_json} \
      --sequence-generator-name BPNet \
      --input-seq-len {params.input_len} --output-len {params.output_len} \
      --shuffle --epochs 100 --early-stopping-patience 10
    """

rule predict:
  input:
    flag = "results/{sample}/bpnet_model/training.done",
    input_json = "results/{sample}/configs/input_data.json"
  output:
    metrics = "results/{sample}/bpnet_model/evaluation/predictions_metrics.json"
  params:
    model_dir = "results/{sample}/bpnet_model/models",
    sizes = "resources/reference/chrom.sizes",
    genome = "resources/reference/genome.fa",
    test_chroms = " ".join(config["bpnet_params"]["splits"]["test"]),
    out_dir = "results/{sample}/bpnet_model/evaluation",
    input_len = config["bpnet_params"]["input_len"],
    output_len = config["bpnet_params"]["output_len"]
  conda: "../envs/bpnet.yaml"
  shell:
    """
    MODEL_PATH=$(find {params.model_dir} -name "model_split*" | head -n 1)
    
    bpnet-predict \
      --model $MODEL_PATH --chrom-sizes {params.sizes} \
      --chroms {params.test_chroms} --test-indices-file None \
      --reference-genome {params.genome} --output-dir {params.out_dir} \
      --input-data {input.input_json} --sequence-generator-name BPNet \
      --input-seq-len {params.input_len} --output-len {params.output_len} \
      --output-window-size {params.output_len} --batch-size 64 \
      --generate-predicted-profile-bigWigs
    """
