
rule generate_input_data_config:
  input:
    plus = "results/{sample}/data/experiment_plus.bw",
    minus = "results/{sample}/data/experiment_minus.bw",
    ctl_plus = "results/{sample}/data/control_plus.bw",
    ctl_minus = "results/{sample}/data/control_minus.bw",
    peaks = "results/{sample}/data/peaks_inliers.bed",
    negatives = "results/{sample}/data/gc_negatives.bed"
  output:
    input_json = "results/{sample}/configs/input_data.json"
  run:
    import json
    ratio = float(config["bpnet_params"]["ratio"])
    input_data = {
      "0": {
        "signal": {
			"source": [input.plus, input.minus]
		},
        "loci": {
			"source": [input.peaks]
		},
        "background_loci": {
			"source": [input.negatives], "ratio": [ratio]
		},
        "bias": {
          "source": [input.ctl_plus, input.ctl_minus],
          "smoothing": [None, None]
        }
      }
    }
    with open(output.input_json, 'w') as f: json.dump(input_data, f, indent=4)

rule compute_counts_weight:
  input:
    input_json = "results/{sample}/configs/input_data.json"
  output:
    weight_file = "results/{sample}/configs/counts_weight.txt"
  conda: "../envs/bpnet.yaml"
  shell:
    """
    bpnet-counts-loss-weight --input-data {input.input_json} > {output.weight_file}
    """

rule generate_training_configs:
  input:
    plus = "results/{sample}/data/experiment_plus.bw",
    minus = "results/{sample}/data/experiment_minus.bw",
    ctl_plus = "results/{sample}/data/control_plus.bw",
    ctl_minus = "results/{sample}/data/control_minus.bw",
    peaks = "results/{sample}/data/peaks_inliers.bed",
    negatives = "results/{sample}/data/gc_negatives.bed",
    weight_file = "results/{sample}/configs/counts_weight.txt" 
  output:
    splits_json = "results/{sample}/configs/splits.json",
    params_json = "results/{sample}/configs/bpnet_params.json"
  run:
    import json
    
    with open(input.weight_file, 'r') as f:
      counts_weight = float(f.read().strip())

    p = config["bpnet_params"]
    model_params = {
      "input_len": p["input_len"],
      "output_profile_len": p["output_len"],
      "motif_module_params": {"filters": [p["filters"]], "kernel_sizes": [p["conv_kernel_size"]], "padding": "valid"},
      "syntax_module_params": {"num_dilation_layers": p["num_dilation_layers"], "filters": p["filters"], "kernel_size": 3, "padding": "valid", "pre_activation_residual_unit": True},
      "profile_head_params": {"filters": 1, "kernel_size": p["profile_kernel_size"], "padding": "valid"},
      "counts_head_params": {"units": [1], "dropouts": [0.0], "activations": ["linear"]},
      "profile_bias_module_params": {"kernel_sizes": [1]},
      "counts_bias_module_params": {},
      "use_attribution_prior": False,
	  "attribution_prior_params": {"frequency_limit": 150, "limit_softness": 0.2, "grad_smooth_sigma": 3, "profile_grad_loss_weight": 200, "counts_grad_loss_weight": 100},
      "loss_weights": [1, counts_weight],
      "counts_loss": "MSE"
    }

    splits_data = {"0": p["splits"]}

    with open(output.splits_json, 'w') as f: json.dump(splits_data, f, indent=4)
    with open(output.params_json, 'w') as f: json.dump(model_params, f, indent=4)

rule train_bpnet:
  input:
    input_json = "results/{sample}/configs/input_data.json",
    splits_json = "results/{sample}/configs/splits.json",
    params_json = "results/{sample}/configs/bpnet_params.json",
    peaks = "results/{sample}/data/peaks_inliers.bed",
    genome = "resources/reference/genome.fa",
    chroms = "resources/reference/chroms.txt",
    sizes = "resources/reference/chrom.sizes"
  output:
    model_dir = directory("results/{sample}/models"),
    flag = touch("results/{sample}/models/training.done")
  params:
    input_len = config["bpnet_params"]["input_len"],
    output_len = config["bpnet_params"]["output_len"]
  conda: "../envs/bpnet.yaml"
  threads: 1
  resources:
    slurm_partition = "akundaje",
    gpu = 1,
    tasks_per_gpu = 0,
    runtime = "6h"
  shell:
    """
    mkdir -p {output.model_dir}
    bpnet-train \
      --input-data {input.input_json} \
      --output-dir {output.model_dir} \
      --reference-genome {input.genome} \
      --chroms $(paste -s -d ' ' {input.chroms}) \
      --chrom-sizes {input.sizes} \
      --splits {input.splits_json} \
      --model-arch-name BPNet \
      --model-arch-params-json {input.params_json} \
      --sequence-generator-name BPNet \
      --model-output-filename model \
      --input-seq-len {params.input_len} \
      --output-len {params.output_len} \
      --shuffle \
      --threads {threads} \
      --epochs 100 \
      --batch-size 64 \
      --reverse-complement-augmentation \
      --early-stopping-patience 10 \
      --reduce-lr-on-plateau-patience 5 \
      --learning-rate 0.001
    """

rule predict_bpnet_test:
  input:
    training_done = "results/{sample}/models/training.done",
    input_json = "results/{sample}/configs/input_data.json",
    genome = "resources/reference/genome.fa",
    sizes = "resources/reference/chrom.sizes"
  output:
    preds_dir = directory("results/{sample}/predictions_test"),
    flag = touch("results/{sample}/predictions_test/predict.done")
  params:
    model_path = "results/{sample}/models/model_split000",
    input_len = config["bpnet_params"]["input_len"],
    output_len = config["bpnet_params"]["output_len"],
    test_chroms = " ".join(config["bpnet_params"]["splits"]["test"])
  conda: "../envs/bpnet.yaml"
  resources:
    slurm_partition = "akundaje",
    gpu = 1,
    tasks_per_gpu = 0,
    runtime = "3h"
  threads: 1
  shell:
    """
    mkdir -p {output.preds_dir}
    
    bpnet-predict \
      --model {params.model_path} \
      --chrom-sizes {input.sizes} \
      --chroms {params.test_chroms} \
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
