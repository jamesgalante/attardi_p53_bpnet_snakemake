rule run_tfmodisco:
  input:
    shap_scores = "results/{sample}/shap/{head}_scores.h5"
  output:
    modisco_results = "results/{sample}/modisco/{head}/modisco_results.h5"
  conda: "../envs/tfmodisco.yaml"
  resources:
    slurm_partition = "akundaje",
    mem = "16G",
    runtime = "4h"
  threads: 8
  shell:
    """
    modisco motifs \
      --max_seqlets 2000 \
      --h5py {input.shap_scores} \
      -o {output.modisco_results} \
      --trim_size 20 \
      --initial_flank_to_add 5 \
      --final_flank_to_add 10 \
      --verbose
    """

rule generate_modisco_reports:
  input:
    counts = "results/{sample}/modisco/counts/modisco_results.h5",
    profile = "results/{sample}/modisco/profile/modisco_results.h5",
    meme_db = "resources/reference/motifs.txt"
  output:
    counts_report = "results/{sample}/modisco/counts/report.html",
    profile_report = "results/{sample}/modisco/profile/report.html"
  conda: "../envs/tfmodisco.yaml"
  resources:
    slurm_partition = "akundaje",
    mem = "8G",
    runtime = "1h"
  shell:
    """
    modisco report \
      -i {input.counts} \
      -o $(dirname {output.counts_report}) \
      -s $(dirname {output.counts_report}) \
      -m {input.meme_db}

    modisco report \
      -i {input.profile} \
      -o $(dirname {output.profile_report}) \
      -s $(dirname {output.profile_report}) \
      -m {input.meme_db}
    """
