rule download_and_prep_references:
  output:
    fasta = "resources/reference/genome.fa",
    fai = "resources/reference/genome.fa.fai",
    sizes = "resources/reference/chrom.sizes",
    chroms_txt = "resources/reference/chroms.txt",
    blacklist = "resources/reference/blacklist.bed"
  params:
    fasta_url = config["reference"]["fasta_url"],
    sizes_url = config["reference"]["chrom_sizes_url"],
    blacklist_url = config["reference"]["blacklist_url"]
  resources:
    mem = "16G",
    time = "2:00:00"
  shell:
    """
    ml load biology samtools

    wget --no-check-certificate {params.fasta_url} -O {output.fasta}.gz
    gunzip -f {output.fasta}.gz
    samtools faidx {output.fasta}

	wget --no-check-certificate -O - {params.sizes_url} | \
	grep -v -e '_' -e 'chrEBV' > {output.sizes}
    
    awk '{{print $1}}' {output.sizes} > {output.chroms_txt}

    wget --no-check-certificate {params.blacklist_url} -O {output.blacklist}.gz
    gunzip -f {output.blacklist}.gz
    """

# See https://bowtie-bio.sourceforge.net/bowtie2/faq.shtml - L=36
rule build_bowtie_index:
  input:
    fasta = "resources/reference/genome.fa"
  output:
    marker = "resources/reference/index/build.done" 
  threads: 8
  resources:
    mem = "16G",
    time = "2:00:00"
  shell:
    """
    ml load biology bowtie/1.2.2
    
    mkdir -p resources/reference/index
    
    # Using bowtie-build for v1 index creation
    bowtie-build --threads {threads} {input.fasta} resources/reference/index/genome && \
    touch {output.marker}
    """

rule download_motifs:
  output: "resources/reference/motifs.txt"
  params: config['reference']['motif_url']
  shell: "wget -O {output} {params}"
