rule trimnami_fastq_to_fasta:
    """Convert the trimmed fastq file to a fasta file"""
    input:
        os.path.join(config["trimnami"]["args"]["output_paths"]["temp"], "{file}.fastq.gz")
    output:
        temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"], "{file}.fasta.gz")),
    params:
        "-" + str(config["trimnami"]["qc"]["compression"])
    resources:
        **config["resources"]["med"]
    threads:
        config["resources"]["med"]["cpu"]
    conda:
        os.path.join("..", "envs", "seqtk.yaml")
    shell:
        "seqtk seq {input} -A "
            "| gzip {params} "
            "> {output}\n\n "
