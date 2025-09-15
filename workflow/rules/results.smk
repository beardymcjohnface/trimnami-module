rule trimnami_save_output:
    """Save the final trimmed output file"""
    input:
        os.path.join(config["trimnami"]["args"]["output_paths"]["temp"], "{file}"),
    output:
        os.path.join(config["trimnami"]["args"]["output_paths"]["results"], "{file}")
    localrule: True
    shell:
        "mv {input} {output}; "
        "ln -s $(pwd)/{output} $(pwd)/{input}"


rule trimnami_unzip_output:
    """if unzipped fastq or fasta are specified for some reason"""
    input:
        os.path.join(config["trimnami"]["args"]["output_paths"]["temp"], "{file}.{ext}.gz"),
    output:
        temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"], "{file}.{ext}"))
    wildcard_constraints:
        ext = r'fasta|fastq'
    localrule: True
    shell:
        "gzip -d -c {input} > {output}"


rule trimnami_init_input_paired_end:
    """Initialise the input files for paired end reads"""
    input:
        r1 = lambda wildcards: config["trimnami"]["samples"]["reads"][wildcards.sample]["R1"],
        r2 = lambda wildcards: config["trimnami"]["samples"]["reads"][wildcards.sample]["R2"],
    output:
        r1 = temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{sample}.R1.fastq.gz")),
        r2 = temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{sample}.R2.fastq.gz")),
        s = temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{sample}.RS.fastq.gz")),
    params:
        s = lambda wildcards: config["trimnami"]["samples"]["reads"][wildcards.sample]["S"],
        is_paired = True
    resources:
        **config["resources"]["med"]
    threads:
        config["resources"]["med"]["cpu"]
    conda:
        os.path.join("..", "envs", "seqtk.yaml")
    script:
        os.path.join("..", "scripts", "copyOrGzip.py")



rule trimnami_init_input_single_end:
    """Initialise the input files for paired end reads"""
    input:
        r1=lambda wildcards: config["trimnami"]["samples"]["reads"][wildcards.sample]["R1"],
    output:
        r1=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{sample}.S.fastq.gz")),
    params:
        is_paired = False
    conda:
        os.path.join("..","envs","seqtk.yaml")
    resources:
        **config["resources"]["med"]
    threads:
        config["resources"]["med"]["cpu"]
    script:
        os.path.join("..", "scripts", "copyOrGzip.py")