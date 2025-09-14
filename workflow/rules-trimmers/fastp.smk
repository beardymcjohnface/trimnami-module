
rule trimnami_fastp_paired_end:
    """Read trimming with fastp for paired reads"""
    input:
        r1=os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.R1.fastq.gz"),
        r2=os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.R2.fastq.gz"),
        s=os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.RS.fastq.gz"),
    output:
        r1=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.fastp.R1.fastq.gz")),
        r2=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.fastp.R2.fastq.gz")),
        s=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.fastp.RS.fastq.gz")),
        s1=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.fastp.S1.fastq.gz")),
        s2=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.fastp.S2.fastq.gz")),
        stats=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.fastp.stats.json")),
        html=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.fastp.stats.html"))
    benchmark:
        os.path.join(config["trimnami"]["args"]["output_paths"]["bench"],"trimnami_fastp_paired_end.{file}.txt")
    log:
        os.path.join(config["trimnami"]["args"]["output_paths"]["log"],"trimnami_fastp_paired_end.{file}.log")
    resources:
        mem_mb=resources["med"]["mem"],
        mem=str(resources["med"]["mem"]) + "MB",
        time=resources["med"]["time"]
    threads:
        resources["med"]["cpu"]
    conda:
        os.path.join("..", "envs","fastp.yaml")
    params:
        fastp=config["trimnami"]["qc"]["fastp"],
        compression=config["trimnami"]["qc"]["compression"],
    shell:
        "fastp "
            "-i {input.r1} "
            "-I {input.r2} "
            "-o {output.r1} "
            "-O {output.r2} "
            "--unpaired1 {output.s1} "
            "--unpaired2 {output.s2} "
            "-z {params.compression} "
            "-j {output.stats} "
            "-h {output.html} "
            "--thread {threads} "
            "{params.fastp} "
            "2> {log}\n\n "
        "if [[ -s {input.s} ]]\n "
        "then "
            "fastp "
            "-i {input.s} "
            "-o {output.s} "
            "-z {params.compression} "
            "-j {output.stats} "
            "-h {output.html} "
            "--thread {threads} "
            "{params.fastp} "
            "2> {log}\n "
        "else "
            "touch {output.s}\n "
        "fi\n\n "
        "cat {output.s1} {output.s2} >> {output.s}\n\n "


rule trimnami_fastp_single_end:
    """Read trimming with fastp for single end reads"""
    input:
        r1=os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.S.fastq.gz"),
    output:
        r1=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.fastp.S.fastq.gz")),
        stats=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.fastp.S.stats.json")),
        html=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.fastp.S.stats.html"))
    benchmark:
        os.path.join(config["trimnami"]["args"]["output_paths"]["bench"],"trimnami_fastp_single_end.{file}.txt")
    log:
        os.path.join(config["trimnami"]["args"]["output_paths"]["log"],"trimnami_fastp_single_end.{file}.log")
    resources:
        mem_mb=resources["med"]["mem"],
        mem=str(resources["med"]["mem"]) + "MB",
        time=resources["med"]["time"]
    threads:
        resources["med"]["cpu"]
    conda:
        os.path.join("..", "envs","fastp.yaml")
    params:
        fastp=config["trimnami"]["qc"]["fastp"],
        compression=config["trimnami"]["qc"]["compression"]
    shell:
        "fastp "
            "-i {input.r1} "
            "-o {output.r1} "
            "-z {params.compression} "
            "-j {output.stats} "
            "-h {output.html} "
            "--thread {threads} "
            "{params.fastp} "
            "2> {log}\n\n "
