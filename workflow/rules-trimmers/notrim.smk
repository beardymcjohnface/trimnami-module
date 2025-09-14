
rule trimnami_notrim_paired_end:
    """Skip read trimming for paired reads"""
    input:
        r1=os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.R1.fastq.gz"),
        r2=os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.R2.fastq.gz"),
        s=os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.RS.fastq.gz"),
    output:
        r1=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.notrim.R1.fastq.gz")),
        r2=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.notrim.R2.fastq.gz")),
        s=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.notrim.RS.fastq.gz")),
    params:
        is_paired = True
    shell:
        """
        mv {input.r1} {output.r1}
        mv {input.r2} {output.r2}
        mv {input.s} {output.s}
        ln -s $(pwd)/{output.r1} $(pwd)/{input.r1}
        ln -s $(pwd)/{output.r2} $(pwd)/{input.r2}
        ln -s $(pwd)/{output.s} $(pwd)/{input.s}
        """


rule trimnami_notrim_single_end:
    """Skip read trimming for single end"""
    input:
        r1=os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.S.fastq.gz"),
    output:
        r1=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.notrim.S.fastq.gz")),
    params:
        is_paired = False
    shell:
        """
        mv {input.r1} {output.r1}
        ln -s $(pwd)/{output.r1} $(pwd)/{input.r1}
        """
