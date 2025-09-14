
rule trimnami_filtlong_single:
    input:
        i=os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.S.fastq.gz")
    output:
        o=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.filtlong.S.fastq.gz")),
    resources:
        **config["resources"]["med"]
    threads:
        config["resources"]["med"]["cpu"]
    conda:
        os.path.join("..", "envs", "filtlong.yaml")
    params:
        config["trimnami"]["qc"]["filtlong"]
    benchmark:
        os.path.join(config["trimnami"]["args"]["output_paths"]["bench"],"trimnami_filtlong_single.{file}.txt")
    log:
        os.path.join(config["trimnami"]["args"]["output_paths"]["log"], "trimnami_filtlong_single.{file}.log")
    shell:
        "filtlong {params} {input.i} 2> {log} "
            "| gzip -1 "
            "> {output.o}; "


rule trimnami_filtlong_paried:
    """
    You probably don't want to be running filtlong on paired reads. The current rule will also mess up the read
    pairing.
    """
    input:
        r1=os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.R1.fastq.gz"),
        r2=os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.R2.fastq.gz"),
        s=os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.RS.fastq.gz"),
    output:
        r1=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.filtlong.R1.fastq.gz")),
        r2=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.filtlong.R2.fastq.gz")),
        s=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.filtlong.RS.fastq.gz")),
    params:
        config["trimnami"]["qc"]["filtlong"]
    resources:
        **config["resources"]["med"]
    threads:
        config["resources"]["med"]["cpu"]
    conda:
        os.path.join("..", "envs", "filtlong.yaml")
    benchmark:
        os.path.join(config["trimnami"]["args"]["output_paths"]["bench"],"trimnami_filtlong_paried.{file}.txt")
    log:
        os.path.join(config["trimnami"]["args"]["output_paths"]["log"], "trimnami_filtlong_paried.{file}.log")
    shell:
        "filtlong {params} {input.r1} 2> {log}"
            "| gzip -1 "
            "> {output.r1}\n\n "
        "filtlong {params} {input.r2} 2> {log}"
            "| gzip -1 "
            "> {output.r2}\n\n "
        "if [[ -s {input.s} ]]\n "
        "then "
            "filtlong {params} {input.s} 2> {log}"
                "| gzip -1 "
                "> {output.s}\n\n "
        "else "
            "touch {output.s}\n "
        "fi\n\n "
