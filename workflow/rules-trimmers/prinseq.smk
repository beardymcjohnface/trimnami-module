
rule trimnami_prinseq_paired:
    input:
        r1=os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.R1.fastq.gz"),
        r2=os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.R2.fastq.gz"),
        s=os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.RS.fastq.gz"),
    output:
        r1=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.prinseq.R1.fastq.gz")),
        r2=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.prinseq.R2.fastq.gz")),
        s=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.prinseq.RS.fastq.gz")),
        s1=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.prinseq.S1.fastq.gz")),
        s2=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.prinseq.S2.fastq.gz")),
    resources:
        **config["resources"]["med"]
    threads:
        config["resources"]["med"]["cpu"]
    conda:
        os.path.join("..", "envs", "prinseq.yaml")
    params:
        config["trimnami"]["qc"]["prinseq"],
    log:
        os.path.join(config["trimnami"]["args"]["output_paths"]["log"], "trimnami_prinseq_paired.{file}.log")
    benchmark:
        os.path.join(config["trimnami"]["args"]["output_paths"]["bench"],"trimnami_prinseq_paired.{file}.txt")
    shell:
        "prinseq++ {params} "
            "-out_gz "
            "-threads {threads} "
            "-out_good {output.r1} "
            "-out_good2 {output.r2} "
            "-out_single {output.s1} "
            "-out_single2 {output.s2} "
            "-out_bad /dev/null "
            "-out_bad2 /dev/null "
            "-fastq {input.r1} "
            "-fastq2 {input.r2}  &> {log}\n\n "
        "if [[ -s {input.s} ]]\n "
        "then "
            "prinseq++ {params} "
                "-out_gz "
                "-threads {threads} "
                "-out_good {output.s} "
                "-out_bad /dev/null "
                "-fastq {input.r1} &> {log}\n\n "
        "else "
            "touch {output.s}\n "
        "fi\n\n "
        "cat {output.s1} {output.s2} >> {output.s}\n\n "


rule prinseq_single:
    input:
        r1=os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.S.fastq.gz"),
    output:
        r1=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.prinseq.S.fastq.gz")),
    resources:
        **config["resources"]["med"]
    threads:
        config["resources"]["med"]["cpu"]
    conda:
        os.path.join("..", "envs","prinseq.yaml")
    params:
        config["trimnami"]["qc"]["prinseq"]
    log:
        os.path.join(config["trimnami"]["args"]["output_paths"]["log"], "prinseq.{file}.log")
    benchmark:
        os.path.join(config["trimnami"]["args"]["output_paths"]["bench"],"prinseq.{file}.txt")
    shell:
        "prinseq++ {params} "
            "-out_gz "
            "-threads {threads} "
            "-out_good {output.r1} "
            "-out_bad /dev/null "
            "-fastq {input.r1} &> {log}\n\n "
