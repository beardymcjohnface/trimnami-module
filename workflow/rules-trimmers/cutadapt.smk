
rule trimnami_cutadapt_paired_end:
    """Skip read trimming for paired reads"""
    input:
        r1=os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.R1.fastq.gz"),
        r2=os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.R2.fastq.gz"),
        s=os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.RS.fastq.gz"),
        adapters=workflow.source_path(
            os.path.join("..", "..", "resources",config["trimnami"]["qc"]["cutadapt"]["adapters"])
        )
    output:
        r1=os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.cutadapt.R1.fastq.gz"),
        r2=os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.cutadapt.R2.fastq.gz"),
        s=os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.cutadapt.RS.fastq.gz"),
    resources:
        **config["resources"]["med"]
    threads:
        config["resources"]["med"]["cpu"]
    params:
        config["trimnami"]["qc"]["cutadapt"]["params"]
    conda:
        os.path.join("..", "envs","cutadapt.yaml")
    benchmark:
        os.path.join(config["trimnami"]["args"]["output_paths"]["bench"],"trimnami_cutadapt_paired_end.{file}.txt")
    log:
        os.path.join(config["trimnami"]["args"]["output_paths"]["log"],"trimnami_cutadapt_paired_end.{file}.log")
    shell:
        ("cutadapt "
            "--cores {threads} "
            "{params} "
            "-b file:{input.adapters} "
            "-B file:{input.adapters} "
            "-o {output.r1} "
            "-p {output.r2} "
            "--fasta "
            "{input.r1} "
            "{input.r2} "
            "&> {log}\n\n "
        "if [[ -s {input.s} ]]\n then "
            "cutadapt "
                "--cores {threads} "
                "{params} "
                "-b {input.adapters} "
                "-o {output.s} "
                "{input.s} "
                "&> {log}\n "
        "else "
            "touch {output.s}\n "
        "fi\n\n ")


rule trimnami_cutadapt_single_end:
    """Skip read trimming for single end"""
    input:
        r1=os.path.join(config["trimnami"]["args"]["output_paths"]["temp"], "{file}.S.fastq.gz"),
        adapters=workflow.source_path(
            os.path.join("..","..","resources",config["trimnami"]["qc"]["cutadapt"]["adapters"])
        )
    output:
        r1=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.cutadapt.S.fastq.gz")),
    resources:
        **config["resources"]["med"]
    threads:
        config["resources"]["med"]["cpu"]
    params:
        config["trimnami"]["qc"]["cutadapt"]
    conda:
        os.path.join("..", "envs","cutadapt.yaml")
    benchmark:
        os.path.join(config["trimnami"]["args"]["output_paths"]["bench"],"trimnami_cutadapt_single_end.{file}.txt")
    log:
        os.path.join(config["trimnami"]["args"]["output_paths"]["log"],"trimnami_cutadapt_single_end.{file}.log")
    shell:
        ("cutadapt "
            "--cores {threads} "
            "{params} "
            "-b file:{input.adapters} "
            "-o {output.r1} "
            "--fasta  "
            "{input.r1} "
            "&> {log}\n\n ")


# rule fasta_to_fastq:
#     """Convert the fasta files to fastq files for cutadapt"""
#     input:
#         os.path.join(dir["cutadapt"], "{file}.fasta")
#     output:
#         temp(os.path.join(dir["cutadapt"],"{file}.fastq.gz"))
#     params:
#         compression = "-" + str(config["qc"]["compression"])
#     conda:
#         os.path.join(dir["env"],"seqtk.yaml")
#     benchmark:
#         os.path.join(dir["bench"],"fasta_to_fastq.{file}.txt")
#     log:
#         os.path.join(dir["log"],"fasta_to_fastq.{file}.log")
#     shell:
#         "seqtk "
#             "seq -F 'B' {input} "
#             "| gzip {params.compression} "
#             "> {output} "
#             "2> {log}\n\n "