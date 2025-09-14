rule trimnami_fastqc:
    input:
        os.path.join(config["trimnami"]["args"]["output_paths"]["temp"], "{file}.fastq.gz"),
    output:
        temp(os.path.join(config["trimnami"]["args"]["output_paths"]["reports"], "{file}_fastqc.zip")),
        temp(os.path.join(config["trimnami"]["args"]["output_paths"]["reports"], "{file}_fastqc.html"))
    params:
        config["trimnami"]["args"]["output_paths"]["reports"],
    conda:
        os.path.join("..", "envs", "fastqc.yaml")
    resources:
        **config["resources"]["med"]
    threads:
        config["resources"]["med"]["cpu"]
    log:
        os.path.join(config["trimnami"]["args"]["output_paths"]["log"], "trimnami_fastqc.{file}.log")
    benchmark:
        os.path.join(config["trimnami"]["args"]["output_paths"]["bench"], "trimnami_fastqc.{file}.txt")
    shell:
        ("fastqc {input} "
            "-t {threads} "
            "--outdir {params} "
            "&> {log}; ")


rule trimnami_multiqc_fastqc:
    input:
        config["trimnami"]["multiqc"]
    output:
        os.path.join(config["trimnami"]["args"]["output_paths"]["results"], "multiqc.html")
    params:
        os.path.join(config["trimnami"]["args"]["output_paths"]["reports"])
    conda:
        os.path.join("..", "envs","multiqc.yaml")
    resources:
        **config["resources"]["med"]
    threads:
        config["resources"]["med"]["cpu"]
    log:
        os.path.join(config["trimnami"]["args"]["output_paths"]["log"],"trimnami_multiqc_fastqc.log")
    benchmark:
        os.path.join(config["trimnami"]["args"]["output_paths"]["bench"],"trimnami_multiqc_fastqc.txt")
    shell:
        "multiqc {params} -n {output} --no-data-dir 2> {log}"
