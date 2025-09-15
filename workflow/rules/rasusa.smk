rule trimnami_rasusa:
    """Note, this will mess up paired reads"""
    input:
        os.path.join(config["trimnami"]["args"]["output_paths"]["temp"], "{file}.{ext}.fastq.gz")
    output:
        temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"], "{file}.subsample.{ext}.fastq.gz")),
    resources:
        **config["resources"]["med"]
    threads:
        config["resources"]["med"]["cpu"]
    conda:
        os.path.join("..", "envs", "rasusa.yaml")
    params:
        config["trimnami"]["qc"]["subsample"]
    benchmark:
        os.path.join(config["trimnami"]["args"]["output_paths"]["bench"], "trimnami_rasusa.{file}.{ext}.txt")
    log:
        os.path.join(config["trimnami"]["args"]["output_paths"]["log"], "trimnami_rasusa.{file}.{ext}.log")
    shell:
        ("if (( $(wc -c {input} | awk '{{print$1}}') > 200 ))\n then "
            "rasusa reads "
                "-o {output} "
                "-O g "
                "{params} "
                "{input} "
                "2> {log}\n "
         "else "
            "touch {output}\n "
         "fi\n\n ")
