import os

rule trimnami_sample_manifest:
    """Print the parsed sample manifest"""
    output:
        tsv = os.path.join(config["trimnami"]["args"]["output"],"samples.tsv")
    params:
        sample_dict = config["trimnami"]["samples"]
    localrule: True
    run:
        from metasnek import fastq_finder
        fastq_finder.write_samples_tsv(params.sample_dict,output.tsv)


rule trimnami_build_env:
    output:
        os.path.join(config["trimnami"]["args"]["output_paths"]["temp"], "{env}.done")
    conda:
        lambda wildcards: os.path.join("..", "envs", wildcards.env)
    shell:
        "touch {output}"
