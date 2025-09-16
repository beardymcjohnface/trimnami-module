import os

from metasnek import fastq_finder


"""
Parse the samples with metasnek
    - also initialise a dictionary for the final output reads files
"""
config["trimnami"]["samples"] = dict()
config["trimnami"]["trimmed"] = dict()
config["trimnami"]["samples"]["reads"] = fastq_finder.parse_samples_to_dictionary(config["trimnami"]["args"]["reads"])
config["trimnami"]["samples"]["names"] = list(config["trimnami"]["samples"]["reads"].keys())


"""
Expand output dir and file paths to include base output dir
"""
for output_file_path in config["trimnami"]["args"]["output_paths"]:
    config["trimnami"]["args"]["output_paths"][output_file_path] =  os.path.join(
        config["trimnami"]["args"]["output"], config["trimnami"]["args"]["output_paths"][output_file_path]
    )


"""
Define target reads files based on steps and samples
"""
config["trimnami"]["targets"] = dict()
config["trimnami"]["targets"]["reads"] = list()

for sample in config["trimnami"]["samples"]["names"]:
    config["trimnami"]["trimmed"][sample] = dict()
    if config["trimnami"]["samples"]["reads"][sample]["R2"] is not None:
        read_pairs = ["R1","R2","RS"]
    else:
        read_pairs = ["S"]
    for read_pair in read_pairs:
        out_file_name = os.path.join(
                config["trimnami"]["args"]["output_paths"]["results"],
                ".".join([
                    sample,
                    ".".join(config["trimnami"]["args"]["steps"]),
                    read_pair,
                    config["trimnami"]["args"]["outfmt"]
                ])
        )
        config["trimnami"]["targets"]["reads"].append(out_file_name)
        config["trimnami"]["trimmed"][sample][read_pair] = out_file_name


"""
Define target fastqc files and inputs based on steps and samples
"""
if config["trimnami"]["args"]["fastqc"]:
    config["trimnami"]["targets"]["reads"].append(
        os.path.join(config["trimnami"]["args"]["output_paths"]["results"], "multiqc.html")
    )

config["trimnami"]["multiqc"] = list()

trimming_stage = []
for operation in config["trimnami"]["args"]["steps"]:
    trimming_stage.append(operation)
    trimming_stage_name = ".".join(trimming_stage)
    for sample in config["trimnami"]["samples"]["names"]:
        if config["trimnami"]["samples"]["reads"][sample]["R2"] is not None:
            read_pairs = ["R1", "R2", "RS"]
        else:
            read_pairs = ["S"]
        for read_pair in read_pairs:
            config["trimnami"]["multiqc"].append(
                os.path.join(
                    config["trimnami"]["args"]["output_paths"]["reports"],
                    ".".join([
                        sample,
                        trimming_stage_name,
                        read_pair
                    ]) + "_fastqc.zip"
                )
            )

"""
Add targets for pre-building the environments
"""
config["trimnami"]["targets"]["envs"] = list()

for filename in os.listdir(os.path.join(workflow.basedir, "envs")):
    if filename.endswith(".yaml") or filename.endswith(".yml"):
        config["trimnami"]["targets"]["envs"].append(
            os.path.join(config["trimnami"]["args"]["output_paths"]["temp"], filename + ".done")
        )
