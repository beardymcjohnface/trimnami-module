import re
import os


def gzip_file(input_file, output_file):
    """
    Check if a file is gzipped and either zip it or copy the file

    Args:
        input_file (str): filepath of input file (gzipped or not)
        output_file (str): filepath of output gzipped file
    """
    fasta_regex = re.compile(r"^[a-zA-Z0-9_-]+\.(fa|fasta|fna|ffn|faa|frn)(\.gz)?$")
    if fasta_regex.match(input_file):
        os.system(f"seqtk -F I {input_file} | gzip -1 -c > {output_file}")
    else:
        os.system(f"seqtk {input_file} | gzip -1 -c > {output_file}")


def main(**kwargs):
    gzip_file(kwargs["input_r1"], kwargs["output_r1"])
    if "input_r2" in kwargs.keys():
        gzip_file(kwargs["input_r2"], kwargs["output_r2"])
    if "input_s" in kwargs.keys():
        open(kwargs["output_s"], "w").close()
        if kwargs["input_s"]:
            if os.path.exists(kwargs["input_s"]) and os.path.getsize(kwargs["input_s"]) > 0:
                gzip_file(kwargs["input_s"], kwargs["output_s"])


if __name__ == "__main__":
    if snakemake.params.is_paired:
        main(
            input_r1=snakemake.input.r1,
            input_r2=snakemake.input.r2,
            input_s=snakemake.params.s,
            output_r1=snakemake.output.r1,
            output_r2=snakemake.output.r2,
            output_s=snakemake.output.s
        )
    else:
        main(
            input_r1=snakemake.input.r1,
            output_r1=snakemake.output.r1,
        )
