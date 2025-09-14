rule trimnami_save_output:
    """Save the final trimmed output file"""
    input:
        os.path.join(config["trimnami"]["args"]["output_paths"]["temp"], "{file}"),
    output:
        os.path.join(config["trimnami"]["args"]["output_paths"]["results"], "{file}")
    localrule: True
    shell:
        "mv {input} {output}"
        "ln -s $(pwd)/{output} $(pwd)/{input}"


rule trimnami_init_input_paired_end:
    """Initialise the input files for paired end reads"""
    input:
        r1 = lambda wildcards: samples["reads"][wildcards.sample]["R1"],
        r2 = lambda wildcards: samples["reads"][wildcards.sample]["R2"],
    output:
        r1 = temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{sample}.R1.fastq.gz")),
        r2 = temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{sample}.R2.fastq.gz")),
        s = temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{sample}.RS.fastq.gz")),
    params:
        s = lambda wildcards: samples["reads"][wildcards.sample]["S"],
    conda:
        os.path.join("..", "envs", "seqtk.yaml")
    shell:
        """
        process_reads_file() {{
            local input_file="$1"
            local output_file="$2"
        
            if gzip -t "$input_file" 2>/dev/null; then
                local catcom="zcat"
            else
                local catcom="cat"
            fi
        
            if "$catcom" "$input_file" | head -1 | grep -Pq "^@"; then
                local seqtkparm=""
            else
                local seqtkparm="-F I"
            fi

            "$catcom" "$input_file" | seqtk "$seqtkparm" - | gzip -1 - > "$output_file"
        }}
        
        process_reads_file {input.r1} {output.r1}
        process_reads_file {input.r2} {output.r2}
        
        if [ -s {params.s} ]; then
            process_reads_file {params.s} {output.s}
        else
            touch {output.s}
        fi
        """


rule trimnami_init_input_single_end:
    """Initialise the input files for paired end reads"""
    input:
        r1=lambda wildcards: samples["reads"][wildcards.sample]["R1"],
    output:
        r1=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{sample}.S.fastq.gz")),
    conda:
        os.path.join("..","envs","seqtk.yaml")
    shell:
        """
        process_reads_file() {{
            local input_file="$1"
            local output_file="$2"

            if gzip -t "$input_file" 2>/dev/null; then
                local catcom="zcat"
            else
                local catcom="cat"
            fi

            if "$catcom" "$input_file" | head -1 | grep -Pq "^@"; then
                local seqtkparm=""
            else
                local seqtkparm="-F I"
            fi

            "$catcom" "$input_file" | seqtk "$seqtkparm" - | gzip -1 - > "$output_file"
        }}

        process_reads_file {input.r1} {output.r1}
        """