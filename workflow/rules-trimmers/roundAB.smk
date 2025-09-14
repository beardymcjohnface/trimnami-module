
rule trimnami_rm_5_primer:
    """Round A/B step 01: Remove 5' primer."""
    input:
        r1=os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.R1.fastq.gz"),
        r2=os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.R2.fastq.gz"),
        s=os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.RS.fastq.gz"),
        primers=workflow.source_path(
            os.path.join("..","..","resources",config["trimnami"]["qc"]["roundAB"]["primerB"])
        )
    output:
        r1=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.R1.s1.fastq")),
        r2=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.R2.s1.fastq")),
        s=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.S.s1.fastq")),
    benchmark:
        os.path.join(config["trimnami"]["args"]["output_paths"]["bench"],"trimnami_rm_5_primer.{file}.txt")
    log:
        os.path.join(config["trimnami"]["args"]["output_paths"]["log"],"trimnami_rm_5_primer.{file}.log")
    resources:
        **config["resources"]["med"]
    threads:
        config["resources"]["med"]["cpu"]
    params:
        config["trimnami"]["qc"]["bbduk"]["rm_5p"],
    conda:
        os.path.join("..", "envs", "bbmap.yaml")
    group:
        "roundAB"
    shell:
        """
        bbduk.sh \
            in={input.r1} \
            in2={input.r2} \
            ref={input.primers} \
            out={output.r1} \
            out2={output.r2} \
            threads={threads} \
            {params} \
            -Xmx{resources.mem_mb}m \
            2> {log}
        if [[ -s {input.s} ]]
        then
            bbduk.sh \
                in={input.s} \
                ref={input.primers} \
                out={output.s} \
                threads={threads} \
                {params} \
                -Xmx{resources.mem_mb}m \
                2> {log}
        else
            touch {output.s}
        fi
        """


rule trimnami_rm_3_contaminant:
    """Round A/B step 02: Remove 3' read through contaminant."""
    input:
        r1=os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.R1.s1.fastq"),
        r2=os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.R2.s1.fastq"),
        s=os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.S.s1.fastq"),
        primers=workflow.source_path(
            os.path.join("..","..","resources",config["trimnami"]["qc"]["roundAB"]["rc_primerB"])
        )
    output:
        r1=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.R1.s2.fastq")),
        r2=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.R2.s2.fastq")),
        s=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.S.s2.fastq")),
    benchmark:
        os.path.join(config["trimnami"]["args"]["output_paths"]["bench"],"trimnami_rm_3_contaminant.{file}.txt")
    log:
        os.path.join(config["trimnami"]["args"]["output_paths"]["log"],"trimnami_rm_3_contaminant.{file}.log")
    resources:
        **config["resources"]["med"]
    threads:
        config["resources"]["med"]["cpu"]
    params:
        config["trimnami"]["qc"]["bbduk"]["rm_3rt"]
    conda:
        os.path.join("..", "envs", "bbmap.yaml")
    group:
        "roundAB"
    shell:
        """
        bbduk.sh \
            in={input.r1} \
            in2={input.r2} \
            ref={input.primers} \
            out={output.r1} \
            out2={output.r2} \
            {params} \
            threads={threads} \
            -Xmx{resources.mem_mb}m \
            2> {log}
        if [[ -s {input.s} ]]
        then
            bbduk.sh \
                in={input.s} \
                ref={input.primers} \
                out={output.s} \
                {params} \
                threads={threads} \
                -Xmx{resources.mem_mb}m \
                2> {log}
        else
            touch {output.s}
        fi
        """


rule trimnami_rm_primer_free_adapter:
    """Round A/B step 03: Remove primer free adapter (both orientations)."""
    input:
        r1=os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.R1.s2.fastq"),
        r2=os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.R2.s2.fastq"),
        s=os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.S.s2.fastq"),
        primers=workflow.source_path(
            os.path.join("..","..","resources",config["trimnami"]["qc"]["roundAB"]["neb_adapters"])
        )
    output:
        r1=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.R1.s3.fastq")),
        r2=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.R2.s3.fastq")),
        s=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.S.s3.fastq")),
    benchmark:
        os.path.join(config["trimnami"]["args"]["output_paths"]["bench"],"trimnami_rm_primer_free_adapter.{file}.txt")
    log:
        os.path.join(config["trimnami"]["args"]["output_paths"]["log"],"trimnami_rm_primer_free_adapter.{file}.log")
    resources:
        **config["resources"]["med"]
    threads:
        config["resources"]["med"]["cpu"]
    params:
        config["trimnami"]["qc"]["bbduk"]["neb"]
    conda:
        os.path.join("..", "envs", "bbmap.yaml")
    group:
        "roundAB"
    shell:
        """
        bbduk.sh \
            in={input.r1} \
            in2={input.r2} \
            ref={input.primers} \
            out={output.r1} \
            out2={output.r2} \
            {params} \
            threads={threads} \
            -Xmx{resources.mem_mb}m \
            2> {log}
        if [[ -s {input.s} ]]
        then
            bbduk.sh \
                in={input.s} \
                ref={input.primers} \
                out={output.s} \
                {params} \
                threads={threads} \
                -Xmx{resources.mem_mb}m \
                2> {log}
        else
            touch {output.s}
        fi
        """


rule trimnami_rm_adapter_free_primer:
    """Round A/B step 04: Remove adapter free primer (both orientations)."""
    input:
        r1=os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.R1.s3.fastq"),
        r2=os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.R2.s3.fastq"),
        s=os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.S.s3.fastq"),
        primers=workflow.source_path(
            os.path.join("..","..","resources",config["trimnami"]["qc"]["roundAB"]["rc_primerB_ad6"])
        )
    output:
        r1=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.R1.s4.fastq")),
        r2=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.R2.s4.fastq")),
        s=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.S.s4.fastq")),
    benchmark:
        os.path.join(config["trimnami"]["args"]["output_paths"]["bench"],"trimnami_rm_adapter_free_primer.{file}.txt")
    log:
        os.path.join(config["trimnami"]["args"]["output_paths"]["log"],"trimnami_rm_adapter_free_primer.{file}.log")
    resources:
        **config["resources"]["med"]
    threads:
        config["resources"]["med"]["cpu"]
    params:
        config["trimnami"]["qc"]["bbduk"]["rm_afp"]
    conda:
        os.path.join("..", "envs", "bbmap.yaml")
    group:
        "roundAB"
    shell:
        """
        bbduk.sh \
            in={input.r1} \
            in2={input.r2} \
            ref={input.primers} \
            out={output.r1} \
            out2={output.r2} \
            {params} \
            threads={threads} \
            -Xmx{resources.mem_mb}m \
            2> {log}
        if [[ -s {input.s} ]]
        then
            bbduk.sh \
                in={input.s} \
                ref={input.primers} \
                out={output.s} \
                {params} \
                threads={threads} \
                -Xmx{resources.mem_mb}m \
                2> {log}
        else
            touch {output.s}
        fi
        """


rule trimnami_rm_vector_contamination:
    """Round A/B step 05: Vector contamination removal (PhiX + NCBI UniVecDB)"""
    input:
        r1=os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.R1.s4.fastq"),
        r2=os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.R2.s4.fastq"),
        s=os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.S.s4.fastq"),
        primers=workflow.source_path(
            os.path.join("..","..","resources",config["trimnami"]["qc"]["roundAB"]["vec_contam"])
        )
    output:
        r1=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.R1.s5.fastq")),
        r2=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.R2.s5.fastq")),
        s=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.S.s5.fastq")),
    benchmark:
        os.path.join(config["trimnami"]["args"]["output_paths"]["temp"]["bench"],"trimnami_rm_vector_contamination.{file}.txt")
    log:
        os.path.join(config["trimnami"]["args"]["output_paths"]["temp"]["log"],"trimnami_rm_vector_contamination.{file}.log")
    resources:
        **config["resources"]["med"]
    threads:
        config["resources"]["med"]["cpu"]
    params:
        config["trimnami"]["qc"]["bbduk"]["rm_vc"]
    conda:
        os.path.join("..", "envs", "bbmap.yaml")
    group:
        "roundAB"
    shell:
        """
        bbduk.sh \
            in={input.r1} \
            in2={input.r2} \
            ref={input.primers} \
            out={output.r1} \
            out2={output.r2} \
            {params} \
            threads={threads} \
            -Xmx{resources.mem_mb}m \
            2> {log}
        if [[ -s {input.s} ]]
        then
            bbduk.sh \
                in={input.s} \
                ref={input.primers} \
                out={output.s} \
                {params} \
                threads={threads} \
                -Xmx{resources.mem_mb}m \
                2> {log}
        else
            touch {output.s}
        fi
        """


rule trimnami_rm_low_qual:
    """Round A/B step 06: Remove remaining low-quality bases and short reads."""
    input:
        r1=os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.R1.s5.fastq"),
        r2=os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.R2.s5.fastq"),
        s=os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.S.s5.fastq"),
    output:
        r1=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.R1.s6.fastq")),
        r2=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.R2.s6.fastq")),
        s=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.S.s6.fastq")),
    benchmark:
        os.path.join(config["trimnami"]["args"]["output_paths"]["bench"],"trimnami_rm_low_qual.{file}.txt")
    log:
        os.path.join(config["trimnami"]["args"]["output_paths"]["log"],"trimnami_rm_low_qual.{file}.log")
    resources:
        **config["resources"]["med"]
    threads:
        config["resources"]["med"]["cpu"]
    params:
        config["trimnami"]["qc"]["bbduk"]["rm_lq"]
    conda:
        os.path.join("..", "envs", "bbmap.yaml")
    group:
        "roundAB"
    shell:
        """
        bbduk.sh \
            in={input.r1} \
            in2={input.r2} \
            out={output.r1} \
            out2={output.r2} \
            threads={threads} \
            {params} \
            -Xmx{resources.mem_mb}m \
            2> {log}
        if [[ -s {input.s} ]]
        then
            bbduk.sh \
                in={input.s} \
                out={output.s} \
                threads={threads} \
                {params} \
                -Xmx{resources.mem_mb}m \
                2> {log}
        else
            touch {output.s}
        fi
        """


rule trimnami_zip_roundAB:
    """Zip the final trimmed reads for Round A/B"""
    input:
        r1=os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.R1.s6.fastq"),
        r2=os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.R2.s6.fastq"),
        s=os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.S.s6.fastq"),
    output:
        r1=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.roundAB.R1.fastq.gz")),
        r2=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.roundAB.R2.fastq.gz")),
        s=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.roundAB.RS.fastq.gz")),
    benchmark:
        os.path.join(config["trimnami"]["args"]["output_paths"]["bench"],"trimnami_zip_roundAB.{file}.txt")
    log:
        os.path.join(config["trimnami"]["args"]["output_paths"]["log"],"trimnami_zip_roundAB.{file}.log")
    resources:
        **config["resources"]["med"]
    threads:
        config["resources"]["med"]["cpu"]
    params:
        compression = config["trimnami"]["qc"]["compression"]
    conda:
        os.path.join("..", "envs", "pigz.yaml")
    group:
        "roundAB"
    shell:
        """
        pigz -p {threads} -{params.compression} -c {input.r1} > {output.r1} 2> {log}
        pigz -p {threads} -{params.compression} -c {input.r2} > {output.r2} 2> {log}
        pigz -p {threads} -{params.compression} -c {input.s} > {output.s} 2> {log}
        """


rule trimnami_roundAB_single_end:
    """Round A/B for single end: This should not occur but this rule is here for testing purposes."""
    input:
        r1=os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.S.fastq.gz"),
    output:
        r1=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.roundAB.S.fastq.gz")),
        tmp=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.roundAB.S.fastq")),
    benchmark:
        os.path.join(config["trimnami"]["args"]["output_paths"]["bench"],"trimnami_roundAB_single_end.{file}.txt")
    log:
        os.path.join(config["trimnami"]["args"]["output_paths"]["log"],"trimnami_roundAB_single_end.{file}.log")
    resources:
        **config["resources"]["med"]
    threads:
        config["resources"]["med"]["cpu"]
    params:
        params = config["trimnami"]["qc"]["bbduk"]["rm_lq"],
        compression = config["trimnami"]["qc"]["compression"]
    conda:
        os.path.join("..", "envs", "bbmap.yaml")
    group:
        "roundAB"
    shell:
        """
        bbduk.sh \
            in={input.r1} \
            out={output.tmp} \
            threads={threads} \
            {params.params} \
            -Xmx{resources.mem_mb}m \
            2> {log}
        gzip -c -{params.compression} {output.tmp} > {output.r1}
        """
