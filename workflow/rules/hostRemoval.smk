rule trimnami_index_host_genome:
    """Pre-index the host genome for mapping with minimap2"""
    input:
        config["trimnami"]["args"]["host"]
    output:
        config["trimnami"]["args"]["host"] + ".idx"
    params:
        config["trimnami"]["qc"]["minimapIndex"]
    resources:
        **config["resources"]["med"]
    threads:
        config["resources"]["med"]["cpu"]
    conda:
        os.path.join("..", "envs","minimap2.yaml")
    benchmark:
        os.path.join(config["trimnami"]["args"]["output_paths"]["bench"],"trimnami_index_host_genome.txt")
    log:
        os.path.join(config["trimnami"]["args"]["output_paths"]["log"],"trimnami_index_host_genome.log")
    shell:
        "minimap2 -t {threads} {params} -d {output} {input} &> {log}\n\n"


rule trimnami_host_rm_mapping_paired:
    """Map reads to host and return unmapped reads"""
    input:
        r1=os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.R1.fastq.gz"),
        r2=os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.R2.fastq.gz"),
        s=os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.RS.fastq.gz"),
        host=config["trimnami"]["args"]["host"] + ".idx"
    output:
        r1=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.host_rm.R1.fastq.gz")),
        r2=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.host_rm.R2.fastq.gz")),
        rs=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"{file}.host_rm.RS.fastq.gz")),
        s=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"rm.{file}_s.host_rm.fastq.gz")),
        o=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"],"rm.{file}_o.host_rm.fastq.gz")),
        O=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"], "rm.{file}_O.host_rm.fastq.gz")),
    params:
        compression=config["trimnami"]["qc"]["compression"],
        minimap_mode=config["trimnami"]["args"]["minimap"],
        flagFilt=config["trimnami"]["qc"]["hostRemoveFlagstat"]
    benchmark:
        os.path.join(config["trimnami"]["args"]["output_paths"]["bench"],"host_removal_mapping.{file}.txt")
    log:
        mm=os.path.join(config["trimnami"]["args"]["output_paths"]["log"],"host_removal_mapping.{file}.minimap.log"),
        sv=os.path.join(config["trimnami"]["args"]["output_paths"]["log"],"host_removal_mapping.{file}.samtoolsView.log"),
        fq=os.path.join(config["trimnami"]["args"]["output_paths"]["log"],"host_removal_mapping.{file}.samtoolsFastq.log"),
    resources:
        **config["resources"]["med"]
    threads:
        config["resources"]["med"]["cpu"]
    conda:
        os.path.join("..", "envs","minimap2.yaml")
    shell:
        "minimap2 "
            "-ax {params.minimap_mode} "
            "-t {threads} "
            "--secondary=no "
            "{input.host} {input.r1} {input.r2} "
            "2> {log.mm} "
        "| samtools view "
            "-h {params.flagFilt} "
            "2> {log.sv} "
        "| samtools fastq "
            "-N -O -c {params.compression} "
            "-1 {output.r1} "
            "-2 {output.r2} "
            "-0 {output.O} "
            "-s {output.rs} "
            "2> {log.fq}\n\n "
        "cat {output.O} >> {output.rs}\n\n "
        "if [[ -s {input.s} ]]\n "
        "then "
            "minimap2 "
                "-ax {params.minimap_mode} "
                "-t {threads} "
                "--secondary=no "
                "{input.host} "
                "{input.s} "
                "2> {log.mm}"
            "| samtools view "
                "-h {params.flagFilt} "
                "2> {log.sv}"
            "| samtools fastq "
                "-n -O -c {params.compression} "
                "-o {output.o} "
                "-0 {output.O} "
                "-s {output.s} "
                "2> {log.fq}\n\n "
            "cat {output.o} {output.O} {output.s} >> {output.rs}\n\n "
        "else "
            "touch {output.o} {output.O} {output.s}\n "
        "fi\n\n "


rule trimnami_host_rm_mapping_single:
    """Map reads to host and return unmapped reads"""
    input:
        r1=os.path.join(config["trimnami"]["args"]["output_paths"]["temp"], "{file}.S.fastq.gz"),
        host=config["trimnami"]["args"]["host"] + ".idx"
    output:
        rs=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"], "{file}.host_rm.S.fastq.gz")),
        r0=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"], "{file}.host_rm.r0.fastq.gz")),
        ro=temp(os.path.join(config["trimnami"]["args"]["output_paths"]["temp"], "{file}.host_rm.rO.fastq.gz")),
    params:
        compression=config["trimnami"]["qc"]["compression"],
        minimap_mode=config["trimnami"]["args"]["minimap"],
        flagFilt=config["trimnami"]["qc"]["hostRemoveFlagstat"]
    benchmark:
        os.path.join(config["trimnami"]["args"]["output_paths"]["bench"],"host_removal_mapping.{file}.txt")
    log:
        mm=os.path.join(config["trimnami"]["args"]["output_paths"]["log"],"host_removal_mapping.{file}.minimap.log"),
        sv=os.path.join(config["trimnami"]["args"]["output_paths"]["log"],"host_removal_mapping.{file}.samtoolsView.log"),
        fq=os.path.join(config["trimnami"]["args"]["output_paths"]["log"],"host_removal_mapping.{file}.samtoolsFastq.log")
    resources:
        **config["resources"]["med"]
    threads:
        config["resources"]["med"]["cpu"]
    conda:
        os.path.join("..", "envs","minimap2.yaml")
    shell:
        "minimap2 "
            "-ax {params.minimap_mode} "
            "-t {threads} "
            "--secondary=no "
            "{input.host} "
            "{input.r1} "
            "2> {log.mm} "
        "| samtools view "
            "-h {params.flagFilt} "
        "2> {log.sv} "
        "| samtools fastq "
            "-n -O -c {params.compression} "
            "-o {output.ro} "
            "-0 {output.r0} "
            "-s {output.rs} "
            "2> {log.fq}\n\n "
        "cat {output.r0} {output.ro} >> {output.rs}\n\n "
