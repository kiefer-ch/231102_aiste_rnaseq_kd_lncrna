rule salmon_prepareDecoys:
    input:
        transcripts = "data/annotation/{}".format(config["TRANSCRIPTS"]),
        genome = "data/annotation/{}".format(config["GENOME"]),
        annotation = "data/annotation/{}".format(config["ANNOTATION"])
    output:
        temp("data/indices/salmon_decoy/gentrome.fa"),
        temp("data/indices/salmon_decoy/decoys.txt")
    log:
        "log/salmon/generateDecoyTranscriptome.log"
    params:
        out_dir = "data/indices/salmon_decoy"
    conda:
        "../envs/mashmap.yaml"
    threads:
        40
    resources:
        runtime = 240,
        mem_mb = 160000
    shell:
        "workflow/scripts/generateDecoyTranscriptome.sh \
            -j {threads} \
            -g {input.genome} \
            -t {input.transcripts} \
            -a {input.annotation} \
            -o {params.out_dir} \
            > {log}"


rule salmon_index:
    input:
        gentrome = "data/indices/salmon_decoy/gentrome.fa",
        decoy = "data/indices/salmon_decoy/decoys.txt"
    output:
        directory("data/indices/salmon_index_{}".format(config["TRANSCRIPTS"]))
    log:
        "log/salmon/salmon_index.log"
    conda:
        "../envs/salmon.yaml"
    threads:
        20
    resources:
        runtime = 10,
        mem_mb = 10000
    shell:
        "salmon index \
            --gencode \
            -t {input.gentrome} \
            -i {output} \
            -d {input.decoy} \
            -p {threads} \
            2> {log}"


rule salmon_quantify:
    input:
        index = directory("data/indices/salmon_index_{}".format(config["TRANSCRIPTS"])),
        fastq_fw = "data/fastq/trimmed/{sample_id}_1_trimmed.fastq.gz",
        fastq_rv = "data/fastq/trimmed/{sample_id}_2_trimmed.fastq.gz"
    output:
        "data/salmon/{sample_id}/quant.sf",
        "data/salmon/{sample_id}/aux_info/meta_info.json",
        "data/salmon/{sample_id}/aux_info/fld.gz"
    log:
        "log/salmon/{sample_id}.log"
    conda:
        "../envs/salmon.yaml"
    threads:
        15
    resources:
        runtime = 10,
        mem_mb = 10000
    shell:
        "salmon quant \
            --gcBias \
            --seqBias \
            --numGibbsSamples 25 \
            -i {input.index} \
            -l A \
            -1 {input.fastq_fw} \
            -2 {input.fastq_rv} \
            -p {threads} \
            --validateMappings \
            -o data/salmon/{wildcards.sample_id} \
            2> {log}"


rule salmon_quantify_all:
    input:
        expand("data/salmon/{sample_id}/quant.sf",
            sample_id = SAMPLE_IDS)
