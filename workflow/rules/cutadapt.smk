rule cutadapt_pe:
    input:
        fastq_fw = "data/fastq/raw/{sample_id}_1_raw.fastq.gz",
        fastq_rv = "data/fastq/raw/{sample_id}_2_raw.fastq.gz"
    output:
        fastq_fw = "data/fastq/trimmed/{sample_id}_1_trimmed.fastq.gz",
        fastq_rv = "data/fastq/trimmed/{sample_id}_2_trimmed.fastq.gz",
        qc = "data/qc/cutadapt/{sample_id}.cutadapt.json"
    log:
        "log/cutadapt/{sample_id}.log"
    threads:
        5
    conda:
        "../envs/cutadapt.yaml"
    resources:
        runtime = 15,
        mem_mb = 500
    shell:
        "cutadapt \
            -j {threads} \
            -q 28 \
            -m 30 \
            -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
            -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
            -o {output.fastq_fw} \
            -p {output.fastq_rv} \
            {input.fastq_fw} {input.fastq_rv} \
            --json={output.qc} \
            > {log}"


rule trim_all:
    input:
        expand("data/fastq/trimmed/{sample_id}_{direction}_trimmed.fastq.gz",
            sample_id=SAMPLE_IDS,
            direction=[1, 2])
    resources:
        runtime = 5
