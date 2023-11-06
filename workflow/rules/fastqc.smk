rule fastqc:
    input:
        fastq = "data/fastq/{state}/{sample_id}_{direction}_{state}.fastq.gz"
    output:
        multiext("data/qc/fastqc/{state}/{sample_id}_{direction}_{state}_fastqc", ".html", ".zip")
    params:
        outputDir = "data/qc/fastqc/{state}"
    conda:
        "../envs/fastqc.yaml"
    resources:
        runtime = 15,
        mem_mb = 2048
    shell:
        "fastqc {input.fastq} \
            --noextract \
            --quiet \
            --memory {resources.mem_mb} \
            -o {params.outputDir}"


rule fastqc_all:
    input:
        expand("data/qc/fastqc/{state}/{sample_id}_{direction}_{state}_fastqc.html",
            sample_id=SAMPLE_IDS,
            state=["raw", "trimmed"],
            direction=[1, 2])
    resources:
        runtime = 5
