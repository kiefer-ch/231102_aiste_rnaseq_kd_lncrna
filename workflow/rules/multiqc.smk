rule multiqc:
    input:
        expand(
            "data/qc/fastqc/{state}/{sample_id}_{direction}_{state}_fastqc.zip",
            sample_id=SAMPLE_IDS,
            direction=["1", "2"],
            state=["raw", "trimmed"]),
        expand("log/cutadapt/{sample_id}.log",
             sample_id = SAMPLE_IDS),
        expand("data/salmon/{sample_id}/",
            sample_id = SAMPLE_IDS)
    output:
        "results/multiqc_data.zip",
        html = "results/multiqc.html"
    log:
        "log/multiqc/multiqc.log"
    conda:
        "../envs/multiqc.yaml",
    resources:
        runtime = 5
    script:
        "../scripts/multiqc.py"
