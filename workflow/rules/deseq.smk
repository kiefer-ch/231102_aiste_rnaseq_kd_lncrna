rule tximport:
    input:
        salmon = expand("data/salmon/{sample_id}/quant.sf",
            sample_id=SAMPLE_IDS),
        annotation = "data/annotation/{}".format(config["ANNOTATION"]),
        runinfo = "config/sample_info.csv"
    output:
        "data/deseq/dds.rds"
    log:
        "log/deseq/tximport.log"
    params:
        design = "~ subject_id + condition"
    conda:
        "../envs/R_4.3.yaml"
    resources:
        runtime = 15,
        mem_mb = 10000
    script:
        "../scripts/tximport.R"


rule export_cm:
    input:
        dds = "data/deseq/dds.rds",
        annotation = "data/annotation/{}".format(config["ANNOTATION"])
    output:
        cts = "results/cm_cts.csv.gz",
        rld = "results/cm_rld.csv.gz",
        tpm = "results/cm_tpm.csv.gz"
    conda:
        "../envs/R_4.3.yaml"
    resources:
        runtime = 15,
        mem_mb = 10000
    script:
        "../scripts/export_cm.R"


rule report:
    input:
        dds = "data/deseq/dds.rds",
        annotation = "data/annotation/{}".format(config["ANNOTATION"])
    output:
        "results/report.html"
    threads:
        6
    script:
        "../scripts/report.Rmd"
