def get_rename_input(wildcards):
    sample_line = SAMPLE_INFO[SAMPLE_INFO["sample_id"] == wildcards.sample_id]
    sample_nr = sample_line["sample_nr"].tolist()[0]
    fastq = "data/fastq/raw/0252_X1_{}_S{}_R{}_001.fastq.gz".format(sample_nr, sample_nr, wildcards.direction)
    return fastq


rule rename_fastq:
    input:
        get_rename_input
    output:
        "data/fastq/raw/{sample_id}_{direction}_raw.fastq.gz"
    wildcard_constraints:
        direction = "1|2"
    resources:
        runtime = 1,
        mem_mb = 500
    shell:
        "ln -sr {input} {output}"
