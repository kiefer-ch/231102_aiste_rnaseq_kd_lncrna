# base url
GENCODE_URL = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_{}/release_{}/".format(
    config["SPECIES"], config["GENCODE_VERSION"])


# combined urls
GENOME_URL = GENCODE_URL + config["GENOME"] + ".gz"
TRANSCRIPTS_URL = GENCODE_URL + config["TRANSCRIPTS"] + ".gz"
ANNOTATION_URL = GENCODE_URL + config["ANNOTATION"] + ".gz"


rule get_transcripts:
    output:
        "data/annotation/{}".format(config["TRANSCRIPTS"])
    params:
        url = TRANSCRIPTS_URL
    resources:
        runtime = 5
    cache:
        True
    shell:
        "wget -q -O \
            - {params.url} \
            | gunzip > {output}"


rule get_genome:
    output:
        "data/annotation/{}".format(config["GENOME"])
    params:
        url = GENOME_URL
    resources:
        runtime = 5
    cache:
        True
    shell:
        "wget -q -O \
            - {params.url} \
            | gunzip > {output}"


rule get_annotation:
    output:
        "data/annotation/{}".format(config["ANNOTATION"])
    params:
        url = ANNOTATION_URL
    cache:
        True
    resources:
        runtime = 5
    shell:
        "wget -q -O \
            - {params.url} \
            | gunzip > {output}"
