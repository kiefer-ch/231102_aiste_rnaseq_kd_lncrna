__author__ = "Christoph Engelhard"
__email__ = "christoph.andreas.engelhard@regionh.dk"

import pandas as pd


# configuration
configfile: "config/config.yaml"

SAMPLE_INFO = pd.read_csv("config/sample_info.csv")
SAMPLE_IDS = SAMPLE_INFO["sample_id"].tolist()

# global wildcard constraints
wildcard_constraints:
    type = "genome|transcriptome",
    state = "raw|trimmed"


# include
include: "workflow/rules/rename_input.smk"
include: "workflow/rules/fastqc.smk"
include: "workflow/rules/cutadapt.smk"
include: "workflow/rules/gencode.smk"
include: "workflow/rules/salmon.smk"
include: "workflow/rules/multiqc.smk"
