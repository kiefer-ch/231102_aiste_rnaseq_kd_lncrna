from os import path
from snakemake.shell import shell
import tempfile


output_dir = path.dirname(snakemake.output.html)
output_name = path.basename(snakemake.output.html)
log = snakemake.log_fmt_shell(stdout=True, stderr=True)


lines = [
    "module_order:",
    "    -salmon",
    "    - fastqc:",
    "        name: 'FastQC (trimmed)'",
    "        info: 'This section of the report shows FastQC results after adapter removal and base quality trimming.'",
    "        target: ''",
    "        path_filters:",
    "            - '*trimmed_fastqc.zip'",
    "    - cutadapt",
    "    - fastqc:",
    "        name: 'FastQC (raw)'",
    "        info: 'This section of the report shows FastQC results before adapter removal and base quality trimming.'",
    "        target: ''",
    "        path_filters:",
    "            - '*raw_fastqc.zip'",
    "fastqc_config:",
    "    fastqc_theoretical_gc: 'hg38_txome'",
    "extra_fn_clean_exts:",
    "    - type: remove",
    "      pattern: _raw"
]


tmp = tempfile.NamedTemporaryFile()
with open(tmp.name, 'w') as f:
    f.write('\n'.join(lines))

shell(
    "(multiqc -f -z "
    "-c {tmp.name} "
    "-o {output_dir} "
    "-n {output_name} "
    "{snakemake.input} "
    ") {log}"
)
