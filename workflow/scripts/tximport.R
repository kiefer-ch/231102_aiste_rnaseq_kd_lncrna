source(".Rprofile")

suppressPackageStartupMessages({
    library("readr")
    library("dplyr")
    library("logger")
})

################################################################################
#
# author: Christoph Engelhard
# email: christoph.andreas.engelhard@regionh.dk
#
################################################################################

log_info("Generating transcript to gene table...")
read_gtf <- function(file) {
    readr::read_delim(file = file,
        delim = "\t",
        comment = "#",
        na = c('.'),
        col_names = c(
            "sequence", "source", "feature", "start", "end", "score",
            "strand", "phase", "attributes"),
        col_types = readr::cols(
            sequence = readr::col_character(),
            source = readr::col_character(),
            feature = readr::col_character(),
            start = readr::col_integer(),
            end = readr::col_integer(),
            score = readr::col_character(),
            strand = readr::col_character(),
            phase = readr::col_character(),
            attributes = readr::col_character()),
        progress = FALSE) %>%
    dplyr::filter(strand %in% c('+', '-')) %>%
    dplyr::mutate(
        feature = as.factor(feature),
        strand = as.factor(strand))
}

separate_transcripts <- function(gtf) {
    gtf %>%
        filter(feature == "transcript") %>%
        tidyr::separate(
            attributes,
            c("gene_id", "transcript_id", "gene_type", "gene_name",
              "transcript_type", "transcript_name", "level", "extra"),
            sep = "; ",
            extra = "merge") %>%
        mutate(
            gene_id = substring(gene_id, 10L, nchar(gene_id) - 1L),
            gene_name = substring(gene_name, 12L, nchar(gene_name) - 1L),
            transcript_id = substring(transcript_id, 16L, nchar(transcript_id) - 1L))
}

tx2g <- read_gtf(snakemake@input$annotation) %>%
    separate_transcripts() %>%
    dplyr::select(transcript_id, gene_id)


log_info("Preparing sample_info...")
salmon <- tibble(path = snakemake@input$salmon) %>%
    mutate(sample_id = basename(dirname(path)))

sample_info <- suppressMessages(read_csv(snakemake@input$runinfo)) %>%
    right_join(salmon, by = "sample_id") %>%
    mutate(across(c("sample_nr", "condition", "treatment", "siRNA", "subject_id"), as.factor))


log_info("Create and output dds...")
outfile <- snakemake@output[[1]]
dir.create(dirname(outfile), showWarnings = FALSE)

tximport::tximport(files = sample_info$path,
        type = "salmon",
        tx2gene = tx2g,
        txOut = FALSE) %>%
    DESeq2::DESeqDataSetFromTximport(.,
        colData = tibble::column_to_rownames(as.data.frame(sample_info), "sample_id"),
        design = as.formula(snakemake@params$design)) %>%
    write_rds(outfile, "gz", compression = 5L)

log_success("Done.")
