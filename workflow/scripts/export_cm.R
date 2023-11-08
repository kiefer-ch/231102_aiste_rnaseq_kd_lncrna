source(".Rprofile")

suppressPackageStartupMessages({
    library("readr")
    library("logger")
    library("dplyr")
    library("DESeq2")
})

################################################################################
#
# author: Christoph Engelhard
# email: christoph.andreas.engelhard@regionh.dk
#
################################################################################

read_gtf <- function(file) {
    readr::read_delim(
        file = file,
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


separate_genes <- function(gtf) {
    gtf %>%
        filter(feature == "gene") %>%
        tidyr::separate(
            attributes,
            c("gene_id", "gene_type", "gene_name", "level", "extra"),
            sep = "; ",
            extra = "merge") %>%
        mutate(
            gene_id = substring(gene_id, 10L, nchar(gene_id) - 1L),
            gene_name = substring(gene_name, 12L, nchar(gene_name) - 1L))
}


log_info("Import data...")
dds <- readRDS(snakemake@input$dds)


# rlog normalisation
log_info("Normalising data...")
rld <- rlog(dds, blind = FALSE)

# write to disc
assay(rld) %>%
    as_tibble(rownames = "gene_id") %>%
    left_join(gene2symbol, by = "gene_id") %>%
    dplyr::select(gene_id, symbol, everything()) %>%
    write_csv(snakemake@output$rld)


# read gtf to get symbols
log_info("Getting gene2symbol table...")
gene2symbol <- read_gtf(snakemake@input$annotation) %>%
    separate_genes() %>%
    dplyr::select(gene_id, gene_name) %>%
    dplyr::rename(symbol = gene_name)


# raw counts
log_info("Exporting raw counts...")
counts(dds, normalized = FALSE) %>%
    as_tibble(rownames = "gene_id") %>%
    left_join(gene2symbol, by = "gene_id") %>%
    dplyr::select(gene_id, symbol, everything()) %>%
    write_csv(snakemake@output$cts)


# tpm
# https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/
fpkmToTpm <- function(fpkm) {
    exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

log_info("Exporting TPM...")
tpm <- fpkm(dds, robust = TRUE)
tpm <- apply(tpm, 2, fpkmToTpm)

tpm %>%
    as_tibble(rownames = "gene_id") %>%
    left_join(gene2symbol, by = "gene_id") %>%
    dplyr::select(gene_id, symbol, everything()) %>%
    write_csv(snakemake@output$tpm)

log_success("Done.")
