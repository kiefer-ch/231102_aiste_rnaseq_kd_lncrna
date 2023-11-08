# Project description

## Purpose

siRNA knockdown of three metabolically relevant lncRNA genes in primary human 
brown adipocytes. 


## Overall design

3 independent differentiations with cells from different subjects. 
Individual control samples for each lncRNA KD.
H19 samples were additionaly treated with norepinephrin.

# Methods

## Growth protocol


## Treatment protocol


## Extraction protocol


## Library construction and sequencing

CBMR mRNAseq: Profiling of the coding mRNA transcriptome. From total RNA, poly-A –captured RNA libraries are prepared and differential gene expression analysis can be performed. Sequencing is performed on NovaSeq6000 52bp paired-end, aiming for 30M reads per sample. 

# Data processing

## Base calling

## Quantification

Using salmon.

## Statistics

## Output files

### Count matrices

There are multiple count matrices to be found in the data/deseq/ folder:

* rld: log2(x + 1) level, between sample normalised, variance stabilised. See regularized log (rlog) transformation
[http://dx.doi.org/10.1186/s13059-014-0550-8]. This is the main table for any plotting and correlation analyses etc.

* tpm: between sample normalised and within sample (between genes) normalised.

* cts: raw counts, count output from tximport as imported from salmon estimated
counts: **No normalisation at all!** This one should be used for statistical
analyses in DESeq or EdgeR, which perform their own normalisation steps.

