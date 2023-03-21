# Libraries ====================================================================
library('tidyverse')
library('magrittr')
library('parallel')
library('optparse')

# Inputs =======================================================================
option_list <- list(
    make_option(
        c("--in_dir"),
        type = "character",
        help = "Target directory. Required."
    ),
    make_option(
        c("--out_dir"),
        type = "character",
        default = './',
        help = "Output directory. Default: [%default]"
    ),
    make_option(
        c("--id"),
        type = "character",
        default = "all_calls",
        help = "Output file name. Default: [%default]"
    ),
    make_option(
        c("--type"),
        type = "character",
        default = "segs",
        help = "Whether to check CNA or SEG output. Default: [%default]"
    ),
    make_option(
        c("--genes"),
        type = "character",
        default = NULL,
        help = "Genes to target. Optional."
    )
)
parse_options <- OptionParser(option_list = option_list)
opt <- parse_args(parse_options)

# Input processing
if (opt$genes %>% is.null) {
    gene_list <- c(
        'DNMT3A', 'ASXL1', 'TET2', 'PPM1D', 'TP53',
        'ATM', 'SRSF2', 'SF3B1', 'JAK2', 'GNB1'
    )
} else {
    gene_list <- opt$genes
}
gene_list %<>% toupper

cna_list <- list.files(opt$in_dir, pattern = 'cna')
cna_list <- paste0(opt$in_dir, '/', cna_list)

seg_list <- list.files(opt$in_dir, pattern = 'seg.txt')
seg_list <- paste0(opt$in_dir, '/', seg_list)

# Load data ====================================================================
# Load target gene info 
gene_ref <- '/storage1/fs1/timley/Active/aml_ppg/tmp/jonathanztang/breakpoint_reader/terra/for_git/data/gene_reference.bed'
gene_targ <- read.table(gene_ref) %>% as_tibble %>% filter(V5 %in% gene_list)
names(gene_targ) <- c('chr', 'start', 'end', 'code', 'gene', 'strand')
gene_targ$chr %<>% str_replace('chr', '')
gene_targ$chr %<>% as.numeric

# Overlaps =====================================================================
# Raw CNA files ----------------------------------------------------------------
call_cnv_cna <- function(cna_path, genes = gene_targ) {
    # Load cna file
    cna <- read.table(cna_path, header = TRUE) %>% as_tibble
    cna_names <- cna %>% names %>% str_replace('NWD[0-9].*\\.', '')
    names(cna) <- cna_names

    # Remove normal events
    cna %<>% filter(
            (chr %in% genes$chr) &
            (event != 'NEUT')
        )

    # Create empty call list
    call_list <- vector('list', nrow(genes))
    names(call_list) <- genes$gene

    # Iteratively call overlaps
    # TODO: Could probably use IRanges to make this more efficient
    for (gene_i in 1:nrow(genes)) {
        gene_row <- genes[gene_i, ]
        call_list[[gene_row$gene]] <- cna %>% filter(
            (chr == gene_row$chr) &
            ((end > gene_row$start) &
            (start < gene_row$end))
        )
        call_list[[gene_row$gene]]$gene <- gene_row$gene
    }
    
    # Collapse and return calls
    calls <- call_list %>% bind_rows
    calls$sample <- cna_path %>% str_extract('NWD[0-9]*')
    return(calls)
}

# Run batch with parallelization
batch_call_cna <- function(cna_list, genes = gene_targ) {
    all_calls <- mclapply(cna_list, call_cnv_cna) %>% bind_rows
    return(all_calls)
}

# Final seg files --------------------------------------------------------------
call_cnv_seg <- function(seg_path, genes = gene_targ) {
    # Load seg file
    seg <- read.table(seg_path, header = TRUE) %>% as_tibble

    # Remove normal events
    seg <- seg[
        (seg$chrom %in% genes$chr) &
        (seg$call != 'NEUT'),
    ]

    # Create empty call list
    call_list <- vector('list', nrow(genes))
    names(call_list) <- genes$gene

    # Iteratively call overlaps
    for (gene_i in 1:nrow(genes)) {
        gene_row <- genes[gene_i, ]
        call_list[[gene_row$gene]] <- seg[
            (seg$chrom == gene_row$chr) &
            (seg$end > gene_row$start) &
            (seg$start < gene_row$end),
        ]
        call_list[[gene_row$gene]]$gene <- gene_row$gene
    }

    # Collapse and return calls
    calls <- call_list %>% bind_rows
    return(calls)
}

# Run batch with parallelization
batch_call_seg <- function(seg_list, genes = gene_targ) {
    all_calls <- mclapply(seg_list, call_cnv_seg) %>% bind_rows
    return(all_calls)
}

# Script =======================================================================
if (opt$type == "segs") {
    output <- batch_call_seg(seg_list)
} else if (opt$type == "cna") {
    output <- batch_call_cna(cna_list)
}

write.table(
    output,
    file = paste0(opt$out_dir, "/", opt$id, ".txt"),
    quote = FALSE,
    row.names = FALSE,
    append = TRUE
)
