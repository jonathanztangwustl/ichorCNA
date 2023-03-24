# ==============================================================================
library('tidyverse')
library('magrittr')
library('DNAcopy')
library('parallel')

# Load gene reference
GENE_REF <- '/storage1/fs1/timley/Active/aml_ppg/tmp/jonathanztang/breakpoint_reader/terra/for_git/data/gene_reference.bed'
ref <- read.table(
    GENE_REF,
    sep = '\t',
    header = FALSE,
    col.names = c('chromosome', 'start', 'end', 'name', 'gene', 'strand')
) %>% as_tibble

# Maybe do some further local segmentation
set_list <- list.files(pattern = '*.correctedDepth.txt')

# Perform DNAcopy CBS on a single sample
seg_cd <- function(set_path) {
    set_id <- set_path %>% str_extract('NWD[0-9]*')
    set_data <- read.table(set_path, header = T) %>% as_tibble
    names(set_data) <- c('chr', 'start', 'end', 'log2')
    set_cna <- CNA(
        set_data$log2,
        chrom = set_data$chr,
        maploc = set_data$start,
        data.type = 'logratio',
        sampleid = set_id
    )
    segs <- segment(
        set_cna,
        verbose = 0,
        alpha = 0.01,
        min.width = 3,
        undo.splits = 'sdund',
        undo.SD = 1.5
    )
    return(segs)
}

# Multiprocess DNAcopy CBS
batch_seg_cd <- function(set_list) {
    cluster <- makeForkCluster(detectCores() - 1)
    seg_list <- clusterMap(
        cluster,
        seg_cd,
        set_list
    )
    stopCluster(cluster)
    return(seg_list)
}

# Pull calls with minimum bin size and absolute segment mean
call_cd <- function(segs, bins = 13, seg_mean = 0.1) {
    seg_tb <- segs %>% lapply(`[[`, 2) %>% bind_rows %>% as_tibble
    calls <- seg_tb[
        (seg_tb$num.mark > bins) &
        (abs(seg_tb$seg.mean) > seg_mean),
    ]
    return(calls)
}

# Pull calls that overlap genes input 
call_gene <- function(calls, gene = 'DNMT3A') {
    gene %<>% toupper
    gene_info <- ref[(ref$gene %in% gene), ]
    gene_info$chromosome %<>% str_replace('chr', '')

    call_list <- vector('list', nrow(gene_info))
    for (gene_i in 1:nrow(gene_info)) {
        call_list[[gene_i]] <- calls[
            (calls$chrom == gene_info$chromosome[gene_i]) &
            (calls$loc.start <= gene_info$end[gene_i]) &
            (calls$loc.end >= gene_info$start[gene_i]),
        ]
        call_list[[gene_i]]$gene <- gene_info$gene[gene_i]
    }
    calls <- call_list %>% bind_rows
    return(calls)
}
