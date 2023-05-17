# INITIALIZATION ===============================================================
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

goi <- c('DNMT3A', 'TET2', 'ASXL1', 'PPM1D', 'TP53', 'ATM', 'SRSF2', 'SF3B1',
    'JAK2', 'GNB1')

# FUNCTIONS ====================================================================
# Chunk median
chunk_med <- function(in_set_data, chunk) {
    return(in_set_data$log2[chunk] %>% median(na.rm = T))
}

# Test waviness
wavy <- function(set_path) {
    set_id <- set_path %>% str_extract('NWD[0-9]*')
    set_data <- read.table(set_path, header = T) %>% as_tibble
    names(set_data) <- c('chr', 'start', 'end', 'log2')
    set_bins <- split(
        1:nrow(set_data),
        ceiling(
            seq_along(1:nrow(set_data))/100
        )
    )
    waviness <- lapply(set_bins, chunk_med, in_set_data = set_data) %>%
        unlist %>% diff # diff helps with large CNVs that have non-0 baselines
    median_dev <- waviness %>% sd
    return(median_dev)
}

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
    segs$output$waviness <- wavy(set_path)
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
# Previous thresholds were bin 13, seg_mean 0.1
# Theoretical thresholds should be bins 10 (163 kb) and seg 0.05 (3.5% VAF)
call_cd <- function(segs, bins = 10, seg_mean = 0.05) {
    seg_tb <- segs %>% lapply(`[[`, 2) %>% bind_rows %>% as_tibble
    calls <- seg_tb[
        (seg_tb$num.mark > bins) &
        (abs(seg_tb$seg.mean) > seg_mean),
    ]
    return(calls)
}

# Pull calls that overlap genes input 
call_gene <- function(calls, gene = goi) {
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
    calls <- call_list %>% bind_rows %>% arrange(ID, chrom, loc.start)
    return(calls)
}

# INPUTS =======================================================================
library('optparse')

inputs_list <- list(
    make_option(
        c('--in_dir', '-i'),
        type = 'character',
        help = 'Input directory to analyze. Required.'
    ),
    make_option(
        c('--out_dir', '-o'),
        type = 'character',
        default = './segmentation/',
        help = 'Output directory. Optional. Default: [%default]'
    ),
    make_option(
        c('--gene', '-g'),
        type = 'character',
        default = 'DNMT3A',
        help = 'Gene or genes to search for, separated by commas (e.g. DNMT3A,TET2,JAK2). Optional. Default: [%default]'
    ),
    make_option(
        c('--bins', '-b'),
        type = 'integer',
        default = 10,
        help = 'Minimum bin threshold for calling segments. Optional. Default: [%default]'
    ),
    make_option(
        c('--mean', '-m'),
        type = 'double',
        default = 0.05,
        help = 'Minimum segment mean change to call segments. Optional. Default: [%default]'
    ),
    make_option(
        c('--save_calls', '-s'),
        action = 'store_true',
        help = 'Add this option to save all calls as well as gene calls.'
    )
)
parser <- OptionParser(option_list = inputs_list)
inputs <- parse_args(parser)

# TESTING =====================================================================
#inputs <- list(
#    in_dir = 'bigcbs',
#    out_dir = './test',
#    gene = 'dnmt3a,tet2',
#    bins = 13,
#    mean = 0.1,
#    save_calls = TRUE
#)

# SCRIPT =======================================================================
# Parse and tweak inputs
print('Parsing inputs.')
set_list <- list.files(inputs$in_dir, pattern = '*.correctedDepth.txt') %>%
    paste0(inputs$in_dir, '/', .)
inputs$gene %<>% toupper %>% str_split(',') %>% .[[1]]
if (!dir.exists(inputs$out_dir)) {
    dir.create(inputs$out_dir, recursive = TRUE)
}

# Run workflow on sets of 200 to reduce memory
sets <- split(set_list, ceiling(seq_along(set_list)/200))
print(paste0('Running workflow, ', length(sets), ' set(s).'))

for (set_i in 1:length(sets)) {
    print('Segmenting...')
    calls <- batch_seg_cd(sets[[set_i]]) %>%
        call_cd(bins = inputs$bins, seg_mean = inputs$mean)
    
    if (inputs$save_calls) {
        print(paste0(
            'Writing all calls, ', set_i, ' of ', length(sets), '...'
        ))
        all_call_path <- paste0(inputs$out_dir, '/all_calls.txt')
        write.table(
            calls,
            all_call_path,
            append = T,
            quote = F,
            row.names = F,
            col.names = !file.exists(all_call_path)
        )
    }
    gene_calls <- calls %>% call_gene(inputs$gene)
    print(paste0(
        'Writing gene calls, ', set_i, ' of ', length(sets), '...'
    ))
    gene_call_path <- paste0(inputs$out_dir, '/gene_calls.txt')
    write.table(
        gene_calls,
        gene_call_path,
        append = T,
        quote = F,
        row.names = F,
        col.names = !file.exists(gene_call_path)
    )
}
print('Complete.')
