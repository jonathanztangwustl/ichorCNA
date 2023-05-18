library('tidyverse')
library('magrittr')
library('DNAcopy')
library('parallel')

# FUNCTIONS ====================================================================
# Perform DNAcopy CBS on a single sample
seg_cd <- function(set_path) {
    set_id <- set_path %>% str_extract('NWD[0-9]*')
    set_data <- read.table(set_path, header = T) %>% as_tibble
    names(set_data) <- c('chr', 'start', 'end', 'log2')
    set_data %<>% filter(chr == 2) # TEMPORARY FOR DNMT3A
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

# Prerequisites for plot
GAPS_BED <- '/storage1/fs1/timley/Active/aml_ppg/tmp/jonathanztang/breakpoint_reader/terra/for_git/data/ichor_gaps.bed'
GAPS <- read.table(GAPS_BED, header = TRUE) %>% as_tibble
names(GAPS) <- c('chr', 'start', 'end')
GAPS$start <- as.numeric(GAPS$start)
GAPS$end <- as.numeric(GAPS$end)

# Plot single sample
GENE_REF <- '/storage1/fs1/timley/Active/aml_ppg/tmp/jonathanztang/breakpoint_reader/terra/for_git/data/gene_reference.bed'
plot_raw <- function(seg_data, gene = 'DNMT3A', out_dir = './plots', cov_type = 'indexcov', mean_thresh = 0.05, bin_thresh = 10) {
    gene_ref <- read.table(GENE_REF, sep = '\t', header = FALSE, col.names = c('chromosome', 'start', 'end', 'name', 'zero', 'strand'))
    names(gene_ref) <- c('chromosome', 'start', 'end', 'label', 'name', 'strand')
    gene <- toupper(gene)
    gene_info <- gene_ref[gene_ref$name == gene, ]
    if (nrow(gene_info) > 1) {stop('Only one gene can be specified.')}

    chrom_lim <- gene_info$chromosome %>% str_replace('chr', '') %>% as.numeric

    points <- seg_data$data %>% as_tibble %>% filter(chrom %in% chrom_lim)
    names(points) <- c('chrom', 'maploc', 'log2')
    segments <- seg_data$output %>% as_tibble %>% filter(chrom %in% chrom_lim)
    set_name <- segments$ID %>% unique

    out_plot <- ggplot(data = points, aes(x = maploc, y = log2)) +
        geom_point(size = 0) +
        geom_segment(data = segments,
            aes(
                x = loc.start,
                xend = loc.end,
                y = seg.mean,
                yend = seg.mean
                ),
            size = 1,
            color = 'red'
            ) +
        ylim(-2, 2) +
        geom_segment(data = (GAPS %>% filter(chr %in% gene_info$chromosome)), aes(x = start, xend = end, y = -0.5, yend = -0.5), color = 'green') +
        geom_text(data = segments, aes(x = loc.start, y = seg.mean, label = ifelse((abs(seg.mean) >= mean_thresh & num.mark >= bin_thresh), '*', '')), color = 'dodgerblue', nudge_y = -0.075) +
        geom_vline(xintercept = gene_info$start, color = 'blue', size = 0.1) +
        geom_vline(xintercept = gene_info$end, color = 'blue', size = 0.1) +
        geom_hline(yintercept = 0, color = 'black', size = 0.1) +
        ggtitle(paste(set_name, gene, cov_type)) + xlab('Base Pair') + ylab('log2')
    ggsave(paste0(out_dir, '/', set_name, '_', gene_info$chromosome, '_', gene_info$name, '_segments.png'), out_plot, device = 'png')
}

# Plot batch of segment data
# plot_raw_batch <- function(deletion_data, inputs) {
#     plot_out_dir <- paste0(inputs$out_dir, '/plots/')
#     if (!file.exists(plot_out_dir)) {
#         dir.create(plot_out_dir, recursive = TRUE)
#     }
#     cluster <- makeForkCluster(detectCores() - 1)
#     clusterMap(cluster,
#         plot_raw,
#         deletion_data$data,
#         MoreArgs = list(
#             gene = inputs$gene,
#             out_dir = plot_out_dir,
#             cov_type = inputs$type,
#             mean_thresh = inputs$mean_threshold,
#             bin_thresh = inputs$bin_threshold
#         ))
#     stopCluster(cluster)
# }


# SCRIPT =======================================================================
samples <- list.files(pattern = 'correctedDepth')

segs <- batch_seg_cd(samples)

for (seg_i in 1:length(segs)) {
    plot_raw(segs[[seg_i]])
}
