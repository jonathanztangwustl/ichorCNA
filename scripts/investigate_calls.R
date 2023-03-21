library('tidyverse')
library('magrittr')

# Look at DNMT3A calls specifically
calls <- read.table('ichor_run_dnmt3a.txt', header = T) %>% as_tibble
hetd <- calls %>% filter(call == 'HETD')

hetd_plot <- ggplot(hetd, aes(x = logR_Copy_Number)) + geom_histogram()
#ggsave('hetd_plot.png', hetd_plot, 'png')

all_plot <- ggplot(calls, aes(x = logR_Copy_Number)) +
    geom_histogram(bins = 100) +
    geom_density(aes(y=100*..density..))
#ggsave('all_plot.png', all_plot, 'png')

# Generating test dataset ------------------------------------------------------
# Let's look at the distribution of number of calls per sample
calls <- read.table('ichor_run_seg.txt', header = T) %>% as_tibble
calls$ID %<>% str_extract('NWD[0-9]*')
counts <- calls %>% count(ID)
count_plot <- ggplot(counts, aes(x = n)) +
    geom_histogram(bins = 30) +
    geom_density(aes(y=2000*..density..)) +
    xlim(0,100)
ggsave('count_plot.png', count_plot, 'png')

# Maybe limit our samples of interest to have < n calls total, arbitrarily?
pl_ids <- counts %>% filter(n < 25) %>% .$ID
pl_calls <- calls %>% filter(ID %in% pl_ids)

gene_list <- 'DNMT3A'
# Load target gene info 
gene_ref <- '/storage1/fs1/timley/Active/aml_ppg/tmp/jonathanztang/breakpoint_reader/terra/for_git/data/gene_reference.bed'
gene_targ <- read.table(gene_ref) %>% as_tibble %>% filter(V5 %in% gene_list)
names(gene_targ) <- c('chr', 'start', 'end', 'code', 'gene', 'strand')
gene_targ$chr %<>% str_replace('chr', '')
gene_targ$chr %<>% as.numeric

pl_dnmt3a <- pl_calls[
    (pl_calls$chrom == gene_targ$chr) &
    (pl_calls$start < gene_targ$end) &
    (pl_calls$end > gene_targ$start),
]

pl_dels <- pl_dnmt3a %>% filter(call == 'HETD')
write.table(pl_dels$ID, 'dnmt3a_25.txt', quote = F, row.names = F)

# Further analysis -------------------------------------------------------------
calls <- read.table('ichor_run_seg.txt', header = T) %>% as_tibble
calls$ID %<>% str_extract('NWD[0-9]*')
counts <- calls %>% count(ID)

# Maybe remove upper outliers with a simple 1.5*IQR?
quan <- quantile(counts$n)
iqr <- quan[4] - quan[2]
upper <- quan[4] + (1.5 * iqr)

pl_ids <- counts %>% filter(n <= upper) %>% .$ID
pl_calls <- calls %>% filter(ID %in% pl_ids)

gene_list <- 'DNMT3A'
# Load target gene info 
gene_ref <- '/storage1/fs1/timley/Active/aml_ppg/tmp/jonathanztang/breakpoint_reader/terra/for_git/data/gene_reference.bed'
gene_targ <- read.table(gene_ref) %>% as_tibble %>% filter(V5 %in% gene_list)
names(gene_targ) <- c('chr', 'start', 'end', 'code', 'gene', 'strand')
gene_targ$chr %<>% str_replace('chr', '')
gene_targ$chr %<>% as.numeric

pl_dnmt3a <- pl_calls[
    (pl_calls$chrom == gene_targ$chr) &
    (pl_calls$start < gene_targ$end) &
    (pl_calls$end > gene_targ$start),
]

# One file with everything
write.table(pl_dnmt3a, 'dnmt3a_54.txt', quote = F, row.names = F)

# One with just deletions
pl_dnmt3a_dels <- pl_dnmt3a %>% filter(Corrected_Call == 'HETD')
write.table(pl_dnmt3a_dels, 'dnmt3a_dels_54.txt', quote = F, row.names = F)
