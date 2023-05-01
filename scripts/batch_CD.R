# Script to parallelize multiple ichor jobs, passing along input arguments

library('optparse')
library('tidyverse')
library('magrittr')
library('parallel')

# Change other settings directly in runIchorCNA.R
option_list <- list(
    make_option(c("--fileList"), type = "character", help = "File list with paths to indexcov files. Required."),
    make_option(c("--outDir"), type = "character", default = "./", help = "Output Directory. Default: [%default]"),
    make_option(c("--usepon"), type = "logical", default = TRUE, help = "Whether to use PON. Default: [%default]")
    #make_option(c("--plotType"), type = "character", default = "png", help = "Output plot file type. Default: [%default]")
)
parseobj <- OptionParser(option_list=option_list)
opt <- parse_args(parseobj)
options(scipen=0, stringsAsFactors=F)

# Test inputs
#opt <- list(
#    fileList = 'test.txt',
#    outDir = '../analysis/ichortest/'
#)

# Parse system calls
# TODO: Add other arguments, including segment inputs and panel of normals
if (!dir.exists(opt$outDir)) {
    dir.create(opt$outDir, recursive = TRUE)
}
sample_list <- readLines(opt$fileList)
sample_names <- sample_list %>% str_extract('NWD[0-9]*')
calls <- paste0(
    "Rscript $psc/scripts/ichorCNA/scripts/run_CD.R --WIG ",
    sample_list,
    "/indexcov.tar.gz --id ",
    sample_names,
    " --outDir ",
    opt$outDir
)
if (opt$usepon) {
    calls <- paste0(calls, " --normalPanel /storage1/fs1/timley/Active/aml_ppg/tmp/jonathanztang/breakpoint_reader/terra/outputs/pharmu/ichor_mutect_pon_median.rds ")
}

cluster <- makeForkCluster(detectCores() - 1)
clusterMap(cluster,
    system,
    calls
)
stopCluster(cluster)

# TODO: Remove extra outputs
