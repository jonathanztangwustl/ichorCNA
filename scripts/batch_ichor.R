# Script to parallelize multiple ichor jobs, passing along input arguments

library('optparse')
library('tidyverse')
library('magrittr')
library('parallel')

# Option list direct from runIchorCNA.R with a directory function
option_list <- list(
    make_option(c("--inDir"), type = "character", help = "Input directory with all samples. Required."),
    # make_option(c("--WIG"), type = "character", help = "Path to tumor WIG file. Required."),
    make_option(c("--NORMWIG"), type = "character", default=NULL, help = "Path to normal WIG file. Default: [%default]"),
    make_option(c("--gcWig"), type = "character", default="/scratch1/fs1/timley/fusions/jonathanztang/references/ichor.fa.gc.wig", help = "Path to GC-content WIG file; Required"),
    make_option(c("--mapWig"), type = "character", default="/scratch1/fs1/timley/fusions/jonathanztang/references/ichor.fa.map.ws_16384.wig", help = "Path to mappability score WIG file. Default: [%default]"),
    make_option(c("--normalPanel"), type="character", default=NULL, help="Median corrected depth from panel of normals. Default: [%default]"),
    make_option(c("--exons.bed"), type = "character", default=NULL, help = "Path to bed file containing exon regions. Default: [%default]"),
    make_option(c("--id"), type = "character", default="test", help = "Patient ID. Default: [%default]"),
    make_option(c("--centromere"), type="character", default="/storage1/fs1/timley/Active/aml_ppg/tmp/jonathanztang/breakpoint_reader/terra/for_git/data/201_gaps.bed", help = "File containing Centromere locations; if not provided then will use hg19 version from ichorCNA package. Default: [%default]"),
    make_option(c("--minMapScore"), type = "numeric", default=0.9, help="Include bins with a minimum mappability score of this value. Default: [%default]."),
    make_option(c("--rmCentromereFlankLength"), type="numeric", default=1e5, help="Length of region flanking centromere to remove. Default: [%default]"),
    make_option(c("--normal"), type="character", default="0.5", help = "Initial normal contamination; can be more than one value if additional normal initializations are desired. Default: [%default]"),
    make_option(c("--scStates"), type="character", default="NULL", help = "Subclonal states to consider. Default: [%default]"),
    make_option(c("--coverage"), type="numeric", default=NULL, help = "PICARD sequencing coverage. Default: [%default]"),
    make_option(c("--lambda"), type="character", default="NULL", help="Initial Student's t precision; must contain 4 values (e.g. c(1500,1500,1500,1500)); if not provided then will automatically use based on variance of data. Default: [%default]"),
    make_option(c("--lambdaScaleHyperParam"), type="numeric", default=3, help="Hyperparameter (scale) for Gamma prior on Student's-t precision. Default: [%default]"),
    #	make_option(c("--kappa"), type="character", default=50, help="Initial state distribution"),
    make_option(c("--ploidy"), type="character", default="2", help = "Initial tumour ploidy; can be more than one value if additional ploidy initializations are desired. Default: [%default]"),
    make_option(c("--maxCN"), type="numeric", default=7, help = "Total clonal CN states. Default: [%default]"),
    make_option(c("--estimateNormal"), type="logical", default=TRUE, help = "Estimate normal. Default: [%default]"),
    make_option(c("--estimateScPrevalence"), type="logical", default=TRUE, help = "Estimate subclonal prevalence. Default: [%default]"),
    make_option(c("--estimatePloidy"), type="logical", default=TRUE, help = "Estimate tumour ploidy. Default: [%default]"),
    make_option(c("--maxFracCNASubclone"), type="numeric", default=0.7, help="Exclude solutions with fraction of subclonal events greater than this value. Default: [%default]"),
    make_option(c("--maxFracGenomeSubclone"), type="numeric", default=0.5, help="Exclude solutions with subclonal genome fraction greater than this value. Default: [%default]"),
    make_option(c("--minSegmentBins"), type="numeric", default=50, help="Minimum number of bins for largest segment threshold required to estimate tumor fraction; if below this threshold, then will be assigned zero tumor fraction."),
    make_option(c("--altFracThreshold"), type="numeric", default=0.05, help="Minimum proportion of bins altered required to estimate tumor fraction; if below this threshold, then will be assigned zero tumor fraction. Default: [%default]"),
    make_option(c("--chrNormalize"), type="character", default="c(1:22)", help = "Specify chromosomes to normalize GC/mappability biases. Default: [%default]"),
    make_option(c("--chrTrain"), type="character", default="c(1:22)", help = "Specify chromosomes to estimate params. Default: [%default]"),
    make_option(c("--chrs"), type="character", default="c(1:22,\"X\")", help = "Specify chromosomes to analyze. Default: [%default]"),
    make_option(c("--genomeBuild"), type="character", default="hg38", help="Geome build. Default: [%default]"),
    make_option(c("--genomeStyle"), type = "character", default = "NCBI", help = "NCBI or UCSC chromosome naming convention; use UCSC if desired output is to have \"chr\" string. [Default: %default]"),
    make_option(c("--normalizeMaleX"), type="logical", default=TRUE, help = "If male, then normalize chrX by median. Default: [%default]"),
    make_option(c("--minTumFracToCorrect"), type="numeric", default=0.1, help = "Tumor-fraction correction of bin and segment-level CNA if sample has minimum estimated tumor fraction. [Default: %default]"), 
    make_option(c("--fracReadsInChrYForMale"), type="numeric", default=0.001, help = "Threshold for fraction of reads in chrY to assign as male. Default: [%default]"),
    make_option(c("--includeHOMD"), type="logical", default=FALSE, help="If FALSE, then exclude HOMD state. Useful when using large bins (e.g. 1Mb). Default: [%default]"),
    make_option(c("--txnE"), type="numeric", default=0.9999999, help = "Self-transition probability. Increase to decrease number of segments. Default: [%default]"),
    make_option(c("--txnStrength"), type="numeric", default=1e7, help = "Transition pseudo-counts. Exponent should be the same as the number of decimal places of --txnE. Default: [%default]"),
    make_option(c("--plotFileType"), type="character", default="pdf", help = "File format for output plots. Default: [%default]"),
    make_option(c("--plotYLim"), type="character", default="c(-2,2)", help = "ylim to use for chromosome plots. Default: [%default]"),
    make_option(c("--outDir"), type="character", default="./", help = "Output Directory. Default: [%default]"),
    make_option(c("--libdir"), type = "character", default=NULL, help = "Script library path. Usually exclude this argument unless custom modifications have been made to the ichorCNA R package code and the user would like to source those R files. Default: [%default]")
)
parseobj <- OptionParser(option_list=option_list)
opt <- parse_args(parseobj)
options(scipen=0, stringsAsFactors=F)

# Test inputs
opt <- list(
    inDir = '.',
    outDir = '../analysis/ichortest/'
)

# Parse system calls
# TODO: Add other arguments, including segment inputs and panel of normals
if (!dir.exists(opt$outDir)) {
    dir.create(opt$outDir, recursive = TRUE)
}
sample_list <- opt$inDir %>% list.dirs %>% str_subset('NWD') %>% str_extract('NWD[0-9]*')
calls <- paste0(
    "Rscript $psc/scripts/ichorCNA/scripts/runIchorCNA.R --WIG ./",
    sample_list,
    "/indexcov.tar.gz --id ",
    sample_list,
    " --outDir ",
    opt$outDir
)

cluster <- makeForkCluster(detectCores() - 1)
clusterMap(cluster,
    system,
    calls
)
stopCluster(cluster)
