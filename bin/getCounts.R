# THIS SCRIPT IS CALLED FROM IDENTIFYSUPERENHANCERS.PY

# 1. get signal from distal TREs (enhancers)
# 2. normalize signal to account for differences in sequencing depth between samples
# 3. output table of counts for each sample at each distal TRE

library(bigWig)
library(parallel)
library(DESeq2)

######################################################################

# 1. get signal from distal TREs (enhancers)

# create function to get counts from plus and minus strand bigwig files
getCounts <- function(plus, minus, intervals, path = "") {
  pl <- load.bigWig(paste(path, plus, sep = ""))
  mn <- load.bigWig(paste(path, minus, sep = ""))
  counts.plus <- bed.region.bpQuery.bigWig(pl, intervals, abs.value = TRUE)
  counts.minus <- bed.region.bpQuery.bigWig(mn, intervals, abs.value = TRUE)
  counts <- (counts.plus + counts.minus)
}

# read in arguments
args <- commandArgs(trailingOnly = TRUE)
bigwig.files <- args[1] # list of bigwig files
enhancers.file <- args[2] # distal TRE coordinates

plus.bw <- read.delim(bigwig.files, header = FALSE)
plus.bw <- as.character(plus.bw$V1)
minus.bw <- gsub("_plus.bw", "_minus.bw", plus.bw)

enhancers <- read.delim(enhancers.file, header = FALSE)
enhancers$coord <- paste(paste(enhancers$V1, enhancers$V2, sep = ":"), enhancers$V3, sep = "-")

if (!all(file.exists(plus.bw)))
	stop(paste("One or more bigwig plus strand files do not exist(", plus.bw, ")"));
if (!all(file.exists(minus.bw)))
	stop(paste("One or more bigwig minus strand files do not exist(", minus.bw, ")"));

# extract counts
cat(as.character(Sys.time()), "Extracting read counts at distal TREs...\n")
enhancer.counts <- cbind(mclapply(1:length(plus.bw), function(x) getCounts(plus.bw[x], minus.bw[x], intervals = enhancers[,1:3]), mc.cores = 2))
enhancer.counts <- data.frame(sapply(enhancer.counts, c))
row.names(enhancer.counts) <- enhancers$coord

# 2. normalize for differences in sequencing depth between samples

# import raw counts into DESeq2
cat(as.character(Sys.time()), "Normalizing read counts at distal TREs with DESeq2...\n") 
colData <- data.frame(sample = c(rep("Sample", ncol(enhancer.counts))))
row.names(colData) <- names(enhancer.counts)
dds <- DESeqDataSetFromMatrix(enhancer.counts, colData = colData, design = ~ 1)
dds <- estimateSizeFactors(dds)
norm.enhancer.counts <- as.data.frame(counts(dds, normalized = TRUE))
norm.enhancer.counts$Mean <- rowMeans(norm.enhancer.counts)
norm.enhancer.counts$Var <- rowVars(as.matrix(norm.enhancer.counts))
norm.enhancer.counts$TRE <- enhancers$coord
outfile <- norm.enhancer.counts[,c('TRE', 'Mean')]

# 3. output table of mean normalized counts at each distal TRE

outfilename = paste(substr(enhancers.file, 1, 57), "_counts.txt", sep = "")

write.table(outfile, outfilename, quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
