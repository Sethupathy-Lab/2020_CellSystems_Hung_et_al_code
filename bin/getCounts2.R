# THIS SCRIPT IS CALLED FROM TRE2GENE.PY

# 1. get signal from specified genomic regions (TREs or gene body)
# 2. normalize signal to account for differences in sequencing depth between samples
# 3. output table of counts for each sample at each distal TRE

library(bigWig)
library(parallel)
library(DESeq2)

######################################################################

# 1. get signal from specified genomic regions

# create function to get counts from plus and minus strand bigwig files (for TREs)
getCounts <- function(plus, minus, intervals, path = "") {
  pl <- load.bigWig(paste(path, plus, sep = ""))
  mn <- load.bigWig(paste(path, minus, sep = ""))
  counts.plus <- bed.region.bpQuery.bigWig(pl, intervals, abs.value = TRUE)
  counts.minus <- bed.region.bpQuery.bigWig(mn, intervals, abs.value = TRUE)
  counts <- (counts.plus + counts.minus)
}

# create function to get counts from plus OR minus strand bigwig files (for gene bodies)
getStrandedCountsGenes <- function(plus, minus, intervals, path = "") {
  plus.intervals <- intervals[intervals[,6] == "+",]
  minus.intervals <- intervals[intervals[,6] == "-",]
  pl <- load.bigWig(paste(path, plus, sep = ""))
  mn <- load.bigWig(paste(path, minus, sep = ""))
  counts.plus <- bed.region.bpQuery.bigWig(pl, plus.intervals, abs.value = TRUE)
  counts.minus <- bed.region.bpQuery.bigWig(mn, minus.intervals, abs.value = TRUE)
  counts <- c(counts.plus, counts.minus)
}

# read in arguments
args <- commandArgs(trailingOnly = TRUE)
bigwig.files <- args[1] # list of bigwig files
coord.file <- args[2] # gene body/TRE coordinates
thing <- args[3] # TRE or gene body

plus.bw <- read.delim(bigwig.files, header = FALSE)
plus.bw <- as.character(plus.bw$V1)
minus.bw <- gsub("_plus.bw", "_minus.bw", plus.bw)

coords <- read.delim(coord.file, header = FALSE, stringsAsFactors = FALSE)
coords <- coords[which(coords$V1 != "chrM"),]
if (ncol(coords) == 4) {
	coords$coord <- paste(coords$V4, paste(paste(coords$V1, coords$V2, sep = ":"), coords$V3, sep = "-"), sep = "|")
} else {
	coords$coord <- paste(paste(coords$V1, coords$V2, sep = ":"), coords$V3, sep = "-")
}

if (!all(file.exists(plus.bw)))
	stop(paste("One or more bigwig plus strand files do not exist(", plus.bw, ")"));
if (!all(file.exists(minus.bw)))
	stop(paste("One or more bigwig minus strand files do not exist(", minus.bw, ")"));

# extract counts
cat(as.character(Sys.time()), "Extracting read counts...\n")

if (thing == "TRE") {
  coord.counts <- cbind(mclapply(1:length(plus.bw), function(x) getCounts(plus.bw[x], minus.bw[x], intervals = coords[,1:3]), mc.cores = 2))
  coord.counts <- data.frame(sapply(coord.counts, c))
  row.names(coord.counts) <-  coords$coord
} else if (thing == "gene") {
  coord.counts <- cbind(mclapply(1:length(plus.bw), function(x) getStrandedCountsGenes(plus.bw[x], minus.bw[x], intervals = coords), mc.cores = 2))
  coord.counts <- data.frame(sapply(coord.counts, c))
  row.names(coord.counts) <- c(coords$V4[coords[,6] == "+"], coords$V4[coords[,6] == "-"])
} else {
  stop("3rd argument of call to getCounts2.R was neither 'TRE' nor 'gene'")
}

# 2. normalize for differences in sequencing depth between samples

# import raw counts into DESeq2
cat(as.character(Sys.time()), "Normalizing read counts with DESeq2...\n") 
colData <- data.frame(sample = c(rep("Sample", ncol(coord.counts))))
row.names(colData) <- names(coord.counts)
dds <- DESeqDataSetFromMatrix(coord.counts, colData = colData, design = ~ 1)
dds <- estimateSizeFactors(dds)
norm.coord.counts <- as.data.frame(counts(dds, normalized = TRUE))
# remove all loci without any reads
norm.coord.counts <- norm.coord.counts[rowSums(norm.coord.counts) > 0,]
# perform log2(norm counts + 1) transformation
log2.norm.coord.counts <- log2(norm.coord.counts + 1)

# 3. output table of mean normalized counts at each distal TRE

outfilename = paste(substr(coord.file, 1, nchar(coord.file) - 4), "Counts.txt", sep = "")

write.table(log2.norm.coord.counts, outfilename, quote = FALSE, sep = "\t", col.names = FALSE, row.names = TRUE)
