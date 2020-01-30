# THIS SCRIPT IS CALLED FROM IDENTIFYSUPERENHANCERS.PY

# 1. sum signal from individual enhancers to stitched enhancers
# 2. rank enhancers based on signal
# 3. call super enhancers
# 4. output table and graphs

library(ggplot2)
library(plyr)
library(dplyr)

######################################################################

calculate_superenhancer_cutoff <- function(signal.vector){
  signal.vector <- sort(signal.vector)
  slope <- (max(signal.vector) - min(signal.vector))/length(signal.vector) # slope/diagonal we want to slide to find tangent line
  xPt <- floor(optimize(numPtsBelowLine, lower = 1, upper = length(signal.vector), input.vector = signal.vector, slope = slope)$minimum) # find the x-axis point where a line passing through it has the minimum number of points below it (i.e. tangent)
  y_cutoff <- signal.vector[xPt] #y-value at the x-point cutoff
}

numPtsBelowLine <- function(input.vector, slope, x){
  yPt <- input.vector[x]
  b <- yPt - (slope*x)
  xPts <- 1:length(input.vector)
  return(sum(input.vector <= (xPts*slope + b)))
}


# 1,2. sum signal from individual enhancers to stitched enhancers & rank super enhancers based on signal

# read in arguments
args <- commandArgs(trailingOnly = TRUE)
signal.file <- args[1] # TRE signal and parentage
prefix <- args[2]

signal.file <- paste0(substr(signal.file, 1, 57), "_parentage.txt")
signal <- read.delim(signal.file)

cat(as.character(Sys.time()), "Summing signal for each stitched enhancer...\n")

stitched <- signal %>%
  group_by(Parent) %>%
  summarize(TotalSignal = sum(Signal),
            NumTRE = n(),
            Children = paste(TRE, collapse = ";")) %>%
  mutate(chrom = sapply(strsplit(as.character(Parent), split = ":", fixed = TRUE), function(x) x[1]),
         start = sapply(strsplit(sapply(strsplit(as.character(Parent), split = ":", fixed = TRUE), function(x) x[2]), split = "-", fixed = TRUE), function(x) x[1]),
         stop = sapply(strsplit(as.character(Parent), split = "-", fixed = TRUE), function(x) x[2]),
         Length = as.numeric(stop) - as.numeric(start)) %>%
  arrange(desc(TotalSignal)) %>%
  mutate(Rank = 1:length(TotalSignal),
         RevRank = length(TotalSignal):1)

# 3. call super enhancers

cat(as.character(Sys.time()), "Calling super enhancers...\n")
SEcutoff <- calculate_superenhancer_cutoff(stitched$TotalSignal)

stitched <- stitched %>%
  mutate(Super = if_else(TotalSignal > SEcutoff, TRUE, FALSE))

# 4. write out files

cat(as.character(Sys.time()), "Writing out files...\n")
graphname <- paste0(prefix, "_SuperEnhancersRanked.png")
png(graphname, height = 5, width = 6, units = "in", res = 300)
ggplot(stitched, aes(x = RevRank, y = TotalSignal)) +
  geom_rect(data = NULL, aes(xmin = min(stitched$RevRank[which(stitched$TotalSignal > SEcutoff)]), xmax = Inf, ymin = -Inf, ymax = Inf), fill = "gray90") +
  geom_line(color = "red", size = 0.25) +
  geom_point(aes(color = Super)) +
  scale_color_manual(values = c("black", "red")) +
  theme_bw() +
  labs(y = "Transcriptional signal",
       x = "Ranked stitched enhancers") +
  theme(panel.grid = element_line(color = "white"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.75),
        axis.text = element_text(size = rel(1.5), color = "black"),
        axis.title = element_text(size = rel(1.5), color = "black"),
        axis.title.y = element_text(margin = margin(0,10,0,0)),
        axis.title.x = element_text(margin = margin(10,0,0,0))) +
  guides(color = FALSE)
dev.off()

stitched <- stitched %>%
  select(Parent, Rank, TotalSignal, NumTRE, Length, Super, Children)


outfilename = paste0(prefix, "_stitched_enhancers_info.txt")

write.table(stitched, outfilename, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
