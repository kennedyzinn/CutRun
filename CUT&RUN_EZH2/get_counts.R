library(GenomicRanges)
# load data

## define conditions of each sample
conditions <- c("R179+T", "R179-T", "WT+T", "WT-T")
replicates <- C("R1", "R2")
sample_names <- unlist(lapply(conditions, function(cond){
  paste0(cond, "-EZH2_", replicates)
}))
condition_files <- lapply(sample_names, paste0, ".seacr.peaks.stringent.bed")

# read peak files
peaks <- lapply(condition_files, read.table, header=FALSE, sep="\t")

# create GRanges objects
granges <- lapply(peaks, function(peak){
  GRanges(seqnames = peak$V1,
          ranges = IRanges(start = peak$V2, end = peak$V3))
})

# combine
all_peaks <- do.call(c, granges)
master_peaks <- reduce(all_peaks) # saved in folder as master_peaks.bed. Not sure when I saved this either, possible I did it in AWS

master_peaks$peak_id <- paste0("peak_", seq_along(master_peaks))

# define alignment files
USB_PATH="/Volumes/EZH2/"
align_path = sapply(USB_PATH, paste0, sample_names)
align_files <- sapply(align_path, paste0, ".target.markdup.sorted.bam")

#get counts
library(chromVAR)
count_matrix <- getCounts(alignment_files = align_files, peaks = master_peaks, paired = TRUE)
write.csv(count_matrix, "dense_countMat.csv") # not sure but pretty sure this R 
# script is how I got the dense_countMat.csv. I think I accidentally save the 
# count_matrix object using the console instead of including this line in the 
# script