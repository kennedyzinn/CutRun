library(GenomicRanges)
# load data

## define conditions of each sample
conditions <- c("R179+T", "R179-T")
replicates <- c("R1", "R2", "R3") 
sample_names <- unlist(lapply(conditions, function(cond){
  paste0(cond, "-EZH2_", replicates)
}))
condition_files <- lapply(sample_names, paste0, ".seacr.peaks.stringent.bed")

# read peak files
peaks <- lapply(condition_files, function(file) {
  if (file.size(file) > 0) {
    read.table(file, header=FALSE)
  } else {
    warning(paste("Skipping empty file:", file))
    return(NULL)
  }
})
peaks <- peaks[!sapply(peaks, is.null)]

# create GRanges objects
granges <- lapply(peaks, function(peak){
  GRanges(seqnames = peak$V1,
          ranges = IRanges(start = peak$V2, end = peak$V3))
})

# combine
all_peaks <- do.call(c, granges)
master_peaks <- reduce(all_peaks)

master_peaks$peak_id <- paste0("peak_", seq_along(master_peaks))

# define alignment files
USB_PATH="/Volumes/EZH2/"
align_path = sapply(USB_PATH, paste0, sample_names)
align_files <- sapply(align_path, paste0, ".target.markdup.sorted.bam")

#get counts
library(chromVAR)
count_matrix <- getCounts(alignment_files = align_files, peaks = master_peaks, paired = TRUE)
count_r179 <- count_matrix
