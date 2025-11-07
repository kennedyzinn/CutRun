# spikein_utils.R
# Utility functions for extracting and processing E. coli spike-in counts
# from CUT&RUN Bowtie2 alignment logs
# 
# Usage in analysis script:
#   source("spikein_utils.R")
#   spikein_data <- get_spikein_counts("path/to/logs/")
#   sizeFactors(dds) <- spikein_data$size_factors

# ============================================================================
# Core Functions
# ============================================================================

#' Parse a single Bowtie2 log file to extract E. coli alignment counts
#'
#' @param log_file Path to a Bowtie2 log file
#' @return Integer count of aligned reads
#' @export
parse_bowtie2_spikein <- function(log_file) {
  if (!file.exists(log_file)) {
    stop("Log file not found: ", log_file)
  }
  
  log_lines <- readLines(log_file)
  
  # Extract alignment counts
  aligned_once_line <- grep("aligned concordantly exactly 1 time", 
                            log_lines, value = TRUE)
  aligned_multi_line <- grep("aligned concordantly >1 times", 
                             log_lines, value = TRUE)
  
  if (length(aligned_once_line) == 0 || length(aligned_multi_line) == 0) {
    stop("Could not parse alignment statistics from: ", log_file)
  }
  
  aligned_once <- as.numeric(gsub("^\\s*(\\d+).*", "\\1", aligned_once_line))
  aligned_multi <- as.numeric(gsub("^\\s*(\\d+).*", "\\1", aligned_multi_line))
  
  total_aligned <- aligned_once + aligned_multi
  
  return(total_aligned)
}


#' Extract E. coli spike-in counts from all log files in a directory
#'
#' @param log_dir Directory containing Bowtie2 log files
#' @param pattern File pattern to match (default: "\\.log$")
#' @param exclude_igg Whether to exclude IgG samples (default: TRUE)
#' @param verbose Print progress messages (default: TRUE)
#' @param strip_suffix Remove common suffixes from sample names (default: TRUE)
#' @return List containing spike-in counts and metadata
#' @export
get_spikein_counts <- function(log_dir, 
                               pattern = "\\.log$", 
                               exclude_igg = TRUE,
                               strip_suffix = TRUE,
                               verbose = TRUE) {
  
  # Validate input
  if (!dir.exists(log_dir)) {
    stop("Directory not found: ", log_dir)
  }
  
  # Get all log files
  log_files <- list.files(log_dir, pattern = pattern, full.names = TRUE)
  
  if (length(log_files) == 0) {
    stop("No log files found in: ", log_dir)
  }
  
  if (verbose) {
    cat("Found", length(log_files), "log files in", log_dir, "\n")
  }
  
  # Parse each log file
  ecoli_counts <- sapply(log_files, function(f) {
    tryCatch(
      parse_bowtie2_spikein(f),
      error = function(e) {
        warning("Failed to parse ", basename(f), ": ", e$message)
        return(NA)
      }
    )
  })
  
  # Clean sample names
  sample_names <- basename(log_files)
  sample_names <- gsub("\\.log$", "", sample_names)
  
  if (strip_suffix) {
    # Remove common pipeline suffixes
    sample_names <- gsub("\\.spikein\\.bowtie2$|\\.spikein$|\\.bowtie2$", "", sample_names)
    sample_names <- gsub("^spikein_|^ecoli_", "", sample_names)
  }
  
  names(ecoli_counts) <- sample_names
  
  # Remove any failed parses
  if (any(is.na(ecoli_counts))) {
    warning("Removing ", sum(is.na(ecoli_counts)), " samples with parsing errors")
    ecoli_counts <- ecoli_counts[!is.na(ecoli_counts)]
  }
  
  if (verbose) {
    cat("\n=== E. coli Spike-in Counts ===\n")
    print(ecoli_counts)
  }
  
  # Identify IgG samples
  igg_idx <- grep("igg|IGG", names(ecoli_counts), ignore.case = TRUE)
  
  result <- list(
    all_samples = ecoli_counts,
    experimental_samples = ecoli_counts,
    igg_samples = NULL,
    excluded_igg = FALSE
  )
  
  if (length(igg_idx) > 0) {
    if (verbose) {
      cat("\n=== IgG Samples Detected ===\n")
      cat("IgG samples:", names(ecoli_counts)[igg_idx], "\n")
    }
    
    if (exclude_igg) {
      if (verbose) {
        cat("Excluding IgG samples from normalization\n")
      }
      result$experimental_samples <- ecoli_counts[-igg_idx]
      result$igg_samples <- ecoli_counts[igg_idx]
      result$excluded_igg <- TRUE
    }
  }
  
  return(result)
}


#' Calculate spike-in normalization size factors
#'
#' @param spikein_counts Named vector or list of spike-in counts
#' @param sample_order Optional vector of sample names to reorder results
#' @param verbose Print progress messages (default: TRUE)
#' @param fuzzy_match Try to match samples even with suffix differences (default: TRUE)
#' @return List containing size factors and QC metrics
#' @export
calculate_spikein_size_factors <- function(spikein_counts, 
                                           sample_order = NULL,
                                           verbose = TRUE,
                                           fuzzy_match = TRUE) {
  
  # Handle list input (from get_spikein_counts)
  if (is.list(spikein_counts) && !is.null(spikein_counts$experimental_samples)) {
    counts <- spikein_counts$experimental_samples
  } else {
    counts <- spikein_counts
  }
  
  # Validate
  if (length(counts) == 0) {
    stop("No spike-in counts provided")
  }
  
  if (any(counts <= 0)) {
    warning("Some samples have zero or negative spike-in counts")
  }
  
  # Reorder if requested
  if (!is.null(sample_order)) {
    missing <- setdiff(sample_order, names(counts))
    
    # Try fuzzy matching if exact match fails
    if (length(missing) > 0 && fuzzy_match) {
      if (verbose) {
        cat("\nAttempting fuzzy matching for samples...\n")
      }
      
      # Create a mapping
      matched_names <- names(counts)
      for (target in sample_order) {
        if (!(target %in% names(counts))) {
          # Try to find a match by removing common suffixes
          pattern <- gsub("_R[123]$", "", target)  # Remove replicate suffix
          candidates <- grep(pattern, names(counts), value = TRUE, fixed = TRUE)
          
          if (length(candidates) == 1) {
            if (verbose) {
              cat("  Matched '", target, "' to '", candidates[1], "'\n", sep = "")
            }
            # Rename in counts
            idx <- which(names(counts) == candidates[1])
            names(counts)[idx] <- target
            matched_names[idx] <- target
          }
        }
      }
      
      # Check again
      missing <- setdiff(sample_order, names(counts))
    }
    
    if (length(missing) > 0) {
      cat("\nAvailable spike-in samples:\n")
      cat("  ", paste(names(counts), collapse = "\n  "), "\n")
      cat("\nRequested samples (from count matrix):\n")
      cat("  ", paste(sample_order, collapse = "\n  "), "\n")
      stop("Missing spike-in data for samples: ", paste(missing, collapse = ", "),
           "\nCheck that sample names match between count matrix and log files.")
    }
    
    extra <- setdiff(names(counts), sample_order)
    if (length(extra) > 0 && verbose) {
      cat("Removing extra samples not in sample_order: ", 
          paste(extra, collapse = ", "), "\n")
    }
    
    counts <- counts[sample_order]
  }
  
  # Calculate size factors
  size_factors <- counts / mean(counts)
  
  # QC metrics
  cv <- sd(counts) / mean(counts)
  fold_range <- max(counts) / min(counts)
  
  if (verbose) {
    cat("\n=== Spike-in Size Factors ===\n")
    print(round(size_factors, 3))
    cat("\n=== QC Metrics ===\n")
    cat("Coefficient of Variation:", round(cv, 3), "\n")
    cat("Fold range:", round(fold_range, 2), "\n")
    cat("Mean spike-in count:", round(mean(counts), 0), "\n")
  }
  
  # Warnings
  if (fold_range > 3) {
    warning("Large variation in spike-in counts (>3-fold). ",
            "Check for technical issues.")
  }
  
  if (cv > 0.5) {
    warning("High coefficient of variation (>0.5). ",
            "Spike-in normalization may be unreliable.")
  }
  
  return(list(
    size_factors = size_factors,
    spikein_counts = counts,
    cv = cv,
    fold_range = fold_range,
    mean_count = mean(counts)
  ))
}


#' Main wrapper function: Get spike-in counts and calculate size factors
#'
#' @param log_dir Directory containing Bowtie2 log files
#' @param sample_order Optional vector to order samples (matching count matrix)
#' @param exclude_igg Exclude IgG samples (default: TRUE)
#' @param strip_suffix Remove pipeline suffixes from sample names (default: TRUE)
#' @param fuzzy_match Try to match samples with different suffixes (default: TRUE)
#' @param verbose Print messages (default: TRUE)
#' @return List with size_factors, spikein_counts, and QC metrics
#' @export
get_spikein_size_factors <- function(log_dir,
                                     sample_order = NULL,
                                     exclude_igg = TRUE,
                                     strip_suffix = TRUE,
                                     fuzzy_match = TRUE,
                                     verbose = TRUE) {
  
  # Step 1: Extract counts from logs
  counts_data <- get_spikein_counts(log_dir, 
                                    exclude_igg = exclude_igg,
                                    verbose = verbose,
                                    strip_suffix = strip_suffix)
  
  # Step 2: Calculate size factors
  sf_data <- calculate_spikein_size_factors(counts_data$experimental_samples,
                                            sample_order = sample_order,
                                            verbose = verbose,
                                            fuzzy_match = fuzzy_match)
  
  # Combine results
  result <- c(sf_data, list(
    igg_samples = counts_data$igg_samples,
    all_samples = counts_data$all_samples,
    excluded_igg = counts_data$excluded_igg
  ))
  
  return(result)
}


# ============================================================================
# Diagnostic Functions
# ============================================================================

#' Generate QC plots for spike-in normalization
#'
#' @param spikein_data Output from get_spikein_size_factors()
#' @param condition Vector of condition labels matching samples
#' @export
plot_spikein_qc <- function(spikein_data, condition = NULL) {
  
  counts <- spikein_data$spikein_counts
  sf <- spikein_data$size_factors
  
  if (is.null(condition)) {
    condition <- rep("Sample", length(counts))
  }
  
  if (length(condition) != length(counts)) {
    stop("Length of condition vector must match number of samples")
  }
  
  # Set up plot layout
  par(mfrow = c(2, 2), mar = c(8, 4, 3, 2))
  
  # 1. Spike-in counts
  barplot(counts, 
          las = 2, 
          col = as.factor(condition),
          main = "E. coli Spike-in Counts",
          ylab = "Aligned Reads",
          cex.names = 0.7)
  legend("topright", unique(as.character(condition)),
         fill = 1:length(unique(condition)), cex = 0.8)
  
  # 2. Size factors
  barplot(sf,
          las = 2,
          col = as.factor(condition),
          main = "Spike-in Size Factors",
          ylab = "Size Factor",
          cex.names = 0.7)
  abline(h = 1, lty = 2, col = "red")
  
  # 3. Boxplot by condition
  if (length(unique(condition)) > 1) {
    boxplot(sf ~ condition,
            main = "Size Factors by Condition",
            ylab = "Size Factor",
            col = 1:length(unique(condition)))
    abline(h = 1, lty = 2, col = "red")
  } else {
    plot(1:length(sf), sf, 
         pch = 19, 
         col = as.factor(condition),
         main = "Size Factor Distribution",
         xlab = "Sample", 
         ylab = "Size Factor")
    abline(h = 1, lty = 2, col = "red")
  }
  
  # 4. Summary text
  plot.new()
  text(0.5, 0.8, "Spike-in QC Summary", cex = 1.5, font = 2)
  text(0.1, 0.6, paste0("Mean count: ", round(spikein_data$mean_count, 0)), 
       pos = 4, cex = 1.1)
  text(0.1, 0.5, paste0("CV: ", round(spikein_data$cv, 3)), 
       pos = 4, cex = 1.1)
  text(0.1, 0.4, paste0("Fold range: ", round(spikein_data$fold_range, 2)), 
       pos = 4, cex = 1.1)
  text(0.1, 0.3, paste0("N samples: ", length(counts)), 
       pos = 4, cex = 1.1)
  
  # Quality assessment
  status <- "GOOD"
  if (spikein_data$cv > 0.5 || spikein_data$fold_range > 3) {
    status <- "CAUTION"
  }
  text(0.1, 0.15, paste0("Status: ", status), 
       pos = 4, cex = 1.2, font = 2,
       col = if(status == "GOOD") "darkgreen" else "orange")
  
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))
}


#' Print detailed spike-in summary
#'
#' @param spikein_data Output from get_spikein_size_factors()
#' @export
print_spikein_summary <- function(spikein_data) {
  cat("\n╔════════════════════════════════════════════════════════╗\n")
  cat("║         SPIKE-IN NORMALIZATION SUMMARY                ║\n")
  cat("╚════════════════════════════════════════════════════════╝\n\n")
  
  cat("Samples included:", length(spikein_data$spikein_counts), "\n")
  
  if (!is.null(spikein_data$igg_samples)) {
    cat("IgG samples excluded:", length(spikein_data$igg_samples), "\n")
    cat("  IgG sample names:", paste(names(spikein_data$igg_samples), collapse = ", "), "\n")
  }
  
  cat("\nSpike-in Statistics:\n")
  cat("  Mean count:", round(spikein_data$mean_count, 0), "\n")
  cat("  Range:", range(spikein_data$spikein_counts), "\n")
  cat("  Coefficient of Variation:", round(spikein_data$cv, 3), "\n")
  cat("  Fold range:", round(spikein_data$fold_range, 2), "\n")
  
  cat("\nSize Factor Statistics:\n")
  cat("  Mean:", round(mean(spikein_data$size_factors), 3), "\n")
  cat("  Range:", round(range(spikein_data$size_factors), 3), "\n")
  
  cat("\nQuality Assessment:\n")
  if (spikein_data$cv > 0.5) {
    cat("  ⚠ WARNING: High CV (>0.5) - check for technical issues\n")
  } else {
    cat("  ✓ CV acceptable (<0.5)\n")
  }
  
  if (spikein_data$fold_range > 3) {
    cat("  ⚠ WARNING: Large fold range (>3) - verify sample quality\n")
  } else {
    cat("  ✓ Fold range acceptable (<3)\n")
  }
  
  cat("\n")
}


# ============================================================================
# Utility Functions
# ============================================================================

#' Preview a single log file
#'
#' @param log_file Path to log file
#' @export
preview_spikein_log <- function(log_file) {
  cat("═══════════════════════════════════════════════════════\n")
  cat("Log file:", basename(log_file), "\n")
  cat("═══════════════════════════════════════════════════════\n\n")
  
  log_lines <- readLines(log_file)
  cat(paste(log_lines, collapse = "\n"))
  
  cat("\n\n═══════════════════════════════════════════════════════\n")
  cat("Parsed E. coli count:", parse_bowtie2_spikein(log_file), "\n")
  cat("═══════════════════════════════════════════════════════\n")
}