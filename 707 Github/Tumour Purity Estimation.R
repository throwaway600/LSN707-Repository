# TPES-style Tumor Purity Estimation (Locallo et al., 2019)


# Load required packages
suppressPackageStartupMessages({
  library(readr)        
  library(dplyr)        
  library(GenomicRanges)
  library(dbscan)       
  library(FNN)          
  library(purrr)        
})

# Load input data

# Somatic SNV file
snv_df <- read_tsv("TARGET-OS.somaticmutation_wxs.tsv.gz")

# ASCAT CNV file 
cnv_df <- read_tsv("TARGET-OS.allele_cnv_ascat2.tsv.gz")

# Pick Eps Kneedle Function - DBSCAN

pick_eps_kneedle <- function(x, k = 2) {
  # If not enough points to estimate neighborhood distance, return NA
  if (length(x) <= k) return(NA_real_)
  
  # Compute k-th nearest neighbor distance for each point
  d <- FNN::knn.dist(matrix(x, ncol = 1), k = k)[, k]
  
  # Sort distances ascendingly
  ds <- sort(d)
  n <- length(ds)
  
  # If very few points, fallback to median distance
  if (n < 5) return(stats::median(ds))
  
  # Normalize x (rank positions) and y (distances) to [0,1]
  xs <- seq_len(n)
  xsn <- (xs - min(xs)) / (max(xs) - min(xs))
  dsn <- (ds - min(ds)) / (max(ds) - min(ds) + .Machine$double.eps)
  
  # Compute difference between normalized distance curve and diagonal
  # The point of maximum difference = "elbow" in k-distance plot
  gap <- dsn - xsn
  idx <- which.max(gap)
  eps <- ds[idx]
  
  # Ensure eps is non-negative
  eps <- max(0, eps)
  eps
}


# Helper function 2: Merge Similar Clusters based on Wilcoxon Rank-Sum Test
# Merges DBSCAN clusters whose VAF distributions are not statistically different

merge_similar_clusters <- function(vafs, clusters, alpha = 0.05) {
  vdf <- data.frame(vaf = vafs, cl = clusters)
  
  repeat {
    clv <- sort(unique(vdf$cl))
    if (length(clv) <= 1) break  # stop if only one cluster remains
    
    # Create all possible cluster pairs
    pairs <- combn(clv, 2, simplify = FALSE)
    
    # Perform Wilcoxon test for each pair of clusters
    pvals <- sapply(pairs, function(p) {
      x <- vdf$vaf[vdf$cl == p[1]]
      y <- vdf$vaf[vdf$cl == p[2]]
      if (length(x) >= 1 && length(y) >= 1)
        suppressWarnings(stats::wilcox.test(x, y)$p.value)
      else 1
    })
    
    # If all clusters are significantly different, stop merging
    if (all(pvals <= alpha, na.rm = TRUE)) break
    
    mp <- pairs[[which.max(pvals)]]
    
    vdf$cl[vdf$cl == mp[2]] <- mp[1]
    
    vdf$cl <- as.integer(factor(vdf$cl))
  }
  vdf$cl
}

# Core function: Estimate purity for smple

estimate_purity_for_sample <- function(snv_s, cnv_s, maxAF = 0.55, minPts = 2) {
  
  # Clean NA chromosomes
  snv_s <- snv_s %>% filter(!is.na(chrom))
  cnv_s <- cnv_s %>% filter(!is.na(Chrom))
  
  # Skip sample entirely if no SNV or CNV data
  if (nrow(snv_s) == 0 || nrow(cnv_s) == 0) {
    return(list(purity = NA_real_, n_neutral = 0, eps = NA_real_, n_clusters = 0))
  }
  
  # Build genomic ranges for SNVs and CNV segments ---
  # These allow us to identify which CNV segment each SNV belongs to
  snv_gr <- GRanges(
    seqnames = snv_s$chrom,
    ranges = IRanges(start = snv_s$start, end = snv_s$start)
  )
  cnv_gr <- GRanges(
    seqnames = cnv_s$Chrom,
    ranges = IRanges(start = cnv_s$Start, end = cnv_s$End),
    Tcn = cnv_s$value
  )
  
  # Map SNVs to CNV segments
  hits <- findOverlaps(snv_gr, cnv_gr, select = "first")
  snv_s$seg_idx <- hits
  
  # Keep only SNVs in copy-number-neutral regions (Tcn == 2)
  snv_s <- snv_s %>%
    filter(!is.na(seg_idx)) %>%
    mutate(Tcn = cnv_s$value[seg_idx]) %>%
    filter(Tcn == 2) %>%
    filter(!is.na(dna_vaf), dna_vaf > 0, dna_vaf < maxAF)
  
  n_neutral <- nrow(snv_s)
  
  # Skip sample if too few SNVs remain for meaningful clustering
  if (n_neutral < 3) {
    return(list(purity = NA_real_, n_neutral = n_neutral,
                eps = NA_real_, n_clusters = 0))
  }
  
  # Cluster variant allele frequencies (VAFs) using DBSCAN ---
  v <- snv_s$dna_vaf
  
  
  # Automatically select DBSCAN eps via k-distance method
  eps <- suppressWarnings(pick_eps_kneedle(v, k = minPts))
  
  # Fallback if eps failed or is too small
  if (is.na(eps) || eps <= 0) {
    eps <- stats::median(FNN::knn.dist(matrix(v, ncol = 1), k = 1)[, 1])
  }
  
  # Perform DBSCAN clustering
  db <- tryCatch(
    dbscan::dbscan(matrix(v, ncol = 1), eps = eps, minPts = minPts),
    error = function(e) return(NULL)
  )
  
  # Skip sample if DBSCAN fails
  if (is.null(db)) {
    return(list(purity = NA_real_, n_neutral = n_neutral,
                eps = eps, n_clusters = 0))
  }
  
  # Handle noise points (cluster 0)
  cl <- db$cluster
  if (any(cl == 0L)) {
    nz <- which(cl == 0L)
    # Assign noise points to new unique cluster IDs
    cl[nz] <- max(cl) + seq_along(nz)
  }
  
  # Merge statistically similar clusters (Wilcoxon p > 0.05) ---
  cl <- merge_similar_clusters(v, cl, alpha = 0.05)
  
  # Identify main (clonal) cluster and calculate purity ---
  means <- tapply(v, cl, mean)
  mean_vaf <- means[[as.character(as.integer(names(which.max(means))))]]
  
  # Compute tumor purity
  purity <- max(0, min(1, 2 * mean_vaf))
  
  # Return results for this sample
  list(purity = purity,
       n_neutral = n_neutral,
       eps = eps,
       n_clusters = length(unique(cl)))
}


# Applies the above per-sample function to the entire cohort

estimate_purity_cohort <- function(snv_df, cnv_df, maxAF = 0.60, minPts = 2) {
  
  # Identify overlapping sample IDs in both SNV and CNV tables
  samples <- intersect(unique(snv_df$sample), unique(cnv_df$sample))
  
  # Apply purity estimation sample-by-sample
  res <- purrr::map_dfr(samples, function(s) {
    out <- estimate_purity_for_sample(
      snv_s = snv_df %>% filter(sample == s),
      cnv_s = cnv_df %>% filter(sample == s),
      maxAF = maxAF, minPts = minPts
    )
    tibble(sample = s,
           purity = out$purity,
           n_neutral_snvs = out$n_neutral,
           eps = out$eps,
           n_clusters = out$n_clusters)
  })
  
  # Compute median purity across all valid samples
  med <- median(res$purity[!is.na(res$purity)], na.rm = TRUE)
  
  # Replace missing purity values (no neutral SNVs) with cohort median
  res %>% mutate(purity = ifelse(is.na(purity), med, purity))
}


# Run purity estimation and export results

# Run the cohort-level purity estimation
purity_tbl <- estimate_purity_cohort(snv_df, cnv_df)

# Save in a tsv
write_tsv(purity_tbl, "TARGET-OS.tpes_purity.tsv")

purity <- read_tsv("TARGET-OS.tpes_purity.tsv")


# Plot histogram displaying distribution
library(ggplot2)
library(dplyr)

low_thresh  <- 0.40
high_thresh <- 0.60

# Categorise with explicit ordering per thresholds
purity_tbl <- purity_tbl %>%
  mutate(
    purity_tier = case_when(
      purity < low_thresh ~ "Low (<0.40)",
      purity < high_thresh ~ "Moderate (0.40–0.60)",
      TRUE ~ "High (≥0.60)"
    ),
    purity_tier = factor(purity_tier,
                         levels = c("Low (<0.40)", "Moderate (0.40–0.60)", "High (≥0.60)"))
  )

# Counts per tier
tier_counts <- purity_tbl %>%
  count(purity_tier, .drop = FALSE) %>%
  mutate(label = paste0(as.character(purity_tier), "\n(n = ", n, ")"))

# Compute label positions per tier (centers of each region)
rng_max <- min(1, max(purity_tbl$purity, na.rm = TRUE))
label_pos <- tibble(
  purity_tier = factor(c("Low (<0.40)", "Moderate (0.40–0.60)", "High (≥0.60)"),
                       levels = levels(purity_tbl$purity_tier)),
  x = c(low_thresh/2, (low_thresh + high_thresh)/2, (high_thresh + rng_max)/2)
)

# Robust y position at 90% of density peak (handles any distribution shape)
dens <- density(purity_tbl$purity[is.finite(purity_tbl$purity)], na.rm = TRUE)
y_lab <- 0.90 * max(dens$y)

tier_ann <- tier_counts %>% left_join(label_pos, by = "purity_tier")

ggplot(purity_tbl, aes(x = purity)) +
  annotate("rect", xmin = -Inf, xmax = low_thresh,  ymin = 0, ymax = Inf, alpha = 0.15, fill = "#fc9272") +
  annotate("rect", xmin = low_thresh, xmax = high_thresh, ymin = 0, ymax = Inf, alpha = 0.10, fill = "#fee08b") +
  annotate("rect", xmin = high_thresh, xmax = Inf, ymin = 0, ymax = Inf, alpha = 0.10, fill = "#a1d99b") +
  geom_histogram(aes(y = ..density..), bins = 25, fill = "grey85", color = "white") +
  geom_density(color = "black", size = 1) +
  geom_vline(xintercept = low_thresh,  color = "#d73027", linetype = "dashed", size = 1.1) +
  geom_vline(xintercept = high_thresh, color = "#1a9850", linetype = "dashed", size = 1.1) +
  geom_text(
    data = tier_ann,
    aes(x = x, y = y_lab, label = label),
    color = "black", size = 4.5, fontface = "bold", hjust = 0.5
  ) +
  labs(
    title = "Tumor Purity Distribution — TARGET-OS Cohort",
    subtitle = "TPES (Locallo et al., 2019) with purity categories and sample counts",
    x = "Estimated Tumor Purity",
    y = "Density"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 11, color = "grey30"),
    panel.grid.minor = element_blank()
  )
