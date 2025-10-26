remove <- c("TARGET-40-PALWWX-01A", "TARGET-40-PASYUK-01A", "TARGET-40-PALECC-01A", "TARGET-40-PANVJJ-01A")
# Low-purity samples discovered through purity estimation

library(readr)
library (DESeq2)
library (tidyverse)

counts <- read_tsv("TARGET-OS.htseq_counts.tsv.gz")
metadata <- read_tsv("TARGET-OS.somaticmutation_wxs.tsv.gz")
phenotype <- read_tsv("TARGET-OS.clinical.tsv.gz")
purity

keep_samples <- purity$sample

counts_filtered <- counts %>%
  select(
    Ensembl_ID,
    all_of(setdiff(intersect(colnames(counts), keep_samples), remove))
  )

metadata_filtered <- metadata %>%
  filter(sample %in% keep_samples) %>%
  filter(!sample %in% remove)

phenotype_filtered <- phenotype %>%
  filter(sample %in% keep_samples) %>%
  filter(!sample %in% remove)


# Create new metadata based on mutation status
tp53_mutations <- metadata_filtered %>%
  filter(gene == "TP53")

rb1_mutations <- metadata_filtered %>%
  filter(gene == "RB1")

atrx_mutations <- metadata_filtered %>%
  filter(gene == "ATRX")


# Binary encoding of mutation status for TP53, ATRX, RB1
metadata_cov <- metadata_filtered %>%
  mutate(
    TP53_status = ifelse(sample %in% tp53_mutations$sample, 1, 0),
    RB1_status  = ifelse(sample %in% rb1_mutations$sample, 1, 0),
    ATRX_status = ifelse(sample %in% atrx_mutations$sample, 1, 0)
  )



metadata_cov <- metadata_cov %>%
  group_by(sample) %>%
  summarise(
    TP53_status = max(TP53_status),
    RB1_status  = max(RB1_status),
    ATRX_status = max(ATRX_status),
    .groups = "drop"
  )

metadata_cov %>%
  summarise(
    TP53 = sum(TP53_status),
    RB1  = sum(RB1_status),
    ATRX = sum(ATRX_status),
    wildtype = sum(TP53_status == 0 & RB1_status == 0 & ATRX_status == 0)
  )




samples_counts <- colnames(counts_filtered)[-1]     # exclude Ensembl_ID
samples_meta   <- metadata_cov$sample
samples_pheno  <- phenotype_filtered$sample
samples_purity <- purity$sample


shared_all <- Reduce(intersect, list(
  samples_counts, samples_meta, samples_pheno, samples_purity
))


unique_counts   <- setdiff(samples_counts, shared_all)
unique_metadata <- setdiff(samples_meta, shared_all)
unique_pheno    <- setdiff(samples_pheno, shared_all)
unique_purity   <- setdiff(samples_purity, shared_all)


cat("Unique to counts_filtered:\n"); print(unique_counts)
cat("\nUnique to metadata_cov:\n"); print(unique_metadata)
cat("\nUnique to phenotype_filtered:\n"); print(unique_pheno)
cat("\nUnique to purity:\n"); print(unique_purity)


metadata_cov <- metadata_cov %>%
  filter(sample != "TARGET-40-PAPVYW-01A")


metadata_cov %>%
  summarise(
    TP53 = sum(TP53_status),
    RB1  = sum(RB1_status),
    ATRX = sum(ATRX_status),
    wildtype = sum(TP53_status == 0 & RB1_status == 0 & ATRX_status == 0)
  )


counts_samples <- colnames(counts_filtered)[-1]  # skip first column if gene IDs

# Extract sample names from metadata_cov
metadata_samples <- metadata_cov$sample

# Check if all counts columns exist in metadata
all(counts_samples %in% metadata_samples)


metadata_cov <- metadata_cov %>%
  arrange(match(sample, colnames(counts_filtered)[-1]))

counts_filtered_df <- as.data.frame(counts_filtered)


# Remove last 5 rows as they're not count data for genes
counts_filtered_df <- counts_filtered_df[1:(nrow(counts_filtered_df) - 5), ]

count_matrix <- counts_filtered_df[, -1]           # remove gene ID column
rownames(count_matrix) <- counts_filtered_df$Ensembl_ID  # set gene IDs as rownames

count_matrix <- as.matrix(count_matrix)
mode(count_matrix) <- # DESeq2 requires integers

ncol(count_matrix)

write.csv(count_matrix, "count_matrix_integer.csv", row.names = FALSE)

metadata_cov <- metadata_cov %>%
  mutate(
    TP53_status = factor(TP53_status, levels = c(0,1), labels = c("WT","MUT")),
    RB1_status  = factor(RB1_status,  levels = c(0,1), labels = c("WT","MUT")),
    ATRX_status = factor(ATRX_status, levels = c(0,1), labels = c("WT","MUT"))
  )

library(edgeR) # Removing low count genes using EdgeR
y <- DGEList(counts = count_matrix)
keep <- filterByExpr(y, design = model.matrix(~ TP53_status + RB1_status + ATRX_status, data = metadata_cov))
count_matrix_filtered <- count_matrix[keep, ]


library(DESeq2)


#DESeq2 Matrix Setup for DEA and PCA
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix_filtered,
  colData = metadata_cov,
  design = ~ TP53_status + RB1_status + ATRX_status
)



dds <- DESeq(dds)


resultsNames(dds)


res_tp53 <- results(dds, name = "TP53_status_MUT_vs_WT", pAdjustMethod = "BH")
res_rb1  <- results(dds, name = "RB1_status_MUT_vs_WT",  pAdjustMethod = "BH")
res_atrx <- results(dds, name = "ATRX_status_MUT_vs_WT", pAdjustMethod = "BH")



library(tibble)

res_tp53_df <- as.data.frame(res_tp53) %>% rownames_to_column("gene_id")
res_rb1_df  <- as.data.frame(res_rb1)  %>% rownames_to_column("gene_id")
res_atrx_df <- as.data.frame(res_atrx) %>% rownames_to_column("gene_id")




library(ggplot2)
library(dplyr)



# PCA of TARGET-OS VST-normalised counts by specific mutation combinations


library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(matrixStats)

# Convert mutation factors to numeric ----
metadata_cov <- metadata_cov %>%
  mutate(
    TP53_num = ifelse(TP53_status == "MUT", 1, 0),
    RB1_num  = ifelse(RB1_status  == "MUT", 1, 0),
    ATRX_num = ifelse(ATRX_status == "MUT", 1, 0)
  )

# Define mutation combination categories
metadata_cov <- metadata_cov %>%
  mutate(
    MUT_status = case_when(
      TP53_num + RB1_num + ATRX_num == 0 ~ "Wildtype",
      TP53_num == 1 & RB1_num == 0 & ATRX_num == 0 ~ "TP53_only",
      RB1_num  == 1 & TP53_num == 0 & ATRX_num == 0 ~ "RB1_only",
      ATRX_num == 1 & TP53_num == 0 & RB1_num == 0 ~ "ATRX_only",
      TP53_num == 1 & RB1_num == 1 & ATRX_num == 0 ~ "TP53+RB1",
      TP53_num == 1 & ATRX_num == 1 & RB1_num == 0 ~ "TP53+ATRX",
      RB1_num  == 1 & ATRX_num == 1 & TP53_num == 0 ~ "RB1+ATRX",
      TP53_num == 1 & RB1_num == 1 & ATRX_num == 1 ~ "TP53+RB1+ATRX",
      TRUE ~ "Other"
    )
  )

# Variance-stabilising transformation
vsd <- vst(dds, blind = TRUE)

# PCA on top 500 most variable genes
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 500)
vst_mat <- assay(vsd)[topVarGenes, ]

pca <- prcomp(t(vst_mat), center = TRUE, scale. = FALSE)

# Computing Variance for each component
percentVar <- round(100 * (pca$sdev^2 / sum(pca$sdev^2))[1:2], 1)

# Build PCA dataframe for plotting
pca_df <- data.frame(
  Sample = colnames(vsd),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
) %>%
  left_join(metadata_cov %>% select(sample, MUT_status),
            by = c("Sample" = "sample"))


# Plot PCA (no sample names)
ggplot(pca_df, aes(PC1, PC2, color = MUT_status)) +
  geom_point(size = 4, alpha = 0.9) +
  theme_minimal(base_size = 14) +
  labs(
    title = "PCA of TARGET-OS (VST counts, top 500 variable genes)",
    x = paste0("PC1 (", percentVar[1], "% variance)"),
    y = paste0("PC2 (", percentVar[2], "% variance)")
  ) +
  scale_color_manual(values = c(
    "Wildtype" = "grey60",
    "TP53_only" = "#D7191C",
    "RB1_only" = "#2C7BB6",
    "ATRX_only" = "#31A354",
    "TP53+RB1" = "#FDAE61",
    "TP53+ATRX" = "#A6CEE3",
    "RB1+ATRX" = "#B2DF8A",
    "TP53+RB1+ATRX" = "#756BB1",
    "Other" = "black"
  )) +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )


library(ggplot2)
library(dplyr)
library(ggrepel)

# GSEA
# Run first for TP53, then manually for RB1 and ATRX


# Install if not already
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "DOSE"),
                     #force = TRUE
                     #)

# Load libraries
library(clusterProfiler)
library(org.Hs.eg.db)   .
library(enrichplot)
library(DOSE)


# Remove version numbers from Ensembl IDs
gene_ids <- gsub("\\..*", "", rownames(res_tp53))

# Map Ensembl → Entrez using org.Hs.eg.db
library(org.Hs.eg.db)
entrez_ids <- mapIds(org.Hs.eg.db,
                     keys = gene_ids,
                     column = "ENTREZID",
                     keytype = "ENSEMBL",
                     multiVals = "first")


# Add mapped Entrez IDs to your DESeq2 results
res_tp53$entrez <- entrez_ids

# Remove genes without mapping
res_tp53 <- res_tp53[!is.na(res_tp53$entrez), ]

# Create ranked list for GSEA based on log fold change
gene_list <- res_tp53$log2FoldChange
names(gene_list) <- res_tp53$entrez

# Sort in descending order
gene_list <- sort(gene_list, decreasing = TRUE)


library(clusterProfiler) # Package for GO-based GSEA
library(org.Hs.eg.db)
library(enrichplot)

# Create a clean ranked gene list again
gene_list <- res_tp53$log2FoldChange
names(gene_list) <- res_tp53$entrez

# Remove any duplicate names
gene_list <- gene_list[!duplicated(names(gene_list))]

# Sort in descending order
gene_list <- sort(gene_list, decreasing = TRUE)


# GO: Biological Processes
gsea_go <- gseGO(
  geneList = gene_list,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = FALSE
)


library(ggplot2)


# Top 15 enriched GO terms
dotplot(gsea_go, showCategory = 15, font.size = 12, label_format = 40) +
  ggtitle("GSEA – GO Biological Processes (TP53 Mutated)") +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 10)
  )



