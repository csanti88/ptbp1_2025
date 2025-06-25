# Load required libraries
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(ggrepel)
library(corrplot)
library(here)

# ============================
# 1. Load Input Data
# ============================
filewithcounts <- here("data","Bulk-hetvsko_and_all","20250612_combined_counts_with_names.tsv")

# Load counts and gene metadata
counts <- read.delim(filewithcounts, check.names = FALSE)
genes <- counts[, c("Geneid", "GeneName")]

# Prepare count matrix
rownames(counts) <- counts$Geneid
counts <- counts[, !(colnames(counts) %in% c("Geneid", "GeneName"))]

# ============================
# 2. Define Sample Groups
# ============================
group <- factor(c(rep("WT", 5), rep("Het", 3), rep("KO", 3), rep("Cell", 8)))
colData <- data.frame(row.names = colnames(counts), group = group)

# ============================
# 3. Run DESeq2
# ============================
dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ group)
dds <- dds[rowSums(counts(dds) >= 6) >= 3, ]
dds <- DESeq(dds)

# ============================
# 4. rlog Transformation
# ============================
rld <- rlog(dds, blind = TRUE)
rld_mat <- assay(rld)
rld_df <- as.data.frame(rld_mat)
rld_df$Geneid <- rownames(rld_df)
rlog_counts <- merge(rld_df, genes, by = "Geneid")

# ============================
# 5. Heatmap of Selected Genes
# ============================
selected_genes <- c("Ptbp1", "Ptbp2", "Fgf15", "Nr2f1", "Foxp1", "Elavl4",
                    "Sfrp2", "Hes1", "Atoh7", "Pou4f2", "Meis2", "Sox8",
                    "Hes5", "Nfix", "Opn1sw", "Rcvrn", "Gfap", "Id1")

heatmap_df <- rlog_counts[match(selected_genes, rlog_counts$GeneName), ]
rownames(heatmap_df) <- heatmap_df$GeneName

# Reorder columns manually
ordered_cols <- c("Geneid", "GeneName",
                  "wtP28_1", "wtP28_2", "bam0_1_wt_1",
                  "bam1_1_het_1", "bam1_2_het_2", "bam1_3_het_3",
                  "bam2_1_ko_1", "bam2_2_ko_2", "bam2_3_ko_3",
                  "wtE16_1", "wtE16_2",
                  "rod", "thy1_1",
                  "astro1", "astro2",
                  "microglia1", "microglia2", "hair_cell", "heart")
heatmap_df <- heatmap_df[, ordered_cols]

heatmap_labels <- c("WT E16 1", "WT E16 2", "WT1*",
                    "Het1*", "Het2*", "Het3*",
                    "KO1*", "KO2*", "KO3*",
                    "WT P28 1", "WT P28 2",
                    "Rod", "Thy1 Neuron", "Astrocyte1", "Astrocyte2",
                    "Microglia1", "Microglia2", "Hair cell", "Cardiomyocyte")

p<-pheatmap(as.matrix(heatmap_df[3:21]),
         color = viridis::viridis(100),
         labels_col = heatmap_labels,
         cluster_rows = FALSE, cluster_cols = FALSE,
         gaps_row = 2, gaps_col = 12,
         scale = "row", legend = TRUE,
         fontsize_row = 10, border_color = NA)
p

# Ensure the results folder exists
if (!dir.exists(here("results"))) dir.create(here("results"))

# Save PCA plot as PDF
ggsave(
  filename = here("results", "HeatMapDGE_alltissue.pdf"),
  plot = p,
  width = 8,
  height = 6,
  units = "in"
)

# ============================
# 6. PCA and Correlation Analysis
# ============================

de_file<-here("data","Bulk-hetvsko_and_all",
              "combined_counts_with_names_7samples.tsv2025_06_2320250612DESeq2_results_with_symbols_log2.1.2_pvalue0.05.csv")
DE_genes <- read.csv(de_file)
ids_DE <- DE_genes$Geneid

# Subset expression matrix
corr_pca_df <- rlog_counts[match(ids_DE, rlog_counts$Geneid), ]
corr_pca_df <- corr_pca_df[, ordered_cols]

# Relabel columns
colnames(corr_pca_df) <- c("Geneid", "GeneName", heatmap_labels)

# Run PCA
pca_mat <- t(corr_pca_df[, 3:21])
pca_res <- prcomp(pca_mat, center = TRUE, scale. = TRUE)
pca_df <- as.data.frame(pca_res$x)
pca_df$Sample <- rownames(pca_df)

# Define group coloring
pca_df$Group <- factor(c("WT", "WT", "WT_1",
                         "Het", "Het", "Het",
                         "KO", "KO", "KO",
                         "WT28", "WT28",
                         "Rod", "Thy1 Neuron",
                         "Astrocyte", "Astrocyte",
                         "Microglia", "Microglia",
                         "Hair cell", "Cardiomyocyte"))

# Plot PCA
pca_plot<-ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3, max.overlaps = 100) +
  labs(title = "PCA Gene Explression for all samples",
       x = paste0("PC1 (", round(100 * summary(pca_res)$importance[2,1], 1), "%)"),
       y = paste0("PC2 (", round(100 * summary(pca_res)$importance[2,2], 1), "%)")) +
  theme_minimal()

# Ensure the results folder exists
if (!dir.exists(here("results"))) dir.create(here("results"))

# Save PCA plot as PDF
ggsave(
  filename = here("results", "PCA_geneexpression_all.pdf"),
  plot = pca_plot,
  width = 8,
  height = 6,
  units = "in"
)


# Correlation matrix
expr_data <- corr_pca_df[, 3:21] %>%
  lapply(as.numeric) %>%
  as.data.frame()
rownames(expr_data) <- corr_pca_df$Geneid

cor_mat <- cor(expr_data, method = "pearson")

# Heatmap of correlations
pheatmap(cor_mat,
         color = viridis::viridis(20),
         border_color = 'white',
         display_numbers = TRUE,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Sample-to-Sample Correlation")

# Open PDF device
pdf(file = here("results", "correlation_plot.pdf"), width = 7, height = 6)

# Generate corrplot
corrplot(cor_mat, method = "color", type = "lower",
         tl.col = "black", addCoef.col = "black",
         number.cex = 0.6,
         col = c(rep("white", 140), viridis::viridis(60)),
         col.lim = c(0.4, 1))

# Close PDF device
dev.off()
