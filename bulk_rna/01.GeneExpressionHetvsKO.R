# Load required libraries
library(DESeq2)
library(tidyverse)
library(dplyr)
library(here)

# ============================
# 1. Load Input Data
# ============================
filewithcounts <- here("data","Bulk-hetvsko","combined_counts_with_names_7samples.tsv")

counts <- read.delim(filewithcounts, check.names = FALSE)
genes <- counts[, c("Geneid", "GeneName")]

# Prepare count matrix
rownames(counts) <- counts$Geneid
counts <- counts[, !(colnames(counts) %in% c("Geneid", "GeneName"))]

# ============================
# 2. Sample Metadata
# ============================
samplenames <- substring(colnames(counts), 8)
group <- factor(c("WT", rep("Het", 3), rep("KO", 3)))

colData <- data.frame(row.names = colnames(counts), group = group)

# ============================
# 3. DESeq2 Analysis
# ============================
dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ group)
dds <- dds[rowSums(counts(dds) >= 6) >= 3, ]  # Filter low count genes
dds <- DESeq(dds)
res <- results(dds, contrast = c("group", "KO", "Het"))

# ============================
# 4. Normalized & Raw Counts
# ============================
norm_counts_df <- as.data.frame(counts(dds, normalized = TRUE))
norm_counts_df$Geneid <- rownames(norm_counts_df)
raw_counts_df <- as.data.frame(counts(dds))
raw_counts_df$Geneid <- rownames(raw_counts_df)

# DE Results
res_df <- as.data.frame(res)
res_df$Geneid <- rownames(res_df)

# Merge all data
merged_df <- norm_counts_df %>%
  left_join(raw_counts_df, by = "Geneid", suffix = c(".n", ".r")) %>%
  left_join(res_df, by = "Geneid") %>%
  left_join(genes, by = "Geneid")

# ============================
# 5. Filter Significant Genes
# ============================
merged_df_clean <- merged_df %>%
  filter(!is.na(log2FoldChange) & !is.na(padj))

merged_df_filtered <- merged_df_clean %>%
  filter(abs(log2FoldChange) >= log2(1.2), padj < 0.05)

# Save filtered results
#output_file <- paste0(filewithcounts, format(Sys.time(), "%Y_%m_%d"), "_DESeq2_filtered.csv")

write.csv(merged_df_filtered, here("results","HetvsKO_DESeq2_filtered.csv"), row.names = FALSE)

# ============================
# 6. Top Up & Down Regulated Genes
# ============================
filter_pseudo <- function(df) {
  df %>%
    filter(!str_detect(GeneName, "ps|\\-ps|\\-rs")) %>%
    filter(!str_starts(GeneName, "Gm"),
           !str_starts(GeneName, "RP"),
           !str_ends(GeneName, "Rik"))
}

# Top upregulated
top_up <- merged_df_filtered %>% arrange(desc(log2FoldChange))
top_up_clean <- filter_pseudo(top_up)
top15up <- top_up_clean$GeneName[1:15]
top15up_withpseudo <- top_up$GeneName[1:15]

# Top downregulated
top_down <- merged_df_filtered %>% arrange(log2FoldChange)
top_down_clean <- filter_pseudo(top_down)
top15down <- top_down_clean$GeneName[1:15]
top15down_withpseudo <- top_down$GeneName[1:15]

plotgenes <- c(top15up, top15down, "Ptbp1", "Ptbp2")
plotgenesnf <- c(top15up_withpseudo, top15down_withpseudo, "Ptbp1", "Ptbp2")

# ============================
# 7. Volcano Plot
# ============================
othergenes <- c("Pde6b", "Rcvrn", "Nr2e3", "Rdh8", "Gngt1", "Samd7", "Tfap2e", "Pitx2",
                "Apoa1", "Rpl29", "Foxp4", "Ptbp1", "Ptbp2", "Foxe3", "H2-Ab1", "Atp6v0c")

p<-EnhancedVolcano(merged_df_clean,
                lab = merged_df_clean$GeneName,
                x = 'log2FoldChange',
                y = 'padj',
                title = 'KO vs Het',
                selectLab = othergenes,
                pCutoff = 0.05,
                FCcutoff = log2(1.2),
                max.overlaps = 100,
                pointSize = 1,
                labSize = 6,
                legendPosition = 'bottom',
                legendLabSize = 17,
                legendIconSize = 4.0,
                boxedLabels = FALSE,
                drawConnectors = TRUE,
                widthConnectors = 1,
                lengthConnectors = 30,
                colConnectors = 'black',
                arrowheads = FALSE,
                gridlines.major = TRUE,
                gridlines.minor = FALSE,
                col = c('gray', 'gray', 'gray', 'steelblue1'),
                ylim = c(0, 100),
                xlim = c(-10, 10)) +
  ylab("-Log10 P-adj")

p

# Ensure the results folder exists
if (!dir.exists(here("results"))) dir.create(here("results"))

# Save PCA plot as PDF
ggsave(
  filename = here("results", "DGE_Het_vs_KO_plot.pdf"),
  plot = p,
  width = 8,
  height = 6,
  units = "in"
)