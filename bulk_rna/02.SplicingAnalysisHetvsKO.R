# ===============================
# Setup: Load Libraries & Set Paths
# ===============================
library(EnhancedVolcano)
library(ggplot2)
library(here)

# Set base directory
base_dir <- here("data")
test_dir <- file.path("data", "Rmats_ab1_het_vs_b2_ko")
test_dir
# ===============================
# Pie Chart Function for Splicing Event Types
# ===============================
count_lines <- function(event, base_dir) {
  up_file <- file.path(base_dir, paste0("up_", event, ".MATS.JC.txt"))
  dn_file <- file.path(base_dir, paste0("dn_", event, ".MATS.JC.txt"))
  
  n_up <- length(readLines(up_file))
  n_dn <- length(readLines(dn_file))
  
  return(n_up + n_dn - 2)  # subtract 2 for headers
}

piechart_function <- function(base_dir) {
  Prop <- c(
    SE = count_lines("SE", base_dir),
    A5SS = count_lines("A5SS", base_dir),
    A3SS = count_lines("A3SS", base_dir),
    MXE = count_lines("MXE", base_dir),
    RI = count_lines("RI", base_dir)
  )
  
  labels <- paste0(names(Prop), ' ', format(Prop, big.mark = ""))
  title <- basename(base_dir)
  
  # Save to PDF
  pdf_path <- file.path(base_dir, paste0("pie_", title, ".pdf"))
  pdf(pdf_path, width = 4, height = 4)
  pie(Prop, labels = labels, clockwise = TRUE, angle = 0, main = title)
  dev.off()
  
  return(pdf_path)
}

# Run pie chart function
piepath <- piechart_function(test_dir)

# ===============================
# Load and Combine Up/Down Events
# ===============================
load_up_dn_data <- function(event, base_dir) {
  up_file <- file.path(base_dir, paste0("up_", event, ".MATS.JC.txt"))
  dn_file <- file.path(base_dir, paste0("dn_", event, ".MATS.JC.txt"))
  
  df_up <- read.table(up_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  df_dn <- read.table(dn_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  
  df_up$direction <- "up"
  df_dn$direction <- "down"
  
  return(rbind(df_up, df_dn))
}

# Load data for each splicing type
df_SE   <- load_up_dn_data("SE", test_dir)
df_A5SS <- load_up_dn_data("A5SS", test_dir)
df_A3SS <- load_up_dn_data("A3SS", test_dir)
df_MXE  <- load_up_dn_data("MXE", test_dir)
df_RI   <- load_up_dn_data("RI", test_dir)

# ===============================
# Bar Plot Summary of Splicing Changes
# ===============================
values <- c(
  RI = nrow(df_RI),
  A5SS = nrow(df_A5SS),
  A3SS = nrow(df_A3SS),
  MXE = nrow(df_MXE),
  SE = nrow(df_SE)
)

labels <- paste0(names(values), ' ', format(values, big.mark = ""))
bar_data <- data.frame(
  group = "Het vs KO",
  SplicingType = factor(labels, levels = labels),
  count = values
)

p<-ggplot(bar_data, aes(x = group, y = count, fill = SplicingType)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Splicing Changes: Het vs KO",
    x = "Group", y = "# Splicing Events"
  ) +
  scale_fill_manual(values = c("khaki1", "#BDC3C7", "grey50", "seagreen1", "#3498DB")) +
  theme_minimal()
print(p)

# Ensure the results folder exists
if (!dir.exists(here("results"))) dir.create(here("results"))

# Save PCA plot as PDF
ggsave(
  filename = here("results", "BarPSI_hetvsKO.pdf"),
  plot = p,
  width = 3,
  height = 6,
  units = "in"
)

# ===============================
# Combine All Filtered Junctions
# ===============================
splicejunctiontypes <- c("A5SS", "A3SS", "RI", "MXE", "SE")

load_eventdf_data <- function(file_prefix, base_dir) {
  file_path <- file.path(base_dir, paste0(file_prefix, ".MATS.JC.txt"))
  read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
}

ab1_het_vs_b2_ko_dfs <- lapply(splicejunctiontypes, function(type) {
  load_eventdf_data(paste0("filtered_", type), test_dir)
})
names(ab1_het_vs_b2_ko_dfs) <- splicejunctiontypes

# Extract relevant columns
subset_ab1_dfs <- lapply(ab1_het_vs_b2_ko_dfs, function(df) {
  df[, c("ID", "IncLevelDifference", "FDR", "geneSymbol", "IncLevel1", "IncLevel2")]
})
combined_df <- do.call(rbind, subset_ab1_dfs)




# ===============================
# Extract Top Spliced Genes
# ===============================
extract_and_label <- function(df, event_type) {
  df_subset <- df[, c("ID", "IncLevelDifference", "FDR", "geneSymbol", "direction", "IncLevel1", "IncLevel2")]
  df_subset$eventType <- event_type
  return(df_subset)
}

# Combine all splicing types
df_all <- rbind(
  extract_and_label(df_SE, "SE"),
  extract_and_label(df_A5SS, "A5SS"),
  extract_and_label(df_A3SS, "A3SS"),
  extract_and_label(df_MXE, "MXE"),
  extract_and_label(df_RI, "RI")
)

write.csv(df_all,here("results","filteredPSIhetvsko.csv"),row.names = FALSE)
write.csv(df_all,here("data","filteredPSIhetvsko.csv"),row.names = FALSE)
# Top positive and negative PSI changes
top_pos_unique <- df_all[order(-df_all$IncLevelDifference), ]
top_pos_unique <- top_pos_unique[!duplicated(top_pos_unique$geneSymbol), ][1:30, ]

top_neg_unique <- df_all[order(df_all$IncLevelDifference), ]
top_neg_unique <- top_neg_unique[!duplicated(top_neg_unique$geneSymbol), ][1:25, ]

# ===============================
# Volcano Plot for Manual Junctions
# ===============================

up_ids <-c("chr9:106444089-106444306|-|106440856|106447530",
      "chr17:22224010-22224094|-|22205352|22225538",
      "chr17:6118295-6118366|+|6107788|6121194",
      "chr4:3569630-3569752|-|3569463|3574703",
      "chr11:84533382-84533456|+|84525748|84534479",
      "chr14:60706713-60707079|+|60692548|60709424"
)
down_ids<-c("chr9:106445525-106445829|-|106444306|106447530",
        "chr16:45774162-45774345|-|45772275|45774943",
        "chr8:69744191-69744252|-|69743610|69744460",
        "chr15:99737546-99737653|-|99736436|99738664",
        "chr1:80503918-80503994|-|80501826|80505324",
        "chr12:73293698-73293799|+|73292233|73310144"
)

selected_ids <- c(up_ids, down_ids)

custom_labels <- ifelse(
  combined_df$ID %in% selected_ids,
  combined_df$geneSymbol,
  ""
)


p<-EnhancedVolcano(combined_df,
                lab = custom_labels,
                x = 'IncLevelDifference',
                y = 'FDR',
                pCutoff = 0.05,
                FCcutoff = 0.1,
                pointSize = 1,
                labSize = 7,
                legendPosition = 'bottom',
                legendLabSize = 5,
                legendIconSize = 3,
                drawConnectors = TRUE,
                max.overlaps = 1000,
                col = c('grey', 'grey', 'grey', '#63B8FF'),
                colAlpha = 1,
                widthConnectors = 0.5,
                xlim = c(-1, 1),
                xlab = "PSI difference",
                ylab = "-log(FDR)"
)
p

# Ensure the results folder exists
if (!dir.exists(here("results"))) dir.create(here("results"))

# Save PCA plot as PDF
ggsave(
  filename = here("results", "VolcanoPSI_hetvsKO.pdf"),
  plot = p,
  width = 6,
  height = 8,
  units = "in"
)

