# Load libraries
library(dplyr)
library(ggplot2)
library(ggrepel)
library(tidyr)
library(purrr)
library(stringr)
library(here)

# Define project paths using `here`
rmats_all_dir <- here("data", "Rmats_all_tissue")
rmats_wte16p28_dir <- here("data", "Rmats_wte16p28_relaxfiltering")

#function
load_eventdf_data <- function(event, base_dir) {
  event_file <- file.path(base_dir, paste0(event, ".MATS.JC.txt"))
  event_df <- read.table(event_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  return(event_df)
}

# Load and process rMATS data for all samples
splicejunctiontypes <- c("A5SS", "A3SS", "RI", "MXE", "SE")
fileevents <- c("up_", "dn_", "bg_", "filtered_")

all_rmats_dfs <- lapply(splicejunctiontypes,
                               function(splicejunctiontype) load_eventdf_data(paste0("filtered_",
                                                                                     splicejunctiontype),
                                                                              rmats_all_dir))

names(all_rmats_dfs) <- fileevents
walk(all_rmats_dfs, ~print(dim(.x)))


all_rmats_combined<-lapply(all_rmats_dfs, function(df) {
  df[, c("ID","IncLevelDifference", "FDR", "geneSymbol", "IncLevel1", "IncLevel2")]
})
all_rmats_combined <- do.call(rbind, all_rmats_combined)

# Define sample groupings
samples_group1 <- c("thy1_1", "thy1_2", "wtE16_1", "wtE16_2", "wtP28_1", "wtP28_2", "hairy", "wt")
samples_group2 <- c("heart", "rod", "astro1", "astro2", "microglia1", "microglia2", "het1", "het2", "het3", "ko1", "ko2", "ko3")

# Expand inclusion levels
df_all_expanded <- all_rmats_combined %>%
  mutate(across(starts_with("IncLevel"), ~str_split(as.character(.), ","))) %>%
  rowwise() %>%
  mutate(
    IncLevel1 = list({ v <- IncLevel1; length(v) <- 8; v }),
    IncLevel2 = list({ v <- IncLevel2; length(v) <- 12; v })
  ) %>%
  ungroup() %>%
  unnest_wider(IncLevel1, names_sep = "_") %>%
  unnest_wider(IncLevel2, names_sep = "_")

colnames(df_all_expanded)[grep("IncLevel1_", colnames(df_all_expanded))] <- samples_group1
colnames(df_all_expanded)[grep("IncLevel2_", colnames(df_all_expanded))] <- samples_group2

# Filter based on differential splicing events
splicing_events_ids_file <- here("results","filteredPSIhetvsko.csv")
splicing_events_ids_file<-read.csv(splicing_events_ids_file)
splicing_events_ids<-splicing_events_ids_file$ID
splicing_events_all <- df_all_expanded %>% filter(ID %in% splicing_events_ids)

# Load and process wtE16 vs wtP28 data


wt16p28_dfs <- lapply(splicejunctiontypes,
                      function(splicejunctiontype) load_eventdf_data(paste0("filtered_",
                                                                            splicejunctiontype),
                                                                     rmats_wte16p28_dir))
names(wt16p28_dfs) <- fileevents
walk(wt16p28_dfs, ~print(dim(.x)))

wt16p28_combined<-lapply(wt16p28_dfs, function(df) {
  df[, c("ID","IncLevelDifference", "FDR", "geneSymbol", "IncLevel1", "IncLevel2")]
})
wt16p28_combined <- do.call(rbind, wt16p28_combined)

wt16p28_expanded <- wt16p28_combined %>%
  mutate(across(starts_with("IncLevel"), ~str_split(as.character(.), ","))) %>%
  rowwise() %>%
  mutate(
    IncLevel1 = list({ v <- IncLevel1; length(v) <- 1; v }),
    IncLevel2 = list({ v <- IncLevel2; length(v) <- 1; v })
  ) %>%
  ungroup() %>%
  unnest_wider(IncLevel1, names_sep = "_") %>%
  unnest_wider(IncLevel2, names_sep = "_")

colnames(wt16p28_expanded)[grep("IncLevel1_", colnames(wt16p28_expanded))] <- "wtE16"
colnames(wt16p28_expanded)[grep("IncLevel2_", colnames(wt16p28_expanded))] <- "wtP28"

splicing_events_wt16p28 <- wt16p28_expanded %>% filter(ID %in% splicing_events_ids)

# Merge both datasets
merged_events <- splicing_events_all %>%
  select(-wtE16_1, -wtE16_2, -wtP28_1, -wtP28_2, -rod) %>%
  left_join(splicing_events_wt16p28 %>% select(ID, wtE16, wtP28), by = "ID")

# Rename and prepare for PCA
colnames(merged_events) <- c("ID", "IncLevelDifference", "FDR", "geneSymbol",
                             "Thy1_1", "Thy1_2", "HairCell", "Wt", "Heart",
                             "Astrocyte1", "Astrocyte2", "Microglia1", "Microglia2",
                             "Het1", "Het2", "Het3", "Ko1", "Ko2", "Ko3", "WtE16", "WtP28")

merged_clean <- merged_events %>% select(-Thy1_2)
sample_cols <- colnames(merged_clean)[5:ncol(merged_clean)]
group_labels <- c("Thy1", "Hair", "WT", "Heart", "Astrocyte", "Astrocyte", "Microglia", "Microglia",
                  rep("Het", 3), rep("KO", 3), rep("WT", 2))

data_pca <- merged_clean %>%
  select(all_of(sample_cols)) %>%
  mutate(across(everything(), as.numeric)) %>%
  na.omit()

pca_result <- prcomp(t(data_pca), scale. = TRUE)
pca_df <- as.data.frame(pca_result$x)
pca_df$Sample <- rownames(pca_df)
pca_df$Group <- group_labels[sample_cols %in% colnames(data_pca)]

# Plot PCA
p <- ggplot(pca_df, aes(x = PC1, y = PC2, label = Sample, color = Group)) +
  geom_point(size = 2) +
  geom_text_repel(size = 4, max.overlaps = 100, box.padding = 0.4) +
  labs(title = "PCA of PSI values (NAs removed)", x = "PC1", y = "PC2") +
  theme_minimal() +
  theme(legend.position = "bottom")
print(p)

# Result
data_pca

# Ensure the results folder exists
if (!dir.exists(here("results"))) dir.create(here("results"))

# Save PCA plot as PDF
ggsave(
  filename = here("results", "PCA_PSI_all_plot.pdf"),
  plot = p,
  width = 8,
  height = 6,
  units = "in"
)
