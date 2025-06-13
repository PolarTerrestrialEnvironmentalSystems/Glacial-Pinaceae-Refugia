#=============================================================================
# Script: 4.2_PLOT_READ_COUNT.R
# Description: Plot Kraken2-assigned Pinaceae read counts for each site.
#   - Part 1: Saves separate plots for Pinaceae and species (uniform & variable Y).
#   - Part 2: Saves combined plot with both uniform & variable Y-axis species plots.
#=============================================================================

# Load required libraries
library(tidyverse)
library(viridis)
library(cowplot)
library(stringr)

# Set base directory (update this path to your setup)
base_dir <- "path-to-your-project"

# Load site name mapping file (folder_name to new_name)
mapping_file_path <- file.path(base_dir, "results/6_7_read_counts", "new_names.txt")
mapping_file <- read.table(
  mapping_file_path, header = TRUE, sep = "\t",
  stringsAsFactors = FALSE, fileEncoding = "UTF-8", fill = TRUE, strip.white = TRUE
)
if (ncol(mapping_file) < 2) {
  mapping_file <- read.table(
    mapping_file_path, header = TRUE, sep = "",
    stringsAsFactors = FALSE, fileEncoding = "UTF-8", fill = TRUE, strip.white = TRUE
  )
}
colnames(mapping_file) <- c("folder_name", "new_name")
mapping_file <- mapping_file[!is.na(mapping_file$folder_name) & !is.na(mapping_file$new_name), ]
site_name_mapping <- setNames(mapping_file$new_name, mapping_file$folder_name)
sites <- as.character(mapping_file$folder_name)

# Define function to process each site
process_site <- function(site_name, site_name_mapping) {
  site_name <- as.character(site_name)
  new_site_name <- site_name_mapping[site_name]
  if (is.na(new_site_name)) new_site_name <- site_name
  
  pinaceae_file <- file.path(base_dir, "results/6_7_read_counts", site_name, "reads_count_pinaceae.csv")
  species_file <- file.path(base_dir, "results/6_7_read_counts", site_name, "reads_count_species.csv")
  meta_file <- file.path(base_dir, "scripts/Pipeline/rename_samples", site_name, paste0(site_name, "_meta_file.txt"))
  if (!file.exists(pinaceae_file) || !file.exists(species_file) || !file.exists(meta_file)) return(NULL)
  
  original_data <- read_csv(pinaceae_file, show_col_types = FALSE)
  colnames(original_data) <- c("ID", "reads")
  merged_data <- original_data %>%
    group_by(ID = sub("_(.+)", "", ID)) %>%
    summarise(Pinaceae = sum(reads)) %>%
    ungroup()
  
  meta <- read.delim(meta_file, header = FALSE, stringsAsFactors = FALSE)
  colnames(meta) <- c("file", "library", "age", "lake", "name")
  merged_data <- merge(merged_data, meta, by.x = "ID", by.y = "library", all.x = TRUE) %>%
    select(-file, -name)
  
  species_data <- read.csv(species_file, stringsAsFactors = FALSE)
  species_data <- species_data %>%
    mutate(
      numeric_range = str_extract(BAM.File, "-\d+-\d+\\.bam$"),
      ID_full = str_extract(BAM.File, "^(.*?)_combined"),
      ID = sub("^(.*?)_(.*?)_combined$", "\\1", ID_full),
      Species = case_when(
        numeric_range == "-1-122588.bam" ~ "Larix",
        numeric_range == "-122589-239988.bam" ~ "Pinus",
        numeric_range == "-239989-361212.bam" ~ "Abies",
        numeric_range == "-361213-485324.bam" ~ "Picea",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(Species)) %>%
    select(ID, Species, Mapped.Reads) %>%
    group_by(ID, Species) %>%
    summarise(Mapped.Reads = sum(Mapped.Reads), .groups = 'drop') %>%
    pivot_wider(names_from = Species, values_from = Mapped.Reads, values_fill = 0)
  
  final_data <- left_join(merged_data, species_data, by = "ID")
  for (s in c("Larix", "Pinus", "Abies", "Picea")) if (!s %in% colnames(final_data)) final_data[[s]] <- 0
  final_data[is.na(final_data)] <- 0
  final_data <- final_data %>%
    mutate(Pinaceae_sum = Larix + Pinus + Abies + Picea) %>%
    mutate(across(c(Larix, Pinus, Abies, Picea), ~round(.x / Pinaceae_sum, 2), .names = "{.col}_rel"))
  
  write.csv(final_data, file.path(base_dir, "results/6_7_read_counts", site_name, paste0(site_name, "_final_merged_data.csv")), row.names = FALSE)
  
  x_var <- ifelse(site_name == "ximen", "ID", "age")
  final_data <- final_data %>%
    mutate(across(all_of(x_var), as.character)) %>%
    mutate(across(all_of(x_var), ~factor(.x, levels = unique(.x))))
  
  long_data <- final_data %>%
    select(all_of(x_var), Pinaceae, Larix, Pinus, Abies, Picea) %>%
    pivot_longer(-all_of(x_var), names_to = "variable", values_to = "value") %>%
    mutate(variable = factor(variable, levels = c("Pinaceae", "Picea", "Pinus", "Larix", "Abies")))
  
  colors <- viridis(5)
  plot_pinaceae <- ggplot(filter(long_data, variable == "Pinaceae"), aes_string(x = x_var, y = "value")) +
    geom_bar(stat = "identity", alpha = 0.6) +
    labs(title = new_site_name, subtitle = "Pinaceae - Absolute count", y = "Reads") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  max_val <- max(filter(long_data, variable != "Pinaceae")$value, na.rm = TRUE)
  plot_species_uniform <- ggplot(filter(long_data, variable != "Pinaceae"), aes_string(x = x_var, y = "value", fill = "variable")) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.6) +
    scale_fill_manual(values = colors) +
    facet_wrap(~variable, ncol = 1, scales = "fixed") +
    scale_y_continuous(limits = c(0, max_val * 1.1)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
  
  plot_species_variable <- ggplot(filter(long_data, variable != "Pinaceae"), aes_string(x = x_var, y = "value", fill = "variable")) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.6) +
    scale_fill_manual(values = colors) +
    facet_wrap(~variable, ncol = 1, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
  
  ggsave(file.path(base_dir, "results/6_7_read_counts", site_name, paste0(site_name, "_combined_count_plot_uniform_y.png")), plot_grid(plot_pinaceae, plot_species_uniform, ncol = 1, rel_heights = c(1, 2)), width = 10, height = 8, dpi = 300)
  ggsave(file.path(base_dir, "results/6_7_read_counts", site_name, paste0(site_name, "_combined_count_plot_variable_y.png")), plot_grid(plot_pinaceae, plot_species_variable, ncol = 1, rel_heights = c(1, 2)), width = 10, height = 8, dpi = 300)
  ggsave(file.path(base_dir, "results/6_7_read_counts", site_name, paste0(site_name, "_combined_count_plot_both_y.png")), plot_grid(plot_pinaceae, plot_species_uniform, plot_species_variable, ncol = 1, rel_heights = c(1, 2, 2)), width = 10, height = 12, dpi = 300)
  ggsave(file.path(base_dir, "results/6_7_read_counts", site_name, paste0(site_name, "_combined_count_plot_both_y.pdf")), plot_grid(plot_pinaceae, plot_species_uniform, plot_species_variable, ncol = 1, rel_heights = c(1, 2, 2)), width = 10, height = 12, dpi = 300, device = cairo_pdf)
}

# Run the process for each site
for (site in sites) process_site(site, site_name_mapping)
