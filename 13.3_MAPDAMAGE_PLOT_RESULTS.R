# Load required libraries
library(ggplot2)
library(dplyr)
library(viridis)
library(cowplot)
library(readr)

# Define the base directory containing lake directories
base_dir <- "ADD_PATH_FOR_directory_containing_mapdamage_results"

# Load the new lake names along with the corresponding numbers
new_names_file <- file.path(base_dir, "new_names.txt")
new_names <- read_delim(new_names_file, delim = "\t", col_names = TRUE)

# Ensure lake directories are in the same order as new_names
lake_dirs <- list.dirs(base_dir, recursive = FALSE)

# Create a mapping of old folder names to new lake names
lake_name_mapping <- setNames(new_names$new_name, new_names$folder_name)

# Define the subdirectory name containing the merged data for Pinaceae
merged_folder_name <- "merged_concatenate_reference_larix_gmelinii_pinus_pumila_abies_sibirica_picea_obovata_bam_files_tax3318_output"

# Function to create custom legend labels with the corresponding numbers for each genus
create_custom_labels <- function(site_row) {
  abies_label <- paste0("Abies (", site_row$Abies, ")")
  larix_label <- paste0("Larix (", site_row$Larix, ")")
  picea_label <- paste0("Picea (", site_row$Picea, ")")
  pinus_label <- paste0("Pinus (", site_row$Pinus, ")")
  
  return(c(
    "merged_abies" = abies_label,
    "merged_larix" = larix_label,
    "merged_picea" = picea_label,
    "merged_pinus" = pinus_label
  ))
}

# Function to process and create the Pinaceae plot
process_pinaceae_plot <- function(lake_dir) {
  lake_name <- basename(lake_dir)
  new_lake_name <- lake_name_mapping[lake_name]
  
  # Get the corresponding row from the new_names data
  site_row <- new_names %>% filter(new_name == new_lake_name)
  
  # Define the path to the merged folder inside the lake directory
  merged_dir <- file.path(lake_dir, merged_folder_name)
  
  # Initialize a data frame to hold the frequency data for Pinaceae
  data_freq <- data.frame()
  
  # Path to the frequency file
  freq_file_path <- file.path(merged_dir, '5pCtoT_freq.txt')
  
  # Check for frequency data
  if (file.exists(freq_file_path)) {
    df_freq <- read.table(freq_file_path, header = TRUE, sep = "\t")
    
    # Filter for the first 25 positions
    c_to_t_freq <- df_freq %>% filter(pos <= 25) %>% select(pos, `X5pC.T`)
    
    # Check if the resulting data frame is not empty
    if (nrow(c_to_t_freq) > 0) {
      data_freq <- c_to_t_freq
    }
  }
  
  # If frequency data is empty, skip processing
  if (nrow(data_freq) == 0) {
    message("No data to process for Pinaceae in lake directory: ", lake_dir)
    return(NULL)
  }
  
  # Rename columns for plotting
  colnames(data_freq) <- c("Position", "Frequency")
  
  # Create the label for Pinaceae with the number of reads
  pinaceae_label <- paste0("Pinaceae (", site_row$Pinaceae, ")")
  
  # Create a color palette for the single Pinaceae line
  color_palette <- c("Pinaceae" = viridis(1)[1])
  
  # Generate the plot for C>T frequency (one line for Pinaceae)
  p_pinaceae <- ggplot(data_freq, aes(x = Position, y = Frequency, color = "Pinaceae")) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    scale_color_manual(values = color_palette, labels = pinaceae_label) +
    scale_x_continuous(breaks = 1:25, labels = 1:25) +
    labs(title = new_lake_name, x = "Position", y = "C>T frequency") +
    theme_bw() +
    theme(
      panel.grid.major = element_line(size = 0.5),
      panel.grid.minor = element_line(size = 0.25),
      legend.position = "bottom",
      legend.title = element_blank(),
      text = element_text(size = 18),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12)
    )
  
  return(p_pinaceae)
}

# Function to process the species-specific plot (Abies, Larix, Picea, Pinus)
process_species_plot <- function(lake_dir) {
  lake_name <- basename(lake_dir)
  new_lake_name <- lake_name_mapping[lake_name]
  
  # Get the corresponding row from the new_names data
  site_row <- new_names %>% filter(new_name == new_lake_name)
  
  # Get the list of sample directories
  sample_dirs <- list.dirs(lake_dir, recursive = FALSE)
  
  # Initialize a data frame to hold the frequency data
  data_freq <- data.frame()
  
  # Extract C>T frequencies for the first 25 positions for each sample
  for (sample_dir in sample_dirs) {
    freq_file_path <- file.path(sample_dir, '5pCtoT_freq.txt')
    
    if (file.exists(freq_file_path)) {
      df_freq <- read.table(freq_file_path, header = TRUE, sep = "\t")
      c_to_t_freq <- df_freq %>% filter(pos <= 25) %>% select(pos, `X5pC.T`)
      
      if (nrow(c_to_t_freq) > 0) {
        parts <- strsplit(basename(sample_dir), "_")[[1]]
        sample_name <- paste(parts[1], parts[2], sep = "_")
        c_to_t_freq$Sample <- sample_name
        data_freq <- rbind(data_freq, c_to_t_freq)
      }
    }
  }
  
  if (nrow(data_freq) == 0) {
    message("No data to process for species in lake directory: ", lake_dir)
    return(NULL)
  }
  
  colnames(data_freq) <- c("Position", "Frequency", "Sample")
  
  color_palette <- c(
    "merged_abies" = "#FF6347",
    "merged_larix" = "#4682B4",
    "merged_picea" = "#32CD32",
    "merged_pinus" = "#FFD700"
  )
  
  custom_labels <- create_custom_labels(site_row)
  
  valid_samples <- c("merged_abies", "merged_larix", "merged_picea", "merged_pinus")
  data_freq <- data_freq %>% filter(Sample %in% valid_samples)
  
  p_species <- ggplot(data_freq, aes(x = Position, y = Frequency, color = Sample, group = Sample)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    scale_color_manual(values = color_palette, labels = custom_labels) +
    scale_x_continuous(breaks = 1:25, labels = 1:25) +
    labs(title = NULL, x = "Position", y = "C>T frequency") +
    theme_bw() +
    theme(
      panel.grid.major = element_line(size = 0.5),
      panel.grid.minor = element_line(size = 0.25),
      legend.position = "bottom",
      legend.title = element_blank(),
      text = element_text(size = 18),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12)
    )
  
  return(p_species)
}

# Process each lake directory and generate combined plots (Pinaceae + Species)
process_lake <- function(lake_dir) {
  p_pinaceae <- process_pinaceae_plot(lake_dir)
  p_species <- process_species_plot(lake_dir)
  
  if (!is.null(p_pinaceae) & !is.null(p_species)) {
    combined_plot <- plot_grid(
      p_pinaceae, 
      p_species, 
      ncol = 1, 
      align = "v", 
      rel_heights = c(1, 1.5)
    )
    return(combined_plot)
  } else {
    return(NULL)
  }
}

# Process each lake directory and generate combined plots, store them in a list
plot_list <- lapply(lake_dirs, process_lake)
plot_list <- plot_list[!sapply(plot_list, is.null)]

total_plots <- length(plot_list)
plots_per_file <- 4
num_files <- ceiling(total_plots / plots_per_file)

for (i in 1:num_files) {
  start_index <- (i - 1) * plots_per_file + 1
  end_index <- min(i * plots_per_file, total_plots)
  plot_subset <- plot_list[start_index:end_index]
  
  combined_plot <- plot_grid(plotlist = plot_subset, ncol = 2, align = "hv")
  
  output_png_file <- file.path(base_dir, paste0("CtoT_frequency_plots_PINACEAE_species_part_", i, ".png"))
  output_pdf_file <- file.path(base_dir, paste0("CtoT_frequency_plots_PINACEAE_species_part_", i, ".pdf"))
  
  ggsave(output_png_file, plot = combined_plot, width = 18, height = 16, units = "in", dpi = 300, bg = "white", limitsize = FALSE)
  ggsave(output_pdf_file, plot = combined_plot, width = 18, height = 16, units = "in", dpi = 300)
}
