#=============================================================================
# Script: 7.2_REARRANGE_DEPTH_COVERAGE_RESULTS.R
# Description: Cleans and reshapes depth and breadth of coverage data
#              from mapping results into a more structured format.
#=============================================================================

# Load necessary libraries
library(tidyr)
library(dplyr)

# Define input and output paths (update these according to your environment)
input_file <- "path/to/region_depth_coverage_results.csv"
output_file <- "path/to/rearranged_depth_coverage_results.csv"

# Read the CSV file
data <- read.csv(input_file)

# Clean 'Breadth of Coverage' column by removing '%' and converting to numeric
data$Breadth.of.Coverage <- as.numeric(gsub("%", "", data$Breadth.of.Coverage))

# Replace region ranges with species names
data$Region <- recode(data$Region,
                      "concatenate_reference_larix_gmelinii_pinus_pumila_abies_sibirica_picea_obovata:1-122588" = "Larix",
                      "concatenate_reference_larix_gmelinii_pinus_pumila_abies_sibirica_picea_obovata:122589-239988" = "Pinus",
                      "concatenate_reference_larix_gmelinii_pinus_pumila_abies_sibirica_picea_obovata:239989-361212" = "Abies",
                      "concatenate_reference_larix_gmelinii_pinus_pumila_abies_sibirica_picea_obovata:361213-485324" = "Picea"
)

# Pivot data so that each species has its own depth and breadth columns
rearranged_data <- data %>%
  pivot_wider(
    names_from = Region,
    values_from = c(Depth, Breadth.of.Coverage),
    names_sep = "_"
  )

# Save the rearranged data
write.csv(rearranged_data, output_file, row.names = FALSE)

cat("Rearranged results saved to", output_file, "\n")
