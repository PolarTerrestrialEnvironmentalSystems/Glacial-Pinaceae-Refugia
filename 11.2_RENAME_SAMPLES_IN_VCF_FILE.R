#===============================================================================
# Script: 11.2_RENAME_SAMPLES_IN_VCF_FILE.R
#
# Description:
# Renames sample IDs in VCF headers for multiple Pinaceae genera (Abies, Larix, 
# Picea, Pinus) using a CSV file containing sample name mappings (old → new).
#===============================================================================

library(data.table)
library(readr)

#--- Load CSV with sample name mappings ---#
# ADD the path to the CSV file containing columns: old_name, new_name
csv_file <- "ADD_PATH_FOR_file_containing_old_name_and_new_name.csv"

sample_mapping <- read_csv2(csv_file, col_types = cols(
  old_name = col_character(),
  new_name = col_character()
), show_col_types = FALSE)

# Create a named vector for mapping old sample names to new ones
mapping_dict <- setNames(sample_mapping$new_name, sample_mapping$old_name)

#--- Define VCF input/output files for each genus ---#
vcf_info <- list(
  abies = list(
    vcf_in  = "ADD_PATH_FOR_vcf_input_file_for_abies.vcf",
    vcf_out = "ADD_PATH_FOR_vcf_output_file_for_abies_renamed.vcf"
  ),
  larix = list(
    vcf_in  = "ADD_PATH_FOR_vcf_input_file_for_larix.vcf",
    vcf_out = "ADD_PATH_FOR_vcf_output_file_for_larix_renamed.vcf"
  ),
  picea = list(
    vcf_in  = "ADD_PATH_FOR_vcf_input_file_for_picea.vcf",
    vcf_out = "ADD_PATH_FOR_vcf_output_file_for_picea_renamed.vcf"
  ),
  pinus = list(
    vcf_in  = "ADD_PATH_FOR_vcf_input_file_for_pinus.vcf",
    vcf_out = "ADD_PATH_FOR_vcf_output_file_for_pinus_renamed.vcf"
  )
)

#--- Loop through each genus and rename VCF header sample IDs ---#
for (genus in names(vcf_info)) {
  vcf_file   <- vcf_info[[genus]]$vcf_in
  output_vcf <- vcf_info[[genus]]$vcf_out
  
  cat(paste0("Processing genus: ", genus, "\n"))
  
  # Read the VCF file line by line
  vcf_lines <- readLines(vcf_file)
  
  # Update the header line with new sample names
  for (i in seq_along(vcf_lines)) {
    if (startsWith(vcf_lines[i], "#CHROM")) {
      header_fields <- unlist(strsplit(vcf_lines[i], "\t"))
      header_fields[10:length(header_fields)] <- mapping_dict[header_fields[10:length(header_fields)]]
      vcf_lines[i] <- paste(header_fields, collapse = "\t")
      break
    }
  }
  
  # Write the updated VCF to output
  writeLines(vcf_lines, con = output_vcf)
}
