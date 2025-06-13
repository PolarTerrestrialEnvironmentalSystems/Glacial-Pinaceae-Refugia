#===============================================================================
# Script: 11_SNP-BASED_SPECIES_IDENTIFICATION.R
#
# Description:
# Generalized R script for SNP-based species identification via visual comparison
# of SNPs ascertained from chloroplast reference genome alignments (Script 9)
# and those called from sediment DNA samples (Script 11).
#
# This version performs SNP-based assignment of reads to species — specifically 
# Picea obovata and Picea abies — as shown in Figure 3A and detailed in 
# Appendix S2: Table S3D (Meucci et al., submitted).
#
# Notes:
# - Appendix S2: Table S3C (Meucci et al., submitted) lists the distant chloroplast references used as references
# - Appendix S2: Table S3D (Meucci et al., submitted) lists the region/genus-specific references used for mapping and SNP call
# - Adjust the following per genus:
#     • Reference VCF file (Script 9 output)
#     • Sample VCF file (Script 11 output)
#     • `color_scheme` and reference species labels
#
#===============================================================================

#===============================================================================
# Load required libraries and define working directory
#===============================================================================
library(tidyverse)

# Set working directory (adjust per genus or analysis batch)
setwd("PATH/TO/YOUR/SAMPLES_DIRECTORY")

theme_set(theme_bw())


#===============================================================================
# Load and process SNP data from reference alignments (Script 9 output)
#===============================================================================

# Read VCF file with ascertained SNPs from reference genomes
refsbcf <- read_delim("reference_snps.vcf", delim = "\t", col_names = TRUE, skip = 61)
names(refsbcf)

# Pivot reference SNPs to long format
refsbcf <- pivot_longer(data = refsbcf,
                        cols = contains("Picea"),
                        names_to = "SAMPLES",
                        values_to ="GT:DP:AD:RO:QR:AO:QA:GL")

# Split VCF INFO fields
refsbcf <- separate(data = refsbcf,
                    col = "GT:DP:AD:RO:QR:AO:QA:GL",
                    into = c("GT","DP","AD","RO","QR","AO","QA","GL"),
                    sep = ":")

# Keep relevant columns
refsbcf <- refsbcf %>% select(c("POS", "REF", "ALT", "SAMPLES", "AD"))

# Separate allele depth (AD) into up to 7 columns
refsbcf <- separate(data = refsbcf,
                    col = "AD",
                    into = c("ref_count_1", "alt_count_1", "alt_count_2", "alt_count_3",
                             "alt_count_4", "alt_count_5", "alt_count_6", "alt_count_7"),
                    sep = ",",
                    extra = "merge")

# Separate ALT bases into up to 7 columns
refsbcf <- separate(data = refsbcf,
                    col = "ALT",
                    into = c("alt_base_1","alt_base_2","alt_base_3","alt_base_4",
                             "alt_base_5","alt_base_6","alt_base_7"),
                    sep = ",",
                    extra = "merge")

# Convert counts to numeric
numcols <- c("ref_count_1", "alt_count_1", "alt_count_2", "alt_count_3",
             "alt_count_4", "alt_count_5", "alt_count_6", "alt_count_7")
refsbcf[numcols] <- sapply(refsbcf[numcols], as.numeric)

# Rename REF column for consistency
refsbcf <- rename(refsbcf, ref_base_1=REF)

# Pivot base and count columns to long format
refsbcf <- pivot_longer(data = refsbcf,
                        cols = contains("_"),
                        names_to = c("genotype", ".value", "number"),
                        names_pattern = "(.+)_(.+)_(.+)")

# Filter out NA and zero counts
refsbcf <- refsbcf %>%
  filter(!is.na(count)) %>%
  filter(count != 0)

# Separate metadata from reference sample IDs
refsbcf <- separate(data = refsbcf,
                    col = SAMPLES,
                    into = c("id", "genus","species"),
                    sep = "_",
                    extra = "drop")

# Generate wide-format reference SNP table
refsbcfwider <- refsbcf %>%
  dplyr::select(-c(id, number)) %>%
  group_by(POS, base, species, genotype) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  group_by(POS, base, genotype) %>%
  summarise(species_counts = paste(n, collapse = ","),
            species_names = paste(species, collapse = ","),
            species_bases = paste(base, collapse = ","))


#===============================================================================
# Load and process SNP data from samples (Script 10 output)
#===============================================================================

# Read VCF of sample SNPs
al <- read_delim("sample_snps.vcf", delim = "\t", col_names = TRUE, skip = 64)
names(al)

# Pivot sample SNPs to long format
al <- pivot_longer(data = al,
                   cols = contains("_"),
                   names_to = "SAMPLES",
                   values_to = "GT:DP:AD:RO:QR:AO:QA:GL")

# Split VCF INFO fields
al <- separate(data = al,
               col = "GT:DP:AD:RO:QR:AO:QA:GL",
               into = c("GT","DP","AD","RO","QR","AO","QA","GL"),
               sep = ":")

# Select relevant columns
al <- al %>%
  dplyr::select(c("POS", "REF", "ALT", "SAMPLES", "AD", "DP"))

# Separate AD and ALT fields into up to 10 columns
al <- separate(data = al,
               col = "AD",
               into = c("ref_count_1", "alt_count_1", "alt_count_2", "alt_count_3",
                        "alt_count_4", "alt_count_5", "alt_count_6", "alt_count_7",
                        "alt_count_8", "alt_count_9", "alt_count_10"),
               sep = ",",
               extra = "merge")

numcols <- c("ref_count_1", "alt_count_1", "alt_count_2", "alt_count_3",
             "alt_count_4", "alt_count_5", "alt_count_6", "alt_count_7",
             "alt_count_8", "alt_count_9", "alt_count_10")
al[numcols] <- sapply(al[numcols], as.numeric)

al <- separate(data = al,
               col = "ALT",
               into = c("alt_base_1", "alt_base_2", "alt_base_3", "alt_base_4",
                        "alt_base_5", "alt_base_6", "alt_base_7", "alt_base_8",
                        "alt_base_9", "alt_base_10"),
               sep = ",",
               extra = "warn")

al <- rename(al, ref_base_1 = REF)

# Pivot base and count columns to long format
al <- pivot_longer(data = al,
                   cols = contains("_"),
                   names_to = c("genotype", ".value", "number"),
                   names_pattern = "(.+)_(.+)_(.+)")

#===============================================================================
# Merge sample and reference SNP tables by position, genotype, and base
#===============================================================================
alref <- full_join(al, refsbcfwider, by = c("POS", "genotype", "base"),
                   suffix = c("sample", "ref")) %>%
  filter(!is.na(count)) %>%
  filter(count != 0)

# Assign categories for missing or shared reference SNPs
alref <- alref %>%
  mutate(species_names = case_when(
    genotype == "ref" & is.na(species_names) ~ "all",
    genotype == "alt" & is.na(species_names) ~ "without ref",
    grepl(",", species_names) ~ "all",
    !grepl(",", species_names) ~ species_names
  ))

#===============================================================================
# Extract lake/age metadata from sample names
#===============================================================================
alref <- separate(data = alref,
                  col = "SAMPLES",
                  into = c("lake", "lib_id", "age"),
                  remove = FALSE,
                  extra = "drop", sep = "_")

alref <- within(alref, lake_age <- paste(lake, age, sep = '_'))

alref$age <- as.numeric(alref$age)
alref <- arrange(alref, lake, age)
alref$lake_age <- factor(alref$lake_age, levels = unique(alref$lake_age),
                         labels = unique(alref$lake_age))

# Convert POS to factor for plotting
alref$POS <- as.factor(alref$POS)

#===============================================================================
# Filter informative SNP positions (species-specific variants only)
#===============================================================================
specpos <- alref %>%
  ungroup() %>%
  filter(!(species_names %in% c("without ref", "all"))) %>%
  select(POS) %>%
  distinct() %>%
  pull(POS)

# Filter sample data by species-informative positions
alrefspecpos <- alref %>%
  filter(POS %in% specpos) %>%
  group_by(SAMPLES, POS) %>%
  mutate(percentage_specpos = count / sum(count))

length(unique(alrefspecpos$POS))

unique(alrefspecpos$species_names)


#===============================================================================
# Load required libraries for plotting
#===============================================================================
library(dplyr)
library(ggplot2)
library(ggh4x)
library(cowplot)
library(scales)  # For pretty_breaks


#===============================================================================
# Define color scheme for plotting
#===============================================================================
turbo_colors <- c(
  "#5B3A8E",  # Darker Slate Blue
  "#4178BE",  # Richer Steel Blue
  "#2D7F5E",  # Deeper Medium Sea Green
  "#B8860B",  # Darker Goldenrod
  "#FFD700",  # Golden Yellow (unchanged)
  "#8B4513",  # Saddle Brown
  "#A04000",  # Burnt Orange (unchanged)
  "#B22222",  # Dark Red (unchanged)
  "#006400"   # Dark Green (unchanged)
)

# Choose 2 different colors for the two species
color_scheme <- c(
  turbo_colors[5],  # Golden Yellow for P. obovata
  turbo_colors[8],  # Dark Red for P. abies
  "#000000"         # Black for Unknown
)
names(color_scheme) <- c("P. obovata", "P. abies", "Unknown")
color_scheme


#===============================================================================
# Format species and lake metadata for plotting
#===============================================================================
alrefspecpos <- alrefspecpos %>%
  mutate("Variation" = case_when(
    species_names == "obovata" ~ "P. obovata",
    species_names == "abies" ~ "P. abies",
    species_names == "without ref" ~ "Unknown"
  ))

speciesgroups <- c("P. obovata", "P. abies", "Unknown")
alrefspecpos$Variation <- factor(alrefspecpos$Variation, levels = speciesgroups)

# Define lake names for plot facets
alrefspecpos <- alrefspecpos %>%
  mutate("lakename" = case_when(
    lake == "hidden" ~ "Hidden",
    lake == "ulu" ~ "Ulu (Siberia, RU)",
    lake == "batagay" ~ "Batagay (Siberia, RU)",
    lake == "btoko" ~ "Bolshoe Toko (Siberia, RU)",
    lake == "emanda" ~ "Emanda (Siberia, RU)",
    lake == "kisi" ~ "Kisi (Siberia, RU)",
    lake == "lama" ~ "Lama (Siberia, RU)",
    lake == "khamra" ~ "Khamra (Siberia, RU)",
    lake == "ximen" ~ "Ximen (Tibet)",
    lake == "naleng" ~ "Naleng (Tibet)",
    lake == "imandra" ~ "Imandra (Kola, RU)",
    lake == "SF" ~ "SF (South Finland)",
    lake == "SR" ~ "SR (South Russia)",
    lake == "GroßGlienicker" ~ "Groß Glienicker (D)",
    lake == "CER" ~ "Cerne Jezero (CZ)",
    lake == "JOM" ~ "Otto Mörtsch Cave (CZ)",
    lake == "Sulsseewli" ~ "Sulsseewli (CH)",
    lake == "montenegro" ~ "Zminje Jezero (ME)",
    TRUE ~ lake
  ))

# Define plotting order for lakes
alrefspecpos$lakename <- factor(alrefspecpos$lakename, levels = c(
  "Zminje Jezero (ME)", "Sulsseewli (CH)", "Otto Mörtsch Cave (CZ)", "Cerne Jezero (CZ)",
  "Groß Glienicker (D)", "SF (South Finland)", "SR (South Russia)", "Imandra (Kola, RU)",
  "Lama (Siberia, RU)", "Khamra (Siberia, RU)", "Bolshoe Toko (Siberia, RU)",
  "Batagay (Siberia, RU)", "Emanda (Siberia, RU)", "Ulu (Siberia, RU)",
  "Kisi (Siberia, RU)", "Ximen (Tibet)", "Naleng (Tibet)", "Hidden"
))

# Add group index for plotting SNP positions
alrefspecpos <- alrefspecpos %>% group_by(POS) %>% mutate(POS2 = cur_group_id())


#===============================================================================
# SNP assignment barplot — individual positions
#===============================================================================
p1 <- ggplot(data = subset(alrefspecpos), aes(x = POS2, y = percentage_specpos, fill = Variation)) +
  scale_fill_manual(values = color_scheme, labels = c(
    expression(italic("P. obovata")),
    expression(italic("P. abies")),
    "Unknown"
  ), name = "Assignment to references") +
  scale_x_continuous(expand = c(0, 0), breaks = seq(25, 150, 25)) +
  theme(legend.text.align = 0) +
  geom_bar(stat = "identity") +
  facet_nested(lakename + round(age, -2) ~ .) +
  theme(panel.spacing = unit(0, "lines")) +
  theme(strip.text.y = element_text(angle = 0)) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme(legend.position = "bottom") +
  xlab("Variable positions on chloroplast genome") +
  ylab("Read counts assigned to reference per position [%]") +
  theme(strip.text.x = element_text(size = 8))

p1


#===============================================================================
# SNP assignment barplot — total counts per species per lake/age
#===============================================================================
alrefspecpos <- alrefspecpos %>%
  ungroup() %>%
  group_by(SAMPLES) %>%
  mutate(percentage_species_names = count / sum(count))

alrefspecpossum <- alrefspecpos %>%
  group_by(lake_age, species_names, age, Variation, lakename) %>%
  summarise(countsum = sum(count))

p2 <- ggplot(data = alrefspecpossum, aes(x = as.factor(round(age, -2)), y = countsum, fill = Variation)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(limits = rev, position = "top", name = "") +
  scale_fill_manual(values = color_scheme, labels = c(
    expression(italic("P. obovata")),
    expression(italic("P. abies")),
    "Unknown"
  ), name = "Assignment to references") +
  facet_grid(rows = vars(lakename), scales = "free", space = "free", switch = "y") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  ylab("Read counts assigned to reference") +
  theme(panel.spacing = unit(0.2, "lines")) +
  theme(strip.text.y.left = element_text(angle = 0), axis.title.y.left = element_blank()) +
  theme(legend.position = "bottom") +
  coord_flip()

p2


#===============================================================================
# SNP assignment barplot — filtered absolute counts
#===============================================================================
alrefspecpossum_filtered <- alrefspecpossum %>%
  filter(!is.na(Variation), 
         !grepl("btoko|naleng|ximen|hidden|kisi", lake_age))

p2_filtered <- ggplot(data = alrefspecpossum_filtered, aes(x = as.factor(round(age, -2)), y = countsum, fill = Variation)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_y_continuous(breaks = pretty_breaks(n = 20), expand = c(0, 0)) +
  scale_x_discrete(limits = rev, position = "top", name = "") +
  scale_fill_manual(values = color_scheme, labels = c(
    expression(italic("P. obovata")),
    expression(italic("P. abies")),
    "Unknown"
  ), name = "Assignment to references") +
  facet_grid(rows = vars(lakename), scales = "free", space = "free", switch = "y") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme(axis.text.y = element_text(size = 7)) +
  ylab("Absolute read counts assigned to reference") +
  theme(panel.spacing = unit(0.2, "lines")) +
  theme(strip.text.y.left = element_text(angle = 0), axis.title.y.left = element_blank()) +
  theme(legend.position = "bottom") +
  coord_flip()

p2_filtered


#===============================================================================
# SNP assignment barplot — relative proportions
#===============================================================================
alrefspecpossum_filtered_rel_count <- alrefspecpossum_filtered %>%
  group_by(lake_age) %>%
  summarise(total_count = sum(countsum)) %>%
  left_join(alrefspecpossum_filtered, by = "lake_age") %>%
  mutate(rel_count = countsum / total_count * 100) %>%
  select(-total_count)

p2_filtered_rel_count <- ggplot(data = alrefspecpossum_filtered_rel_count, aes(x = as.factor(lake_age), y = rel_count, fill = Variation)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(limits = rev, position = "top", name = "") +
  scale_fill_manual(values = color_scheme, labels = c(
    expression(italic("P. obovata")),
    expression(italic("P. abies")),
    "Unknown"
  ), name = "Assignment to references") +
  facet_grid(rows = vars(lakename), scales = "free", space = "free", switch = "y") +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 6)) +
  ylab("Relative read counts assigned to reference (%)") +
  theme(panel.spacing = unit(0.2, "lines")) +
  theme(strip.text.y.left = element_text(angle = 0), axis.title.y.left = element_blank()) +
  theme(legend.position = "bottom") +
  coord_flip()

p2_filtered_rel_count


#===============================================================================
# Combine absolute and relative plots into one figure with shared legend
#===============================================================================
combined_plot <- plot_grid(
  p2_filtered + theme(legend.position = "none"), 
  p2_filtered_rel_count + theme(legend.position = "none", axis.text.y = element_blank()), 
  labels = c("A", "B"),
  ncol = 2, 
  align = "h",
  axis = "tb",
  rel_widths = c(1, 1)
)

legend <- get_legend(p2_filtered + theme(legend.position = "bottom"))

final_plot <- plot_grid(
  combined_plot,
  legend,
  ncol = 1,
  rel_heights = c(1, 0.2)
) +
  draw_label("Sample age", x = 0.5, y = 0.03, hjust = 0.5, vjust = 1, angle = 0, size = 12)

# Save combined figure
ggsave("eurasia/picea_snp_eurasia_ref_pungens_combined.png", plot = final_plot, width = 16, height = 8, dpi = 300)
ggsave("eurasia/picea_snp_eurasia_ref_pungens_combined.pdf", plot = final_plot, width = 16, height = 8, dpi = 300)

final_plot
