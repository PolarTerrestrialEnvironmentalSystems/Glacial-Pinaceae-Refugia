# Load necessary libraries
library(ggplot2)
library(reshape2)
library(viridis)

# File Paths (adapt these for each genus)
admixture_results <- "ADD_PATH_TO/GENUS/GENUS_admixture_results.K.Q"
sample_names_file <- "ADD_PATH_TO/GENUS/GENUS_sample_names_for_admixture_plot.txt"
group_names_file <- "ADD_PATH_TO/GENUS/group_names_GENUS.csv"
output_plot_png <- "ADD_PATH_TO/GENUS/GENUS_merged_admixture_plot_KQ.png"
output_plot_pdf <- "ADD_PATH_TO/GENUS/GENUS_merged_admixture_plot_KQ.pdf"

# Read data
tbl <- read.table(admixture_results, header = FALSE, stringsAsFactors = FALSE)
sample_names <- read.table(sample_names_file, header = FALSE, stringsAsFactors = FALSE)
group_names <- read.csv(group_names_file, header = TRUE, stringsAsFactors = FALSE, sep = "\t")

# Ensure the encoding of the group names is UTF-8
group_names$new_name <- enc2utf8(group_names$new_name)

# Extract the groups and order numbers from sample names
sample_names$Group <- sapply(strsplit(sample_names$V1, "_"), `[`, 1)
sample_names$Order <- sapply(strsplit(sample_names$V1, "_"), `[`, 3)

# Replace NA orders with "B"
sample_names$Order[is.na(sample_names$Order)] <- "B"

# Map old group names to new names
group_map <- setNames(group_names$new_name, group_names$old_name)
sample_names$Group <- ifelse(sample_names$Group %in% names(group_map), group_map[sample_names$Group], sample_names$Group)

# Create the modified sample names as "NewGroup_Order"
sample_names$ModifiedName <- paste(sample_names$Group, sample_names$Order, sep = "_")

# Combine the data
data_combined <- cbind(sample_names, tbl)

# Rename the columns for clarity (adapt number of ancestries if needed)
colnames(data_combined) <- c("Sample", "Group", "Order", "ModifiedName",
                             paste0("Ancestry", seq_len(ncol(tbl))))

# Add blank rows as separators between groups
unique_groups <- unique(data_combined$Group)
separator_rows <- data.frame(Sample = NA, Group = unique_groups, Order = NA, 
                             ModifiedName = paste(unique_groups, "separator", sep = "_"),
                             matrix(0, nrow = length(unique_groups), ncol = ncol(tbl)),
                             stringsAsFactors = FALSE)
colnames(separator_rows)[5:ncol(separator_rows)] <- paste0("Ancestry", seq_len(ncol(tbl)))

# Combine the original data with separator rows
data_combined <- rbind(data_combined, separator_rows)

# Reshape the data into long format for plotting
data_melted <- melt(data_combined, id.vars = c("Sample", "Group", "Order", "ModifiedName"))

# Order the groups based on the CSV file
data_melted$Group <- factor(as.character(data_melted$Group), levels = group_names$new_name)

# Convert Order to numeric where possible for correct sorting
data_melted$Order <- suppressWarnings(as.numeric(as.character(data_melted$Order)))

# Order samples within each group
data_melted <- data_melted[order(data_melted$Group, data_melted$Order, data_melted$ModifiedName), ]

# Set x-axis order
data_melted$ModifiedName <- factor(data_melted$ModifiedName, levels = unique(data_melted$ModifiedName))

# Ensure UTF-8 encoding for labels
levels(data_melted$ModifiedName) <- enc2utf8(levels(data_melted$ModifiedName))

# Create stacked bar plot
p <- ggplot(data = data_melted, aes(x = ModifiedName, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "stack", color = "white", size = 0.5) +
  labs(title = "GENUS", x = "Samples", y = "Proportion of Ancestry") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
    legend.title = element_blank(),
    legend.position = "bottom"
  ) +
  scale_fill_viridis_d(option = "viridis")

# Remove x-axis text for separators
p <- p + scale_x_discrete(labels = function(x) ifelse(grepl("separator", x), "", x))

# Display the plot
print(p)

# Save the plot
ggsave(output_plot_png, p, width = 12, height = 8, dpi = 300)
ggsave(output_plot_pdf, p, width = 12, height = 8)
