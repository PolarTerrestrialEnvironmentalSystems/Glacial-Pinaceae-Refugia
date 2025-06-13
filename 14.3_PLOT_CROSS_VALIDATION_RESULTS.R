# Load required libraries
library(ggplot2)
library(patchwork)  # For combining plots

### PICEA
picea_data <- read.csv("ADD_PATH_TO/picea/cross_validation_results.csv")
picea_data$K <- factor(picea_data$K, levels = unique(picea_data$K))
picea_plot <- ggplot(picea_data, aes(x = K, y = Cross.validation)) +
  geom_point() +
  geom_line(aes(group = 1), color = "blue", linetype = "solid") +
  labs(title = "Picea",
       x = "Number of Clusters (K)",
       y = "Cross-Validation error") +
  theme_minimal() +
  theme(panel.grid.major = element_line(color = "gray", linetype = "dashed"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))

### PINUS
pinus_data <- read.csv("ADD_PATH_TO/pinus/cross_validation_results.csv")
pinus_data$K <- factor(pinus_data$K, levels = unique(pinus_data$K))
pinus_plot <- ggplot(pinus_data, aes(x = K, y = Cross.validation)) +
  geom_point() +
  geom_line(aes(group = 1), color = "blue", linetype = "solid") +
  labs(title = "Pinus",
       x = "Number of Clusters (K)",
       y = "Cross-Validation error") +
  theme_minimal() +
  theme(panel.grid.major = element_line(color = "gray", linetype = "dashed"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))

### LARIX
larix_data <- read.csv("ADD_PATH_TO/larix/cross_validation_results.csv")
larix_data$K <- factor(larix_data$K, levels = unique(larix_data$K))
larix_plot <- ggplot(larix_data, aes(x = K, y = Cross.validation)) +
  geom_point() +
  geom_line(aes(group = 1), color = "blue", linetype = "solid") +
  labs(title = "Larix",
       x = "Number of Clusters (K)",
       y = "Cross-Validation error") +
  theme_minimal() +
  theme(panel.grid.major = element_line(color = "gray", linetype = "dashed"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))

### ABIES
abies_data <- read.csv("ADD_PATH_TO/abies/cross_validation_results.csv")
abies_data$K <- factor(abies_data$K, levels = unique(abies_data$K))
abies_plot <- ggplot(abies_data, aes(x = K, y = Cross.validation)) +
  geom_point() +
  geom_line(aes(group = 1), color = "blue", linetype = "solid") +
  labs(title = "Abies",
       x = "Number of Clusters (K)",
       y = "Cross-Validation error") +
  theme_minimal() +
  theme(panel.grid.major = element_line(color = "gray", linetype = "dashed"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))

# Combine plots in a 2x2 grid
combined_plot <- (picea_plot | pinus_plot) / (larix_plot | abies_plot)

# Save combined plot
ggsave("ADD_PATH_TO/combined_cross_validation_plot.png", combined_plot, width = 16, height = 12, units = "in", dpi = 300)
ggsave("ADD_PATH_TO/combined_cross_validation_plot.pdf", combined_plot, width = 16, height = 12)

# Display the plot
print(combined_plot)
