library(ggplot2)
library(tidyr)
library(dplyr)

# Import the Data
data <- read.table("~/qbb2024-answers/Quant_Bio_Lab/Homework/Week_10/week_10_signals.txt", header = TRUE, sep = "\t")
view(data) 

# Convert Formatting to have Gene, Measure, and Value for easier plotting
data_long <- data %>%
  gather(key = "Measure", value = "Value", Mean_NascentRNA, Mean_PCNA, Ratio)
view(data_long)


# Plot nascent RNA
nascent_RNA <- ggplot(subset(data_long, Measure == "Mean_NascentRNA"),
                     aes(x = Gene, y = Value, fill = Gene)) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  labs(title = "Gene-Specific Nascent RNA Signal",
       x = "Gene",
       y = "Nascent RNA Signal") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45))
plot(nascent_RNA)
ggsave("~/qbb2024-answers/Quant_Bio_Lab/Homework/Week_10/nascent_RNA_Plot.png", plot = nascent_RNA, bg = "white")

# Plot PCNA
pcna_plot <- ggplot(subset(data_long, Measure == "Mean_PCNA"), aes(x = Gene, y = Value, fill = Gene)) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  labs(title = "PCNA Expression Across Genes",
       x = "Gene",
       y = "PCNA Signal") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45))
plot(pcna_plot)
ggsave("~/qbb2024-answers/Quant_Bio_Lab/Homework/Week_10/pcna_plot.png", plot = pcna_plot, bg = "white")


# Plot log2 ratio 
log_2_ratio_plot <- ggplot(subset(data_long, Measure == "Ratio"), aes(x = Gene, y = Value, fill = Gene)) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  labs(title = "Log2 Ratio of Nascent RNA to PCNA Across Genes",
       x = "Gene",
       y = "Log2 Ratio (Nascent RNA / PCNA)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45))
plot(log_2_ratio_plot)
ggsave("~/qbb2024-answers/Quant_Bio_Lab/Homework/Week_10/log_2_ratio_plot.png", plot = log_2_ratio_plot, bg = "white")
