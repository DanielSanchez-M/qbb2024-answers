# Load necessary library
library(ggplot2)
library(tidyverse)

# Import the data as snp_data_frame
snp_data_frame <- read.table("~/qbb2024-answers/Quant_Bio_Lab/Homework/Week_2/snp_counts.txt", header=TRUE, sep="\t")

#view(snp_data_frame)

# Convert MAF to a factor to ensure correct plotting order
snp_data_frame$MAF <- as.factor(snp_data_frame$MAF)

# Create the ggplot
snp_enrichment_plot <- ggplot(snp_data_frame, aes(x=MAF, y=Enrichment, color=Feature, group=Feature)) +
  geom_line() +
  geom_point() +
  labs(title="SNP Enrichment Across MAF Levels",
       x="Minor Allele Frequency (MAF)",
       y="log2(Enrichment)") +
  theme_minimal() +
  scale_y_continuous(trans = "log2")

print(snp_enrichment_plot)

# save the plot as a pdf
ggsave("final_snp_enrichment.pdf", plot=snp_enrichment_plot)
