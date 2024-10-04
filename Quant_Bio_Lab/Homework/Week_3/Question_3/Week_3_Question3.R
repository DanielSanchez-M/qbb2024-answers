# Load necessary library
library(ggplot2)
library(tidyverse)

#Import the Allele Frequency and Read Depth Distriubtion Text Files and save them as Data Frames
Allele_Frequencies_DF <- read_tsv("~/qbb2024-answers/Quant_Bio_Lab/Homework/Week_3/Question_3/AF.txt")
Read_Depth_DF <- read_tsv("~/qbb2024-answers/Quant_Bio_Lab/Homework/Week_3/Question_3/DP.txt")

# Allele Frequency Spectrum Histogram
# Create a histogram using ggplot with bins set to 11 to avoid binning artifacts
Allele_Frequency_Spectrum <- ggplot(data = Allele_Frequencies_DF, aes(x = `Allele Frequencies`)) +
  geom_histogram(bins = 11, fill = "orange") +
  labs(x = "Allele Frequency", y = "Number of Variations", title = "Allele Frequency Spectrum")
print(Allele_Frequency_Spectrum)

# save the Allele Frequency Histogram as a pdf
ggsave("~/qbb2024-answers/Quant_Bio_Lab/Homework/Week_3/Question_3/Allele_Frequency_Spectrum.pdf", plot=Allele_Frequency_Spectrum)


# Read Depth Distribution Histogram
# Create a histogram using ggplot with bins set to 21 and x-limitor between 0,20 to make it legible 
Read_Depth_Distribution <- ggplot(data = Read_Depth_DF, aes(x = `Read Depth of Each Variant`)) +
  geom_histogram(bins = 21, fill = "orange") +
  xlim(0, 20) +
  labs(x = "Coverage", y = "Number of Variations", title = "Number of Variations Covered")
print(Read_Depth_Distribution)
  
# save the  Read Depth Distribution Histogram as a pdf
ggsave("~/qbb2024-answers/Quant_Bio_Lab/Homework/Week_3/Question_3/Read_Depth_Distribution.pdf", plot=Read_Depth_Distribution)
