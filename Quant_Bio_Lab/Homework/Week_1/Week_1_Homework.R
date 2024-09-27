library(tidyverse)
library(ggthemes)

# 3x Coverage
# Import Data
genome_coverage <- read.delim("~/qbb2024-answers/Quant_Bio_Lab/Homework/Week_1/genome_coverage_3x.txt")

# Rename Data Frame Column
colnames(genome_coverage) <- c("Coverage")

# Calculate Frequencies
coverage_freq <- genome_coverage %>%
  group_by(Coverage) %>%
  summarize(frequency = n())

# Extract Coverage values for Poisson and Normal calculations
coverage_values <- coverage_freq$Coverage

# Calculate Poisson PMF
poisson_pmf <- dpois(coverage_values, lambda = 3)

# Calculate Normal PDF
normal_pdf <- dnorm(coverage_values, mean = 3, sd = sqrt(3))

# Make a histogram plot for genome coverage with poisson and normal distribution overlayed
x3_plot <- ggplot() +
  geom_histogram(data = genome_coverage, aes(x = Coverage, fill = "Genome 3x Coverage"), binwidth = 1, color = "black", alpha = 0.5) + # Adds Genome coverage as a histogram
  geom_line(data = coverage_freq, aes(x = Coverage, y = poisson_pmf * nrow(genome_coverage), color = "Poisson Distribution"), size = 1) +  # Adds Poisson Distribution as a line plot
  geom_line(data = coverage_freq, aes(x = Coverage, y = normal_pdf * nrow(genome_coverage), color = "Normal Distribution"), size = 1) +  # Adds Normal Distribution as a line plot
  labs(title = "Genome 3x Coverage Distribution with Poisson \n and Normal Distribution",
       x = "Coverage",
       y = "Frequency") +
  scale_fill_manual(name = "Data", values = c("Genome 3x Coverage" = "blue")) +  # Adds histogram to the legend
  scale_color_manual(name = "Distributions", values = c("Poisson Distribution" = "red", "Normal Distribution" = "green")) +  # Adds the distribution lines to the legend
  theme_minimal()

ggsave(filename = "~/qbb2024-answers/Quant_Bio_Lab/Homework/Week_1/ex1_3x_cov.jpg", plot = x3_plot, width = 10, height = 6, dpi = 300)



# For 10x Coverage
# Import Data
genome_coverage <- read.delim("~/qbb2024-answers/Quant_Bio_Lab/Homework/Week_1/genome_coverage_10x.txt")

# Rename Data Frame Column
colnames(genome_coverage) <- c("Coverage")

# Calculate Frequencies
coverage_freq <- genome_coverage %>%
  group_by(Coverage) %>%
  summarize(frequency = n())

# Extract Coverage values for Poisson and Normal calculations
coverage_values <- coverage_freq$Coverage

# Calculate Poisson PMF
poisson_pmf <- dpois(coverage_values, lambda = 10)

# Calculate Normal PDF
normal_pdf <- dnorm(coverage_values, mean = 10, sd = sqrt(10))

# Make a histogram plot for genome coverage with poisson and normal distribution overlayed
x10_plot <- ggplot() +
  geom_histogram(data = genome_coverage, aes(x = Coverage, fill = "Genome 10x Coverage"), binwidth = 1, color = "black", alpha = 0.5) + # Adds Genome coverage as a histogram
  geom_line(data = coverage_freq, aes(x = Coverage, y = poisson_pmf * nrow(genome_coverage), color = "Poisson Distribution"), size = 1) +  # Adds Poisson Distribution as a line plot
  geom_line(data = coverage_freq, aes(x = Coverage, y = normal_pdf * nrow(genome_coverage), color = "Normal Distribution"), size = 1) +  # Adds Normal Distribution as a line plot
  labs(title = "Genome 10x Coverage Distribution with Poisson \n and Normal Distribution",
       x = "Coverage",
       y = "Frequency") +
  scale_fill_manual(name = "Data", values = c("Genome 10x Coverage" = "blue")) +  # Adds histogram to the legend
  scale_color_manual(name = "Distributions", values = c("Poisson Distribution" = "red", "Normal Distribution" = "green")) +  # Adds the distribution lines to the legend
  theme_minimal()

ggsave(filename = "~/qbb2024-answers/Quant_Bio_Lab/Homework/Week_1/ex1_10x_cov.jpg", plot = x10_plot, width = 10, height = 6, dpi = 300)



# For 30x Coverage
# Import Data
genome_coverage <- read.delim("~/qbb2024-answers/Quant_Bio_Lab/Homework/Week_1/genome_coverage_30x.txt")

# Rename Data Frame Column
colnames(genome_coverage) <- c("Coverage")

# Calculate Frequencies
coverage_freq <- genome_coverage %>%
  group_by(Coverage) %>%
  summarize(frequency = n())

# Extract Coverage values for Poisson and Normal calculations
coverage_values <- coverage_freq$Coverage

# Calculate Poisson PMF
poisson_pmf <- dpois(coverage_values, lambda = 30)

# Calculate Normal PDF
normal_pdf <- dnorm(coverage_values, mean = 30, sd = sqrt(30))

# Make a histogram plot for genome coverage with poisson and normal distribution overlayed
x30_plot <- ggplot() +
  geom_histogram(data = genome_coverage, aes(x = Coverage, fill = "Genome 30x Coverage"), binwidth = 1, color = "black", alpha = 0.5) + # Adds Genome coverage as a histogram
  geom_line(data = coverage_freq, aes(x = Coverage, y = poisson_pmf * nrow(genome_coverage), color = "Poisson Distribution"), size = 1) +  # Adds Poisson Distribution as a line plot
  geom_line(data = coverage_freq, aes(x = Coverage, y = normal_pdf * nrow(genome_coverage), color = "Normal Distribution"), size = 1) +  # Adds Normal Distribution as a line plot
  labs(title = "Genome 30x Coverage Distribution with Poisson \n and Normal Distribution",
       x = "Coverage",
       y = "Frequency") +
  scale_fill_manual(name = "Data", values = c("Genome 30x Coverage" = "blue")) +  # Adds histogram to the legend
  scale_color_manual(name = "Distributions", values = c("Poisson Distribution" = "red", "Normal Distribution" = "green")) +  # Adds the distribution lines to the legend
  theme_minimal()

ggsave(filename = "~/qbb2024-answers/Quant_Bio_Lab/Homework/Week_1/ex1_30x_cov.jpg", plot = x30_plot, width = 10, height = 6, dpi = 300)
