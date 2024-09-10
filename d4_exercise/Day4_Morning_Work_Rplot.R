#Load pre-installed libraries
library(tidyverse)
library(ggthemes)

#Open the expression data and save it as an element
Expression_Data <- read_tsv(file = "~/qbb2024-answers/d4_exercise/dicts_expr.tsv")

#Mutate the data to merge the tissue and GeneIDs' together and create a new column and another for the log of the expression values
Expression_Data <- Expression_Data %>% 
  mutate(Tissue_Gene=paste(Tissue, " ", GeneID)) %>%
  mutate(Log2_Expr = log2(Expr + 1))

#Create a violin plot to observe the data set
ggplot(data = Expression_Data, mapping = aes(x = Tissue_Gene, y = Log2_Expr)) +
  geom_violin() +
  coord_flip() +
  xlab("Tissue and Gene") +
  ylab("Log-transformed Expression Level") +
  theme_classic()

#Q# Given the tissue specificity and high expression level of these genes, are you surprised by the results?
#A# Not necessarily surprised. As some of the tissue types that require several complex process that would involve several expressing genes show high variability. While some of the tissue types that are required to be more tightly controlled and regulated show more minimal variability.

#Q# What tissue-specific differences do you see in expression variability? 
#A# Severe of the genes from the Pancrease show tight regulation of specific genes (low variability), while others may show higher variability, like those in the stomach and small intestine.
#Q# Speculate on why certain tissues show low variability while others show much higher expression variability.
#A# This variability can possibly reflect the diverse cell types or different levels of gene regulation present in those tissue types.
