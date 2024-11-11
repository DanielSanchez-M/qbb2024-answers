# Exercise 1: Perform principal component analysis
#____________________________________________________________

# Step 1.1: Loading data and importing libraries
library(tidyverse)
library(broom)
library(DESeq2)   #Statistical background runs

#Set working directory to where files are stored
setwd("~/qbb2024-answers/Quant_Bio_Lab/Homework/Week_7/")

#load the blood counts as a dataframe
counts_df <- read_delim("Week_7_Data/gtex_whole_blood_counts_downsample.txt")
counts_df[1:5,] # Check the data

# move the gene names to row names
counts_df <- column_to_rownames(counts_df, var = "GENE_NAME")
counts_df[1:5,] # Check the data

# load the metadata as a dataframe and change the SUBJECT ID
metadata_df <- read_delim("Week_7_Data/gtex_metadata_downsample.txt")
metadata_df[1:5,] # Check the data

# move the subject IDs to row names
metadata_df <- column_to_rownames(metadata_df, var = "SUBJECT_ID")
metadata_df[1:5,] # Check the data

# Check column from counts data and row from metadata match
table(colnames(counts_df) == rownames(metadata_df))
#____________________________________________________________
# Step 1.2: Create a DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts_df,
                              colData = metadata_df,
                              design = ~ SEX + DTHHRDY + AGE)
#____________________________________________________________
# Step 1.3: Normalization and PCA
# apply VST transformation
vsd <- vst(dds)

PCA_SEX <- plotPCA(vsd, intgroup = "SEX") + 
  ggtitle("PCA of Gene Expression by Sex") +
  labs(color = "Sex")

PCA_DTHHRDY <- plotPCA(vsd, intgroup = "DTHHRDY") + 
  ggtitle("PCA of Gene Expression by Death Hardiness") +
  labs(color = "Death Hardiness")

PCA_AGE <- plotPCA(vsd, intgroup = "AGE") + 
  ggtitle("PCA of Gene Expression by Age Group") +
  labs(color = "Age Group")

ggsave("PCA_SEX.pdf", plot=PCA_SEX)
ggsave("PCA_DTHHRDY.pdf", plot=PCA_DTHHRDY)
ggsave("PCA_AGE.pdf", plot=PCA_AGE)

## Question: 1.3.3: What proportion of variance in the gene expression data is explained by each of the first two principal components? 
## Which principal components appear to be associated with which subject-level variables? 
## Interpret these patterns in your own words and record your answers as a comment in your code.

###### CORRECTED ANSWER ###########
### Answer: The first two principal components have a proportion of variance of either 48% for PC1 and 7% for PC2. After looking at each of the different plots, PC1 (which explains a significant 48% of the variance) may be capturing differences related to death hardiness, as there is notable separation in the PCA plot by this variable. This suggests that a large proportion of the variation in the dataset could be related to factors associated with death hardiness or conditions prior to death. While PC2 (7% variance) seems to be more associated with age group, as there is some degree of separation in that dimension with "older age" being more central with minimal vairance and those of younger age having a larger degree of variance, indicating that age may have a more subtle but still detectable effect on gene expression.

#____________________________________________________________
# Exercise 2: Perform differential expression analysis
#____________________________________________________________

# Step 2.1: Perform a “homemade” test for differential expression between the sexes
vsd_df <- assay(vsd) %>%
  t() %>%
  as_tibble()
vsd_df <- bind_cols(metadata_df, vsd_df)

## 2.1.1: Does WASH7P show significant evidence of sex-differential expression (and if so, in which direction)? Explain your answer.
hist(vsd_df$WASH7P) # Visualize the WASH7P expression distribution

lm(data = vsd_df, formula = WASH7P ~ SEX + DTHHRDY + AGE) %>%
  summary() %>%
  tidy()
### Answer: Taking a look at the p-value of sex_male that is 0.279, this is above the 0.05 threshold indicating the association between sex and the gene's expression is not statistically significant. Although with the estimate being 0.119, idicating a slight shift towards male preference, due to the data being insignificant, it can not be conclusive that the shift is reliable.

#____________________________________________________________
## 2.1.2: Now repeat your analysis for the gene MAP7D2. Does this gene show evidence of sex-differential expression (and if so, in which direction)? Explain your answer.
hist(vsd_df$SLC25A47) # Visualize the MAP7D2 expression distribution

lm(data = vsd_df, formula = SLC25A47 ~ SEX + DTHHRDY + AGE) %>%
  summary() %>%
  tidy()
### Answer: Since the p-value for this gene is 0.0257, which is below the 0.05 p-value threshold, then the data is significant and one can suggest that this gene has sex-differential expression. Along with the estimated value of 0.518, this indicates that the sex-differential is driven towards males rather than females.

#____________________________________________________________
# Step 2.2: Perform differential expression analysis “the right way” with DESeq2
dds <- DESeq(dds)

#____________________________________________________________
# Step 2.3: Extract and interpret the results for sex differential expression
res_SEX <- results(dds, name = "SEX_male_vs_female") %>%
  as_tibble(rownames = "GENE_NAME")

res_SEX <- res_SEX %>%
  filter(!is.na(padj)) %>%
  arrange(padj)

data <- read_delim("Week_7_Data/gene_locations.txt")

res_SEX <- left_join(res_SEX, data, by = "GENE_NAME")

## 2.3.2: How many genes exhibit significant differential expression between males and females at a 10% FDR?
### Answer: We can filter through the results_SEX from 2.3 and count the number of rows that are associated to the gene count.
significant_genes <- res_SEX %>%
filter(padj < 0.1)
num_significant_genes <- nrow(significant_genes)
num_significant_genes
### This outputs to a value of 262. Therefore there are 262 genes that exhibit significant differential expression between males and females at a 10% FDR.

## 2.3.3: Examine your top hits. Which chromosomes encode the genes that are most strongly upregulated in males versus females, respectively? 
## Are there more male-upregulated genes or female-upregulated genes near the top of the list? Interpret these results in your own words.
### Answer: After examining the top hits, it can be observed that the Y-chromosome encodes for more genes that are strongly upregulated. After filtering for the most significantly expressed genes by observing the padjusted values, it can be seen that a majority of those genes are associated/correlated to the Y-chromosome which is associated with males.


## 2.3.4: Examine the results for the two genes (WASH7P and SLC25A47) that you had previously tested with the basic linear regression model in step 2.1. Are the results broadly consistent?
### Answer: Overall yes, the results that were outputted for these genes were correlated well. As WASH7P that did not show any significance towards sex-differential expression is associated with chromosome 1, something shared between males and females. While the other gene, SLC25A47 is associated with the X chromosome which would have quite significance Sex association.

#____________________________________________________________
# Step 2.4: Extract and interpret the results for differential expression by death classification
res_DTHHRDY <- results(dds, name = "DTHHRDY_ventilator_case_vs_fast_death_of_natural_causes") %>%
  as_tibble(rownames = "GENE_NAME")

res_DTHHRDY <- res_DTHHRDY %>%
  filter(!is.na(padj)) %>%
  arrange(padj)

res_DTHHRDY <- left_join(res_DTHHRDY, data, by = "GENE_NAME")
## 2.4.1: How many genes are differentially expressed according to death classification at a 10% FDR?
### Answer: There are approxiametly 16069 genes that are differentially expressed according to death classification at a 10% FDR?
DTHHRDY_significant_genes <- res_DTHHRDY %>%
  filter(padj < 0.1)
DTHHRDY_num_significant_genes <- nrow(DTHHRDY_significant_genes)
DTHHRDY_num_significant_genes
## 2.4.2: Interpret this result in your own words. Given your previous analyses, does it make sense that there would be more genes differentially expressed based on type of death compared to the number of genes differentially expressed according to sex?
### Answer: Yes, it would make since for there to be more variance in expression for genes associated with the type of death rather than sex, as there is greater complexity in terms of biology for death and more variability associated compared to sex differentiation between a Male and Female.

#____________________________________________________________
# Exercise 3: Visualization
#____________________________________________________________

###### CORRECTED PLOT w/ ADJ P-VALUE ###########

Sex_Volcano_Plot <- ggplot(data = res_SEX, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = (abs(log2FoldChange) > 2 & padj < 1e-20))) +
  geom_text(data = res_SEX %>% filter(abs(log2FoldChange) > 2 & padj < 1e-50),
            aes(x = log2FoldChange, y = -log10(padj) + 10, label = GENE_NAME), size = 3) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("gray", "coral")) +
  labs(y = expression(-log[10]("adjusted p-value")), x = expression(log[2]("fold change"))) +
  ggtitle("Volcano Plot of Differential Gene Expression by Sex")

Sex_Volcano_Plot

ggsave("Volcano Plot of Differential Gene Expression by Sex.pdf", plot=Sex_Volcano_Plot)
