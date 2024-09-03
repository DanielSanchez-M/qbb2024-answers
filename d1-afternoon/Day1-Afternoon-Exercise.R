#Q1
library(tidyverse)
library(ggthemes)
df <- read_delim("~/Data/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")

#Q2
df
glimpse(df)

#Q3
RNA_seq <- df %>%
  filter(SMGEBTCHT == "TruSeq.v1")

#Q4
ggplot(data = RNA_seq) +
  geom_bar(mapping = aes(x = SMTSD)) +
  xlab("Tissue Type") +
  ylab("Frequency") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#Q5
ggplot(data = RNA_seq) +
  geom_histogram(mapping = aes(x = SMRIN), bins = 10) +
  xlab("RNA Integrity Number") +
  ylab("Cell Frequency")
##Answers to Q5
  #Best plot for visualization = Histogram
  #Shape of distribution = Unimodal, slightly skew to the right

#Q6
ggplot(data = RNA_seq, mapping = aes(x = SMTSD, y = SMRIN)) +
  geom_violin() +
  xlab("Tissue Type") +
  ylab("RNA Integrity Number") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
##Answer to Q6
# A violin plot would be best to show contrasting distributions across multiple groups
# Compared to other Tissue Types, "Cells" have a collective higher distribution having a high RNA Integrity Number, can be due to their controlled environment
# Outliers = Kidney Medulla, reported a much smaller distribution with an RIN around a 6, compared to the rest of the samples that have a larger distributed RIN.

#Q7
ggplot(data = RNA_seq, mapping = aes(x = SMTSD, y = SMGNSDTC)) +
  geom_violin() +
  xlab("Tissue Type") +
  ylab("Number of Genes") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
##Answer for Q7
#Nearly all genes have a similar distribution for the number of genes present in the tissue sample. However, certain tissue types like the Cervix, Fallopian Tubes, and Kidneys have a smaller distribution of the number of genes but are still concentrated within the same range of the number of genes as the other tissue types.
#One major outlier is the Testis. This is due to them expressing the largest number of genes of any mammalian organ

#Q8
ggplot(data = RNA_seq, mapping = aes(x = SMTSISCH, y = SMRIN)) +
  geom_point(size = 0.5, alpha = 0.5) +
  xlab("Ischemic Time") +
  ylab("RNA Integrity Number") +
  facet_wrap("SMTSD") +
  geom_smooth(method = "lm")
##Answer to Q8
#For most tissue types, RIN does not change over the Ischemic Time, however some tissue types do go down in RIN such as Liver and Fallopian Tubes
#The relationship does depend on the tissue type

#Q9
ggplot(data = RNA_seq, mapping = aes(x = SMTSISCH, y = SMRIN)) +
  geom_point(size = 0.5, alpha = 0.5, aes(color = SMATSSCR)) +
  xlab("Ischemic Time") +
  ylab("RNA Integrity Number") +
  facet_wrap("SMTSD") +
  geom_smooth(method = "lm")
##Answer to Q9
#The longer the Ischemic Time, the lower the RNA Integrity Number, and the Higher the Autolysis
#Yes, the relationship does depend on tissue as there are some that are more stable than others. Ex: Muscular Tissue is relatively stable with low autolysis
