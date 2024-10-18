# Exercise 1: Checking fastq quality
## Step 1.1: Open fastqc and load both sample files from the File menu.
### Question: Can you think of a reason why this sample does not match the expected GC content distribution and base content at the beginning of sequences?
#### Answer: 

## Step 1.2: Overrepresented sequences
### Question: What is the origin of the most overrepresented sequence in this sample? Does it make sense?
#### Answer: 

# Exercise 2: Using MultiQC to check processed data quality
### Question: If you were to reject any samples with the percentage of unique reads less than 45%, how many samples would you keep?
#### Answer: 

### Question: Can you see the blocks of triplicates clearly? Try adjusting the min slider up. Does this suggest anything about consistency between replicates?
#### Answer:

# Exercise 3: Normalization and clustering

## Step 3.3: PCA analysis
### Question: What does the PCA plot suggest to you about the data? Why?
#### Answer: The PCA plot prior to corrections shows a lot of variance with the data as there are only a few tissue samples that cluster well, while those like "A2-3, Cu, Fe, LFC-Fe" are spread out. There for these samples for these tissues showed too much variance for them to cluster as effectively as the other tissue samples. After making adjustments and using the vst that stabilizes the variance, then a much better definded clustering is observed between each of the different tissue types and their replicates. Therefor the samples must have distinct genetic profiles to recognize them as distinct clusters and associated to that specific tissue type.

## Step 3.6: Gene ontology enrichment analysis
### Question: Do the categories of enrichments make sense? Why?
#### Answer: Yes, they do make sense. As cluster 1 is observed to be correlated with proteolysis which is the hydrolysis of proteins by cleavage, this makes sense when considering the digestive track that was being analyzed.