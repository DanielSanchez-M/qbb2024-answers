# Exercise 3
# Step 3.1: Loading and filtering the data

# Install Packages and Libraries
BiocManager::install("DESeq2")
BiocManager::install("vsn")
install.packages("hexbin")
install.packages("ggfortify")

library(DESeq2)
library(vsn)
library(matrixStats)
library(readr)
library(dplyr)
library(tibble)
library(hexbin)
library(ggfortify)

sessionInfo()

# Step 3.1: Loading and filtering the data
# Load the salmon.merged.gene_counts.tsv and save it in a data frame
dataset = readr::read_tsv("salmon.merged.gene_counts.tsv")

# Switch gene_name from column to row and remove column_id
dataset = column_to_rownames(dataset, var="gene_name")
dataset = dataset %>% dplyr::select(-gene_id)

# Convert the numbers to integers to remove decimals
dataset = dataset %>% dplyr::mutate_if(is.numeric, as.integer)

# Filter to remove rows with less than 100 reads
dataset = dataset[rowSums(dataset) > 100,]

# Select only the “narrow region” samples, i.e. those covering midgut sections A1, A2-3, Cu, LFC-Fe, Fe, P1, and P2-4.
narrow_set = dataset %>% dplyr::select("A1_Rep1":"P2-4_Rep3")

# Step 3.2: Creating DESeq2 model and batch-correction
# Create a metadata tibble, containing two columns: one with sample names and the other with replicate number
narrow_metadata = tibble(tissue=c("A1", "A1", "A1",
                           "A2-3", "A2-3", "A2-3",
                           "Cu", "Cu", "Cu",
                           "LFC-Fe", "LFC-Fe", "LFC-Fe",
                           "Fe", "Fe", "Fe",
                           "P1", "P1", "P1",
                           "P2-4", "P2-4", "P2-4"),
                  rep=c("Rep1", "Rep2", "Rep3",
                        "Rep1", "Rep2", "Rep3",
                        "Rep1", "Rep2", "Rep3",
                        "Rep1", "Rep2", "Rep3",
                        "Rep1", "Rep2", "Rep3",
                        "Rep1", "Rep2", "Rep3",
                        "Rep1", "Rep2", "Rep3"))
# Create a DESeq2 object using the command DESeqDataSetFromMatrix
ddsNarrow = DESeqDataSetFromMatrix(countData = narrow_set, colData = narrow_metadata,
                                  design = ~tissue)

# Use the function vst (variance stabilizing transformation)
vstNarrow = vst(ddsNarrow)
# Plot your data using meanSdPlot, convert your data inside the function using assay
meanSdPlot(assay(vstNarrow))

# Step 3.3: PCA Analysis
# Perform PCA on the corrected data and plot it. Save the plot to turn in.
#Look at mean by variance
meanSdPlot(assay(ddsNarrow))
#log transform of data
logNarrow = normTransform(ddsNarrow)
# look at log mean by variance
meanSdPlot(assay(logNarrow))
# Perform PCA
log_Narrow_PCA_data = plotPCA(logNarrow,
                     intgroup=c('rep', 'tissue'), returnData=TRUE)
Narrow_PCA_Plot = ggplot(log_Narrow_PCA_data, aes(PC1, PC2, color=tissue, shape=rep)) +
  geom_point(size=5) +
  ggtitle("PCA of Narrow Regions of the Drosophila melanogaster midguts")
ggsave("Narrow_PCA_Data.pdf", plot=Narrow_PCA_Plot)

# Fix the dataset due to variance issues (Use the vst narrow dataset instead of log Narrow dataset)
log_Narrow_PCA_data = plotPCA(vstNarrow,
                              intgroup=c('rep', 'tissue'), returnData=TRUE)
Narrow_PCA_Plot = ggplot(log_Narrow_PCA_data, aes(PC1, PC2, color=tissue, shape=rep)) +
  geom_point(size=5) +
  ggtitle("PCA of Narrow Regions of the Drosophila melanogaster midguts")
ggsave("Fixed_Narrow_PCA_Data.pdf", plot=Narrow_PCA_Plot)

# Step 3.4: Filtering genes by variance
# Convert the data into a matrix with as.matrix(assay(data))
vst_matrix = as.matrix(assay(vstNarrow))

combined = vst_matrix[,seq(1,21,3)] #Start position 1 go to 21, step by 3
combined = combined + vst_matrix[,seq(2,21,3)] #Start position 2 go to 21, step by 3
combined = combined + vst_matrix[,seq(3,21,3)] #Start position 3 go to 21, step by 3
combined = combined / 3 #Average across the replicated

filtered = rowSds(combined) > 1

# Step 3.5: K-means clustering genes
# Create a matrix with the narrow filtered sample data from 3.4
Narrow_matrix = vst_matrix[filtered,]
# Create the heatmap from the narrow matrix filtered data
heatmap(Narrow_matrix)
# remove the top tree branche
heatmap(Narrow_matrix, Colv=NA)

# Refine the heatmap
# Set seed to 42 for random generator start point
set.seed(42)
# Assign and set clusters as well as ordering them
k = kmeans(Narrow_matrix, centers=12)$cluster
ordering = order(k)
k = k[ordering]
Narrow_matrix = Narrow_matrix[ordering,]
# Create the heatmap along with color coded cluster assignments, save it as a pdf
pdf("Narrow_Heatmap.pdf") #opens a blank pdf file to write the heatmap code into
Narrow_heatmap = heatmap(Narrow_matrix, Rowv=NA, Colv=NA, RowSideColors = RColorBrewer::brewer.pal(12,"Paired")[k])
dev.off() #closes file editor


# Step 3.6: Gene ontology enrichment analysis
genes = rownames(Narrow_matrix[k==1,])
write.table(genes, 'narrow_cluster1.txt', sep="\n", quote=FALSE, row.names=FALSE,col.names=FALSE)
# In README.md File and downloaded GO Analysis Report



