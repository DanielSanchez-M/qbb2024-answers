# Step 1: Load Data and Packages
# 1.1: Load Libraries
BiocManager::install( "zellkonverter" )
library("zellkonverter")
library( "DropletUtils" )
library( "scater" )
library( "scran" )
library("ggplot2")

setwd("~/qbb2024-answers/Quant_Bio_Lab/Homework/Week_8")

# 1.2: Load and inspect data
gut <- readH5AD("v2_fca_biohub_gut_10x_raw.h5ad")

## Question 1: Inspect the gut SingleCellExperiment object
assayNames(gut) <- "counts"
gut <- logNormCounts(gut)
gut

### How many genes are quantitated (should be >10,000)?
gut_dimensions <- dim(gut)
num_rows <- dimensions[1] # dimensions[1] would provide the number of rows = number of genes
num_rows
### Answer: There are 13407 genes being quantitated

### How many cells are in the dataset?
num_cols <- dimensions[2] # dimensions[2] would provide the number of columns = number of cells
num_cols
### Answer: There are 11788 cells in the dataset

### What dimension reduction datasets are present?
reducedDimNames(gut)
### Answer: The  datasets of dimension reduction are PCA, TSNE, and UMAP


# 1.3: Inspect cell metadata
## Question 2: Inspect the available cell metadata
### How many columns are there in colData(gut)?
num_col_colData <- ncol(colData(gut))
num_col_colData
### Answer: There are 39 columns in colData(gut)

### Which three column names reported by colnames() seem most interesting? Briefly explain why.
colnames(colData(gut))
### Answer: The three most interesting column names to me are: Tissue, Age, and Sex as these can provide detailed information to drosophila aging or sex dependent tissue development. 


# 1.3.1: Plot cells according to X_umap using plotReducedDim() and colouring by broad_annotation
plotReducedDim(gut, "X_umap", colour_by = "broad_annotation")

#------------------------------------------------------------------------------
# Step 2: Explore data

# 2.1: Explore gene-level statistics
## Question 3: Explore the genecounts distribution
gene_counts <- rowSums(assay(gut))
### What is the mean and median genecount according to summary()? What might you conclude from these numbers?
summary(gene_counts)
### Answer: The mean is 3185 and median is 254. From these numbers, we can conclude that there is a wide distribution of gene counts with potential outliers as the mean is very far from the median.
### What are the three genes with the highest expression after using sort()? What do they share in common?
head(sort(gene_counts, decreasing = TRUE))
### Answer: lncRNA:Hsromega, pre-rRNA:CR45845, and lncRNA:roX1. These are all non-coding RNAs.

# 2.2: Explore cell-level statistics

## Question 4a: Explore the total expression in each cell across all genes
gene_counts <- colSums(assay(gut))
hist(gene_counts)
### What is the mean number of counts per cell?
summary(gene_counts)
### Answer: The mean number of the counts per cell is 3622.
### How would you interpret the cells with much higher total counts (>10,000)?
### Answer: These cells must be very active as they may be actively transcribing and translating genes.

## Question 4b: Explore the number of genes detected in each cell
## Create a vector named celldetected using colSums() but this time on assay(gut)>0
celldetected <- colSums(assay(gut)>0)
## Create a histogram of celldetected using hist()
hist(celldetected)
### What is the mean number of genes detected per cell?
summary(celldetected)
### Answer: The mean number of genes detected per cell is 1059.
### What fraction of the total number of genes does this represent?
### Answer: 1059 out of the total 13407 represents about 8% of total genes.

# 2.3: Explore mitochondrial reads
## Create a vector named mito of mitochondrial gene names using grep() to search rownames(gut) for the pattern ^mt: and setting value to TRUE
mito <- grep("^mt:", rownames(gut), value = TRUE)
## Create a DataFrame named df using perCellQCMetrics() specifying that subsets=list(Mito=mito)
df <- perCellQCMetrics(gut, subsets = list(Mito = mito))
## Confirm that the mean sum and detected match your previous calculations by converting df to a data.frame using as.data.frame() and then running summary()
df <- as.data.frame(df)
summary(df)
## Add metrics to cell metadata using colData(gut) <- cbind( colData(gut), df )
colData(gut) <- cbind(colData(gut), df)
colData(gut)

## Question 5: Visualize percent of reads from mitochondria
## Use plotColData() to plot the subsets_Mito_percent on the y-axis against the broad_annotation on the x-axis rotating the x-axis labels using theme( axis.text.x=element_text( angle=90 ) ) and submit this plot
Mito_percent_per_cell_type <- plotColData(gut, y = "subsets_Mito_percent", x = "broad_annotation") +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(y = "Mitochondrial Percent", x = "Cell Type")
Mito_percent_per_cell_type
ggsave("Mito_percent_per_cell_type.pdf", plot = Mito_percent_per_cell_type)
### Which cell types may have a higher percentage of mitochondrial reads? Why might this be the case?
### Answer: Epithelial cells have the highest percent of mitochondrial reads, this may be due to their rapid turnover in the gut to provide constant new layers in the extreme gut environment.

#------------------------------------------------------------------------------
# Step 3: Identify marker genes
# 3.1: Analyze epithelial cells

## Question 6a: Subset cells annotated as “epithelial cell”
## Create an vector named coi that indicates cells of interest where TRUE and FALSE are determined by colData(gut)$broad_annotation == "epithelial cell"
coi <- colData(gut)$broad_annotation == "epithelial cell"
## Create a new SingleCellExperiment object named epi by subsetting gut with [,coi]
epi <- gut[,coi]
## Plot epi according to X_umap and colour by annotation and submit this plot
epi_umap <- plotReducedDim(epi, "X_umap", colour = "annotation") +
  labs(title = "Epithelial Cells UMAP")
epi_umap
ggsave("Epithelial Cells UMAP.pdf", plot = epi_umap)

## Identify marker genes in the anterior midgut
## Create a list named marker.info that contains the pairwise comparisons between all annotation categories using scoreMarkers( epi, colData(epi)$annotation )
marker.info <- scoreMarkers(epi, colData(epi)$annotation)
## Identify the top marker genes in the anterior midgut according to mean.AUC using the following code
chosen <- marker.info[["enterocyte of anterior adult midgut epithelium"]]
ordered <- chosen[order(chosen$mean.AUC, decreasing=TRUE),]
head(ordered[,1:4])

## Question 6b: Evaluate top marker genes

## What are the six top marker genes in the anterior midgut? Based on their functions at flybase.org, what macromolecule does this region of the gut appear to specialize in metabolizing?
### Mal-A6 carbohydrate metabolic process
### Men-b glucose homeostasis
### vnd transcription factor
### betaTry digestive enzyme
### Mal-A1 carbohydrate metabolic process
### Nhe2 increases intracellular pH
### Answer: Based on their functions, it appears that this region of the gut is specializes in metabolism of carbohydrates.

## Plot the expression of the top marker gene across cell types using plotExpression() and specifying the gene name as the feature and annotation as the x-axis and submit this plot
top_marker_gene_plot <- plotExpression(epi, c("Mal-A6","Men-b","vnd","betaTry","Mal-A1","Nhe2"), x ="annotation") + 
  labs(x = "Cell Type Annotation", y = "Expression Log Fold Change") +
  labs(title = "Expression Levels of Top Marker Genes in the Anterior Midgut") +
  theme(axis.text.x=element_text(angle=90))
top_marker_gene_plot

ggsave("top_marker_gene_plot.png", plot = top_marker_gene_plot)

# 3.2: Analyze somatic precursor cells
## Subset cells with the broad annotation somatic precursor cell and Identify marker genes for intestinal stem cell
somatic_precursors_genes <- colData(gut)$broad_annotation == "somatic precursor cell"
somatic_precursors_genes <- gut[, somatic_precursors_genes]
somatic_precursors_genes

marker.info <- scoreMarkers(somatic_precursors_genes, colData(somatic_precursors_genes)$annotation)
marker.info

chosen <- marker.info[["intestinal stem cell"]]
ordered <- chosen[order(chosen$mean.AUC, decreasing = TRUE),]

## Question 7: Evaluate top marker genes
## Create a vector goi that contains the names of the top six genes of interest by rownames(ordered)[1:6]
goi <- rownames(ordered)[1:6]
goi
## Plot the expression of the top six marker genes across cell types using plotExpression() and specifying the goi vector as the features and submit this plot
somatic_precursors_genes_plot <- plotExpression(somatic_precursors_genes, features = goi, x = "annotation") +
  labs(title = "Expression Levels of Top Marker Genes in Somatic Precursors",
       x = "Cell Type Annotation",
       y = "Expression Level") +
  theme(axis.text.x = element_text(angle = 90))
somatic_precursors_genes_plot

ggsave("somatic_precursors_genes_plot.png", plot = somatic_precursors_genes_plot)

### Which two cell types have more similar expression based on these markers?
### Answer: Enteroblast and Intenstinal Stem Cells

### Which marker looks most specific for intestinal stem cells?
### Answer: DI looks the most specific to intestinal stem cells