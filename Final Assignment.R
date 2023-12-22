# Load necessary libraries
library(DESeq2)            # For differential expression analysis
library(survival)          # For survival analysis
library(glmnet)            # For fitting generalized linear models
library(ggplot2)           # For data visualization
library(clusterProfiler)   # For pathway enrichment analysis
library(factoextra)        # For extracting and visualizing the results of multivariate data analyses
library(AnnotationDbi)     # For mapping gene symbols to Entrez IDs
library(org.Hs.eg.db)      # For mapping gene symbols to Entrez IDs

# Set directory (modify as needed)
path  = "C:/Users/Nelso/OneDrive/UCD/Bio Principles/R"            # Change this to your own directory
file_name = "brca_tcga_pan_can_atlas_2018.tar.gz"                 # Enter file name
untar(file_name)                                                  # Extract the files into folders
setwd(paste(getwd() , "/brca_tcga_pan_can_atlas_2018", sep = "")) # change directory to the extracted folders

# Reading the data files
clinical <- read.delim("data_clinical_patient.txt") # Reads patient clinical data into 'clinical'
rnaseq <- read.delim("data_mrna_seq_v2_rsem.txt")   # Reads RNASeq data into 'rnaseq'
cna <- read.delim("data_cna.txt")                   # Reads Copy Number Aberrations data into 'cna'

# Processing RNASeq data: Deleting genes with more than one Hugo Symbol
keep <- !duplicated(rnaseq[,1])                     # Identifies unique genes in 'rnaseq'
rnaseq <- rnaseq[keep,]                             # Keeps only unique genes
rownames(rnaseq) <- rnaseq[,1]                      # Sets gene names as rownames in 'rnaseq'

# Finding ERBB2 in CNA data and matching with RNASeq data
erbb2_indx <- which(cna[,1] == 'ERBB2')             # Finds the index of ERBB2 in 'cna'
rna_cna_id <- which(colnames(rnaseq[-c(1,2)]) %in% colnames(cna[-c(1,2)])) # Matches RNASeq and CNA patient IDs
rna_cna_sub <- rnaseq[, 2 + rna_cna_id]             # Subsets RNASeq data to matched patient IDs

# Creating DESeq2 dataset
countData <- round(as.matrix(rna_cna_sub))          # Converts RNASeq data to a matrix and rounds it
colnames(countData) <- colnames(rna_cna_sub)        # Sets the column names of countData to match 'rna_cna_sub'
erbb2_cna_levels <- as.numeric(cna[erbb2_indx, colnames(countData)]) # Extracts ERBB2 CNA levels
metadata <- data.frame(row.names = colnames(countData), 
                       erbb2_amplified = erbb2_cna_levels > 0) # Creates a metadata dataframe with ERBB2 amplification status
metadata$erbb2_amplified <- as.factor(metadata$erbb2_amplified) # Converts ERBB2 amplification status to a factor
dds <- DESeqDataSetFromMatrix(countData = countData, colData = DataFrame(metadata), design = ~ erbb2_amplified) # Creates a DESeq2 dataset

# DESeq2 normalization and obtaining differentially expressed genes
dds <- DESeq(dds)                                   # Runs DESeq2 normalization
res <- results(dds, contrast = c("erbb2_amplified", "TRUE", "FALSE")) # Extracts results for the specified contrast
resOrdered <- res[order(res$padj), ]                # Most significantly differentially expressed genes
top10 <- head(resOrdered, 10)                       # Extracts top 10 significantly differentially expressed genes
resOrdered_up <- res[order(res$log2FoldChange, decreasing = TRUE),] # Orders results in decreasing order by log2 fold change
top10_genes_up <- head(resOrdered_up, 10)                 # Extracts top 10 upregulated differentially expressed genes
resOrdered_down <- res[order(res$log2FoldChange, decreasing = FALSE),] # Orders results in increasing order by log2 fold change
top10_genes_down <- head(resOrdered_down, 10)               # Extracts top 10 downregulated differentially expressed genes

## Save DESeq2 result
write.csv(top10, "top10_mostsignificant.csv")
write.csv(top10_genes_up, "top10_upregulated.csv")
write.csv(top10_genes_down, "top10_downregulated.csv")


# Pathway Enrichment Analysis
gene_symbols <- rownames(resOrdered)[which(resOrdered$padj < 0.05)] # Extracts Entrez IDs of significant genes
entrez_ids <- mapIds(org.Hs.eg.db,           # Map gene symbols to Entrez IDs
                     keys = gene_symbols,
                     column = "ENTREZID",
                     keytype = "SYMBOL",
                     multiVals = "first")
entrez_ids <- entrez_ids[!is.na(entrez_ids)] # Clean up NA values that result from genes that could not be mapped
enrichResult <- enrichKEGG(gene = entrez_ids, organism = 'hsa', pvalueCutoff = 0.05) # Performs KEGG pathway enrichment analysis

## Save Pathway Enrichment Result
write.csv(enrichResult, "PathwayEnrichmentResult.csv")

# PCA Plot
vsd <- vst(dds, blind = FALSE)                      # Performs variance stabilizing transformation
pcaResult <- prcomp(t(assay(vsd)))                  # Computes PCA on transformed data
pcaData <- as.data.frame(pcaResult$x[, 1:2])        # Extracts first two principal components
sampleDists <- dist(t(assay(vsd)))                  # Calculates the distance matrix for clustering
sampleClustering <- hclust(sampleDists, method = "average") # Performs hierarchical clustering
clusters <- cutree(sampleClustering, k = 5)         # Cuts the dendrogram to create clusters
pcaData$cluster <- factor(clusters)                 # Adds cluster assignments to PCA data
ggplot(pcaData, aes(x = PC1, y = PC2, color = cluster)) + geom_point() + theme_minimal() + ggtitle("PCA Plot with Clustering") # Plots PCA with clusters


# Differential Expression Analysis Between Clusters
dds$cluster <- factor(clusters)                     # Adds cluster information to DESeq2 object
design(dds) <- formula(~ cluster)                   # Updates design formula for differential expression analysis
dds <- DESeq(dds)                                   # Runs DESeq2 with new design
resClusterComparison <- results(dds, contrast=c("cluster", "1", "2")) # Extracts differential expression results for cluster comparison
top_genes_cluster <- head(resClusterComparison[order(resClusterComparison$pvalue), ]) # Extracts top genes based on p-value in cluster comparison

## Save Differential Expression Analysis Between Clusters Result
write.csv(top_genes_cluster, "topgenescluster.csv")
