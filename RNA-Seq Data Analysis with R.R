##################### RNA-Seq Data Analysis with R #########################

##################### Installing and Loading Packages #####################

# Install required packages

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinks") 

# TCGAbiolinks is installed from Bioconductor because it is not available on CRAN.

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

# DESeq2 is installed from Bioconductor because it is not available on CRAN.

install.packages("dplyr")
install.packages("pheatmap")

# Load required libraries

library(TCGAbiolinks)
library(DESeq2)
library(dplyr)
library(pheatmap)

############ Exploring TCGA-LIHC Project Information #######################

TCGAbiolinks::getProjectSummary("TCGA-LIHC")

# Retrieve summary information about the TCGA-LIHC project.

############ Querying, Downloading, and Preparing RNA-Seq Data ################

## Query RNA-Seq Data ##

# GDCquery() is used to query RNA-seq data from the GDC portal.

query_TCGA = GDCquery(
  project = "TCGA-LIHC",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts")

# Retrieve query results.

lihc_res = getResults(query_TCGA) 

# Display column names of the query results.

colnames(lihc_res) 

# Check available sample types.

unique(lihc_res$sample_type) 

# Summarize sample types in the dataset.

table(lihc_res$sample_type)   

# Re-query using only primary tumor and normal tissue samples.

query_TCGA = GDCquery(
  project = "TCGA-LIHC",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal"))

## Download RNA-Seq Data ## 

# Download queried RNA-seq data.

GDCdownload(query = query_TCGA)

## Prepare RNA-Seq Data ## 

# Prepare and load downloaded RNA-seq data.

tcga_data = GDCprepare(query_TCGA)

# Check dimensions of the prepared dataset.

dim(tcga_data)

## Count Data (Gene Expression Matrix) ## 

# Extract count matrix from the prepared dataset. 

count_data <- as.data.frame(assay(tcga_data))

# Check dimensions of count data.

dim(count_data)

## Sample Information ##  

# Extract clinical/sample metadata.

sample_info <- as.data.frame(colData(tcga_data))

# Check survival status distribution.

table(sample_info$vital_status)

# Check gender distribution.

table(sample_info$gender)

# Check race distribution.

table(sample_info$race)

## Gene Annotation Information ##

# Extract gene annotation information.

gene_mapping <- as.data.frame(rowData(tcga_data))

# Display the first six gene annotations.

head(gene_mapping)

######################## Saving Data #######################################

## Save Specific Object ##

# Save tcga_data object as an RDS file.

saveRDS(object = tcga_data,
        file = "tcga_data.RDS",
        compress = FALSE)

# Load saved RDS object.

tcga_data = readRDS(file = "tcga_data.RDS")

## Save Entire Workspace ## 

# Save all current R objects.

save.image(file = 'TCGA.RData')

## Save Script ## 

# Save the current R script.

# File < Save As 

# CTRL + S

############### Differential Expression Analysis with DESeq2 ############### 

## Preparing DESeq2 Object ## 

# Check whether sample names match between count data and metadata.

all(colnames(count_data) %in% rownames(sample_info))

# Check whether sample names are in the same order.

all(colnames(count_data) == rownames(sample_info))

# Match column names and row names if necessary.

match(colnames(count_data), rownames(sample_info))

# Convert grouping variable into a factor.

sample_info$shortLetterCode <- factor(sample_info$shortLetterCode)

## Create DESeqDataSet Object ## 

# Create DESeqDataSet object for DESeq2 analysis.

dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = sample_info,
                              design = ~shortLetterCode)

dds

# Filter out low-count genes.

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Set reference level for comparison.

dds$shortLetterCode <- relevel(dds$shortLetterCode, ref = "NT")

################# Running Differential Expression Analysis ###########

# Run DESeq2 analysis.

dds <- DESeq(dds)

saveRDS(object = dds,
        file = "dds_res.RDS",
        compress = FALSE) 

# Save DESeq2 object.

# Extract DESeq2 results.

DESeq_result <- results(dds)

# Convert DESeq2 results into a dataframe.

DESeq_result <- as.data.frame(DESeq_result)
head(DESeq_result)

# Sort results by p-value.

DESeq_result_ordered <- DESeq_result[order(DESeq_result$pvalue),]
head(DESeq_result_ordered)

# Filter significantly differentially expressed genes.

filtered <- filter(DESeq_result, padj < 0.05, abs(log2FoldChange) > 1)

# Separate upregulated and downregulated genes.

up_regulated <- subset(filtered, log2FoldChange > 1)
down_regulated <- subset(filtered, log2FoldChange < -1)

# Map Ensembl IDs to gene symbols.

filtered_name <- rownames(filtered)
filtered_name <- subset(gene_mapping, gene_id %in% filtered_name, select = gene_name)

filtered_up_name <- rownames(up_regulated)
filtered_up_name <- subset(gene_mapping, gene_id %in% filtered_up_name, select = gene_name)

rownames(up_regulated) <- filtered_up_name$gene_name
rownames(up_regulated)<- make.names(filtered_up_name$gene_name, unique=TRUE)

# Save filtered results to a text file.

write.table(filtered, file = "filtered.txt", sep = "\t",row.names = TRUE)

#################### PCA Plot Generation ##############################

# Perform variance stabilizing transformation (VST).

vsd <- vst(dds, blind = FALSE)

# Generate PCA plot using transformed data.

plotPCA(vsd, intgroup = c("shortLetterCode"))


######################### Heatmap Visualization ##############################

# Select top 10 significant genes.

top_hits <- DESeq_result[order(DESeq_result$padj),][1:10,]
top_hits <- rownames(top_hits)
top_hits

# Retrieve gene symbols for top genes.

top_hits_gename <- subset(gene_mapping, gene_id %in% top_hits, select = gene_name)

# Match gene names with selected genes.

top_hits_gename <- top_hits_gename[match(top_hits, rownames(top_hits_gename)),]

# Generate heatmap using normalized expression values.

pheatmap(assay(vsd)[top_hits,1:6],cluster_rows = FALSE, show_names=TRUE, cluster_cols=FALSE, labels_row = top_hits_gename)

# Generate clustered heatmap.

pheatmap(assay(vsd)[top_hits,1:6], cluster_cols = FALSE)

# Add sample annotations to heatmap.

annot_info <- as.data.frame(colData(dds)[,c('shortLetterCode','gender')])

pheatmap(assay(vsd)[top_hits,1:6],cluster_rows = FALSE, show_rownames = TRUE, 
         cluster_cols=FALSE, annotation_col = annot_info)

#########################################################################################
#########################################################################################





