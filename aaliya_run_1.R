# load packages----
library(rhdf5) #provides functions for handling hdf5 file formats (kallisto outputs bootstraps in this format)
library(tidyverse) # provides access to Hadley Wickham's collection of R packages for data science, which we will use throughout the course
library(tximport) # package for getting Kallisto results into R
library(ensembldb) #helps deal with ensembl
library(EnsDb.Hsapiens.v86) #replace with your organism-specific database package
targets <- read_tsv("design_1.txt")
path <- file.path(targets$sample, "abundance.tsv")
all(file.exists(path)) 
Tx <- transcripts(EnsDb.Hsapiens.v86, columns=c("tx_id", "gene_name"))
Tx <- as_tibble(Tx)
Tx <- dplyr::rename(Tx, target_id = tx_id)

Tx <- dplyr::select(Tx, "target_id", "gene_name")
Txi_gene <- tximport(path, 
                     type = "kallisto", 
                     tx2gene = Tx, 
                     txOut = FALSE, #determines whether your data represented at transcript or gene level
                     countsFromAbundance = "lengthScaledTPM",
                     ignoreTxVersion = TRUE)
library(DESeq2)
dds <- DESeqDataSetFromTximport(txi = Txi_gene,
                                colData = targets,
                                design = ~ condition)
dds <- DESeq(dds)
results <- results(dds)
filteredResults <- results[which(abs(results$log2FoldChange) > 1 & results$pvalue < 0.05),]

write.csv(as.data.frame(filteredResults_1), file = "DESeq2_filtered_results_1.csv")
install.packages("openxlsx")
library(openxlsx)
filteredResults <- read.csv("DESeq2_filtered_results.csv")
write.xlsx(filteredResults, file = "DESeq2_filtered_results.xlsx")

#now PCA analysis 
library(DESeq2)
dds <- DESeq(dds) # assuming 'dds' is your DESeqDataSet
vst_data <- vst(dds, blind=FALSE)

pca_res <- plotPCA(vst_data, intgroup=c("condition"), returnData=TRUE)

library(ggplot2)
ggplot(pca_res, aes(x=PC1, y=PC2, color=condition)) + geom_point() + 
  theme_minimal() + ggtitle("PCA of Treated vs Control")

#to calculate % variance 

# Assuming `vst_data` is your variance stabilized transformation data from DESeq2
pca <- prcomp(assay(vst_data))

# Calculate the explained variance
explained_variance <- pca$sdev^2 / sum(pca$sdev^2) * 100

# Create a data frame for plotting
explained_df <- data.frame(PC = 1:length(explained_variance), Variance = explained_variance)

# Now, `explained_df` should have the correct data for plotting or analysis.

library(ggplot2)

ggplot(explained_df, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  theme_minimal() +
  labs(title = "PCA Explained Variance", x = "Principal Component", y = "Variance Explained (%)")




