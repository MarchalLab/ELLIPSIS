#libraries
requiredPackages = c('dplyr', 'Seurat', 'fossil', 'pheatmap', 'optparse', 'readxl', 'plyr', 'stringr')
packagesToInstall = requiredPackages[! requiredPackages %in% installed.packages()[,"Package"]]
if (length(packagesToInstall) > 0)
  install.packages(packagesToInstall, repos="https://ftp.belnet.be/mirror/CRAN/")

library(dplyr)
library(Seurat)
library(fossil)
library(pheatmap)
library(optparse)
library(readxl)
library(plyr)
library(stringr)

#command line options
option_list = list(
  make_option(c("-c","--countFile"), type="character", default=NULL, 
              help="file containing the gene level counts per sample", metavar="character"), 
  make_option(c("-g","--gtfFile"), type="character", default=NULL,
              help="path/to/gtfFile", metavar="character"),
  make_option(c("-o","--outFile"), type="character", default=NULL,
              help="output file to which KNN neighbors are written", metavar="character"),
  make_option(c("-n","--numNB"), type="integer", default=20,
              help="number of neighbors per cell [default = 20]", metavar="character"), 
  make_option(c("-m", "--minGenes"), type="integer", default=200,
              help="minimum number of expressed genes in order to consider a cell [default = 200]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#make sure all required options are filled in
if (is.null(opt$countFile)){
  stop("Please provide the count file", call.=FALSE)
}
if (is.null(opt$gtfFile)){
  stop("Please provide the gtf file", call.=FALSE)
}
if (is.null(opt$outFile)){
  stop("Please provide the output file name", call.=FALSE)
}

#input
countFile <- opt$countFile
gtfFile <- opt$gtfFile
outFile <- opt$outFile
numNB <- opt$numNB
minGenes <- opt$minGenes

#for reproducability
set.seed(42)

# Load raw counts
counts <- read.table(countFile)
counts = counts[!startsWith(rownames(counts),"__"),] #remove non-gene rows

# Initialize the Seurat object with the raw (non-normalized data).
# already remove genes expressed in < 3 cells & cells having expression in < 200 genes
glio <- CreateSeuratObject(counts = counts, min.cells = 3, min.features = minGenes)

print(paste0("Total number of cells: ", ncol(counts)))
print(paste0("Removing cells with expression in < ", minGenes, " genes: ", ncol(glio), " cells remaining"))

# mitochondrial mapping
gtf = read.table(gtfFile, 
                 col.names = c("chr", "source", "type", "start", "stop", ".A", "strand", ".B", "info"), 
                 sep = "\t")
gtf = gtf[(gtf$chr == "chrM") & (gtf$type == "gene"), ]
gtf$geneID = str_replace(gtf$info, ".*gene_id ", "")
gtf$geneID = str_replace(gtf$geneID, ";.*", "")
mitGenes <- gtf$geneID[gtf$geneID %in% rownames(glio)]
glio[["percent.mt"]] <- PercentageFeatureSet(glio, features = mitGenes)

print("No additional filtering on % mitochrondrial reads, nr of genes espressed or total nr of reads")
print("uncomment lines 77-84 and change to preferred filter thresholds.")
# # Visualize QC metrics as a violin plot
# VlnPlot(glio, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# 
# # filtering -> based on previous plots
# glio <- subset(glio, subset = percent.mt < 25)
# print(paste0("Removing cells with more than 25% mitochondrial reads: ", ncol(glio), " cells remaining"))
# glio <- subset(glio, subset = nFeature_RNA > 500 & nCount_RNA > 10000)
# print(paste0("Removing cells with less than 500 expressed genes or less than 10 000 reads: ", ncol(glio), " cells remaining"))

# normalization
glio <- NormalizeData(glio, normalization.method = "LogNormalize", scale.factor = 10000)

# identify highly variable features
glio <- FindVariableFeatures(glio, selection.method = "vst", nfeatures = 2000)

# scaling (mean expression accross cells 0 & variance 1)
glio <- ScaleData(glio, vars.to.regress = c("percent.mt"))

#run PCA
glio <- RunPCA(glio, npc = 30)

glio_neighborhood <- FindNeighbors(glio, reduction = "pca", return.neighbor=TRUE)

#write top 20 neighbors to file
top20Neighbors <- NULL
for (cell in colnames(glio)){
  top20Neighbors <- rbind(top20Neighbors, TopNeighbors(glio_neighborhood[["RNA.nn"]] ,cell = cell, n=numNB))
}
write.table(top20Neighbors, 
            file = outFile,  
            row.names = F, col.names = F, sep = "\t", quote = F)


