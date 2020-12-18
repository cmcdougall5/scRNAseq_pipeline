####ReadME
#Script prompts user for a file, then reads the scRNAseq counts in this file, performs qc, transforms and clusters the cells.


#############Housekeeping
#First clear the environment;
rm(list=ls())
#clear all plots
dev.off(dev.list()["RStudioGD"])
#Then the cache
gc()

###########Call libraries
library(svDialogs)
library(dplyr) 
library(Seurat)
library(patchwork)
library(sctransform)
library(ggplot2)

##########Set seed and working directories
#set seed so plots same each time
set.seed(1234)

#set the working directory to directory wish to read data from and save analysis results to
data.dir <- setwd(here::here())


#########Load data
#####Load in the cluster information from the .rds file
data <- readRDS("FletcherFemales.rds")

########QC data
#Find percent mitochondrial dna
#data <- PercentageFeatureSet(data, pattern = "^mt-", col.name = "percent.mt")
#Fletcher uses "percent.Mt"
data <- PercentageFeatureSet(data, pattern = "^mt-", col.name = "percent.MT")
#visualise the QC metrics - use a violin plot 
#VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), ncol = 3, pt.size=0)
#make the plots
#plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.MT")
plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#then plot
plot1 + plot2

#`Can subset data at this point if desired
#data <- subset(data, subset = nFeature_RNA > 200  & percent.mt < 5)

#transform the data (normalise, scale & regress out mt)
#data <- SCTransform(data,vars.to.regress = "percent.mt")#, verbose = FALSE )
data <- SCTransform(data,vars.to.regress = "percent.MT")#, verbose = FALSE )
##########Dimensional reduction

#identify principal components
data <- RunPCA(data,verbose= FALSE) 

#examine which may be dropped
ElbowPlot(data, ndims=30)

#perform dimensional reduction
data <- RunUMAP(data, dims=1:30, UMAP.method = "UMAP-learn", metric = "correlation", verbose=FALSE)


###########Cluster the cells

#Find nearest neighbours     #change dims!
data <- FindNeighbors(data, dims=1:30, verbose=FALSE)
#CLuster based on neighbour distances
data <- FindClusters(data, verbose=FALSE)#, resolution = 0.7 

#visualise clusters
DimPlot(data, label=TRUE) + NoLegend()
#save plot
png("clusters.png")
#visualise the clusters-make plot
DimPlot(data, label=TRUE) + NoLegend()
#save plot 
dev.off()

#########Export the cluster information
saveRDS(data, file = "Clusters.rds")
