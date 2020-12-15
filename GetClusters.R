####ReadME
#Script prompts user for a file, then reads the scRNAseq counts in this file, performs qc, transforms and clusters the cells.


#############Housekeeping
#First clear the environment;
rm(list=ls())
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
#use a dialog for user to input the folder containint the data
user.data <- dlgInput("Enter name of folder containing data to be analysed", Sys.info()["user"])$res

#load the data
#load raw data using Read10x, this will be relative to your path, name desired folder e.g.  "Cheung Data"
raw.data <- Read10X(user.data)

#initialise a Seurat object with the raw (not normalised) data. This is a count matrix
data <- CreateSeuratObject(counts=raw.data, project = "data", min.cells = 3, min.features = 200) 


########QC data
#Find percent mitochondrial dna
data <- PercentageFeatureSet(data, pattern = "^mt-", col.name = "percent.mt")
#Fletcher uses "percent.Mt"
#data <- PercentageFeatureSet(data, pattern = "^mt-", col.name = "percent.MT")
#visualise the QC metrics - use a violin plot 
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)
#make the plots
plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#then plot
plot1 + plot2

#`Can subset data at this point if desired
#data <- subset(data, subset = nFeature_RNA > 200  & percent.mt < 5)

#transform the data (normalise, scale & regress out mt)
data <- SCTransform(data,vars.to.regress = "percent.mt")#, verbose = FALSE )

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

#save plot
png("clusters.png")
#visualise the clusters-make plot
DimPlot(data, label=TRUE) + NoLegend()
#save plot 
dev.off()

#########Export the cluster information
saveRDS(data, file = "../scRNAseq_pipeline/Clusters.rds")
