# scRNAseq data analysis

############## precursor, set the directories to save to and read data from:

#set the working directory to directory wish to save analysis results to
setwd("/Users/craigmcdougall/Documents/Project 2/Data analysis/Cheung")

#check working directory
getwd()


#define the data directory where the data is held
data.dir <- "/Users/craigmcdougall/Documents/Project 2/Data/Cheung"

#list all the files in the data directory: barcodes, genes, matrix
list.files(data.dir)
#have now seen that the files are there

############################# 1. Load libraries
#libraries required throughought the script
library(dplyr)
library(Seurat)
library(patchwork)
library(sctransform)
library(ggplot2)
library(scales)

############################ 2. Load the data  

#load raw data using Read10x 
raw.data <- Read10X(data.dir=data.dir)

#initialise a Seurat object with the raw (not normalised) data. This is a count matrix
data <- CreateSeuratObject(counts=raw.data, project = "data", min.cells = 3, min.features = 200) 

######################### 3. Examine confounding sources of variation (mitochondrial mapping)

#QC: selecting cells for analysis
  #- # genes per cell (id low quality cells or multiple reads)
  #- total # molecules in a cell
  #- % of reads that map to mitochondrial genome

#store mitochoindrial % in the Seurat object meta data
#Seurat recommends MT- prefix for mitochondria data
#cheung, ho, mayran use mt
data <- PercentageFeatureSet(data, pattern = "^mt-", col.name = "percent.mt")
#Fletcher uses Mt-
#data <- PercentageFeatureSet(data, pattern = "^mt-", col.name = "percent.MT")

#use QC metrics to filter the cells
#visualise the QC metrics - use a violin plot 
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#make the plots
plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#then plot
plot1 + plot2

############################4. Apply SCTransform normalisation (Filter, normalise, regress & detect variable genes)
#SCTransform replaces Normalizedata, ScaleData and FindVariableFeature functions in older revisions of script

#use QC plots to inform which cells should be excluded below
#SCTransform parameters:
  #Filter out cells which have less than this many genes expressed [200]
  #Filter out cells which have higher unique gene count [2500]
  #Filter out cells which have higher mitochondrial transcript percentage [5]
  #Regress out cell cycle differences [no]
  #Number of variable features to return [3000]


#run sctransform, remove mitochondrial data by regressing it out
data <- SCTransform(data,vars.to.regress = "percent.mt")#, verbose = FALSE )


##################### 5. Dimensional reduction 

#PCA to identify PCs
#run PCA following new Seurat standard steps:https://satijalab.org/seurat/v3.1/sctransform_vignette.html
data <- RunPCA(data,verbose= FALSE) 

#Do we need to keep all PCs, or is majority of data captured by a certain point?
#ElbowPlot(data, ndims=30)

#UMAP dimensional reduction of the identified PCs
data <- RunUMAP(data, dims=1:30, umap.method = "umap-learn", metric = "correlation", verbose=FALSE) #note because SCTransform passes 3000 features can run more PCs


#################### 6. Clustering
#Find nearest neighbours
data <- FindNeighbors(data, dims=1:30, verbose=FALSE)

#CLuster based on neighbour distances
data <- FindClusters(data, verbose=FALSE)

#plot the clusters
DimPlot(data, label=TRUE) + NoLegend()

#######################   7. Finding deferentially expressed features (biomarkers for each)

#find the markers that define clusters via differential expression 
#what makes each cluster different?

#use findAllMarkers, 25% de threshold, only +ve ones
#data.markers <-FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#can look to see top 2 DE genes in each cluster
#increase display optiions to show all
#options(dplyr.print_max=1e9)
#data.markers %>% group_by(cluster) %>% top_n(n=2, wt= avg_logFC)


#########################  8. attributing cell type
#how well do the clusters approximate the cell types?
#This is done manually by examination of feature plot 
#Done by cell expression, classification criteria (only really care about corticotrophs, but other plots usefull to exclude cluster)
#gene information from FLetcher & cheung papers:
# Thyrotrophs (Pou1f1 + Tshb),
# Somatotrophs(Pou1f1 + Gh)
# Lactotrophs (Pou1f1 + Prl)
# Gonadotrophs (Lhb)
# Melanotrophs (Pomc + Pcsk2 + Pax7)
#Corticotrophs (Pomc + Crhr1 + Avpr1b - Pcsk2 - Pax7)

#Visualise the canonical marker genes 

#first exclude Pou1f1 expressing cells (Thyrotrophs, SOmatotrophs, Lactotrophs)
FeaturePlot(data, features = c("Pou1f1","Tshb", "Gh","Prl"), pt.size = 0.2)&scale_color_viridis_c() 

#then exclude the gonadotrophs
FeaturePlot(data, features = c("Lhb"), pt.size = 0.2)&scale_color_viridis_c() 

#then exclude the melanotrophs
FeaturePlot(data, features = c("Pomc","Pax7","Pcsk2"), pt.size = 0.2)&scale_color_viridis_c() 

#then id the corticotrophs
FeaturePlot(data, features = c("Pomc", "Crhr1","Avpr1b","Gpc5"), pt.size = 0.2)&scale_color_viridis_c() 

#for Cheung data 
# Melanotrophs cluster 16
#Corticotrophs  cluster 14

# can put names on the featureplot 
#enter the cotricoph cell ids, found in visualisation. Note change to correct cluster #s or will miss-match
#new.cluster.ids <- c("0", "1", "2", "3", "4", "5", "6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22")
new.cluster.ids <- c("0", "1", "2", "3", "4", "5", "6","7","8","9","10","11","12","13","cort","15","16","17","18","19","20","21","22")

 
names(new.cluster.ids) <- levels(data)
 
data <- RenameIdents(data, new.cluster.ids)
#then put on the plot to confirm slected the correct one (sanity check)
DimPlot(data, reduction ="umap", label = TRUE, pt.size = 0.5)+NoLegend()

####################### 9. pull out the corticotrophs

#or pull out a  seurat of cort only
cortico <- subset(data, idents = "cort")

##################### 10. Re-cluster the corticotrophs only

#dimensional analysis
#PCA to identify PCs
cortico <- RunPCA(cortico,verbose= FALSE) 

#dimensional reduction
#UMAP dimensional reduction of the identified PCs
cortico <- RunUMAP(cortico, dims=1:30, verbose=FALSE) #note because SCTransform passes 3000 features can run more PCs

#Find nearest neighbours
cortico <- FindNeighbors(cortico, dims=1:30, verbose=FALSE)

#CLuster based on neighbour distances
cortico <- FindClusters(cortico, verbose=FALSE)

#plot the clusters
DimPlot(cortico, pt.size = 5, label=TRUE) + NoLegend()

#Visualise the canonical marker genes 
FeaturePlot(cortico, features = c("Pomc"), pt.size = 5, label=TRUE)& scale_color_viridis_c() 

#examine DE genes in the  corticotroph cluster
#find the markers for EVERY CLUSTER when compared to remaining cells, only +ve ones
cortico.markers <-FindAllMarkers(cortico, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#examine top 2 deferentially expressed genes in each corticotroph cluster
cortico.markers %>% group_by(cluster) %>% top_n(n=2, wt= avg_logFC)





