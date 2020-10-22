

## Step 1. Call Libraries & Set Directories

library(dplyr) 
library(Seurat)
library(patchwork)
library(sctransform)
library(ggplot2)
library(scales)



#set the working directory to directory wish to read data from and save analysis results to
data.dir <- setwd(here::here())
#check working directory
getwd()

#list all the files in the data directory, to make sur this is where we expect to be
list.files(data.dir)
#have now seen that the files are there and ok to proceed

#load raw data using Read10x, this will be relative to your path, name desired folder e.g.  "Cheung Data"
raw.data <- Read10X("Cheung Data")

#initialise a Seurat object with the raw (not normalised) data. This is a count matrix
data <- CreateSeuratObject(counts=raw.data, project = "data", min.cells = 3, min.features = 200) 


## Step 2. Quality control plots
#Seurat recommends percent.MT label for mitochondria data, however;
#Cheung, Ho & Mayran et.al. use "percent.mt"
data <- PercentageFeatureSet(data, pattern = "^mt-", col.name = "percent.mt")
#Fletcher uses "percent.Mt"
#data <- PercentageFeatureSet(data, pattern = "^mt-", col.name = "percent.MT")
#select appropriate method and comment out other option

#use QC metrics to filter the cells
#visualise the QC metrics - use a violin plot 
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)

#make the plots
plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#then plot
plot1 + plot2


#The QC plots to inform which cells should be excluded on quality basis (low quality, double counts). In the case of the Cheung data, all cells below 200 counts are dropped and cells with >5% mitochondrial DNA are also dropped.


data <- subset(data, subset = nFeature_RNA > 200  & percent.mt < 5)


## Step 3.  Transform the data


data <- SCTransform(data,vars.to.regress = "percent.mt")#, verbose = FALSE )


## Step 4. Dimensional reduction


data <- RunPCA(data,verbose= FALSE) 



#Do we need to keep all PCs, or is majority of data captured by a certain point?
ElbowPlot(data, ndims=30)




#UMAP dimensional reduction of the identified PCs, only first 20 PCs,   #change dims!
data <- RunUMAP(data, dims=1:20, umap.method = "umap-learn", metric = "correlation", verbose=FALSE)




## Step 5. Clustering


#Find nearest neighbours     #change dims!
data <- FindNeighbors(data, dims=1:20, verbose=FALSE)

#CLuster based on neighbour distances
data <- FindClusters(data, verbose=FALSE)


#plot the clusters
DimPlot(data, label=TRUE) + NoLegend()

## Step 6. Cell Type assessment

#Visualise the canonical marker genes 

#first exclude Pou1f1 expressing cells (Thyrotrophs, Somatotrophs, Lactotrophs)
FeaturePlot(data, features = c("Pou1f1","Tshb", "Gh","Prl"), pt.size = 0.2)&scale_color_viridis_c()

#then exclude the gonadotrophs
FeaturePlot(data, features = c("Lhb"), pt.size = 0.2)&scale_color_viridis_c()

#then exclude the melanotrophs
FeaturePlot(data, features = c("Pomc","Pax7","Pcsk2"), pt.size = 0.2)&scale_color_viridis_c() 

#then id the corticotrophs
FeaturePlot(data, features = c("Pomc", "Crhr1","Avpr1b","Gpc5"), pt.size = 0.2)&scale_color_viridis_c() 


#These plots would suggest cluster 12 is the most likely corticotroph cluster.



## Step 7. Isolate desired cell cluster


#enter the cotricoph cell ids, found in visualisation. Note change to correct cluster #s or will miss-match
#new.cluster.ids <- c("0", "1", "2", "3", "4", "5", "6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22")
new.cluster.ids <- c("0", "1", "2", "3", "4", "5", "6","7","8","9","10","11","cort","13","14","15","16","17","18")

names(new.cluster.ids) <- levels(data)
#Rename the cluster
data <- RenameIdents(data, new.cluster.ids)
#then put on the plot to confirm slected the correct one (sanity check)
DimPlot(data, reduction ="umap", label = TRUE, pt.size = 0.5)+NoLegend()


#differential analysis on corticotroph cluster
#find all markers for cort cluster that make it different from all other clusters
cort.markers <- FindMarkers(data, ident.1="cort", min.pct=0.25)

#display top 10 markers differentially expressed by the corticotroph cluster
head(cort.markers, n=10)



#finally pulll out the corticotroph cluster
cortico <- subset(data, idents = "cort")


## Step 8. Re-cluster corticotroph cluster and study cluster homogeneity

#perform process again on the corticotroph cluster
#dimensional analysis
#PCA to identify PCs
cortico <- RunPCA(cortico,verbose= FALSE) 


#Do we need to keep all PCs, or is majority of data captured by a certain point?
ElbowPlot(cortico, ndims=30)


#dimensional reduction
#UMAP dimensional reduction of the identified PCs, change dims!
cortico <- RunUMAP(cortico, dims=1:12, verbose=FALSE) #note because SCTransform passes 3000 features can run more PCs

#Find nearest neighbours, change dims!
cortico <- FindNeighbors(cortico, dims=1:12, verbose=FALSE)

#CLuster based on neighbour distances
cortico <- FindClusters(cortico, verbose=FALSE)



#A dimensional reduction plot can then reveal any heterogeneity in the cluster;
#plot the clusters
DimPlot(cortico, pt.size = 5, label=TRUE) + NoLegend()


#A featureplot can be used to examine gene expression levels across these clusters within the corticotroph data;
#Visualise the canonical marker genes 
FeaturePlot(cortico, features = c("Pomc"), pt.size = 5, label=TRUE)& scale_color_viridis_c()


#examine DE genes in the  corticotroph cluster
#find the markers for EVERY CLUSTER when compared to remaining cells, only +ve ones
cortico.markers <-FindAllMarkers(cortico, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#examine top 2 deferentially expressed genes in each corticotroph cluster
cortico.markers %>% group_by(cluster) %>% top_n(n=2, wt= avg_logFC)








