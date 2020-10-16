# scRNAseq data analysis

############## precursor, set the directories to save to and read data from:

#set the working directory to directory wish to save analysis results to
setwd("/Users/craigmcdougall/Documents/Project 2/Data analysis/Fletcher")

#check working directory
getwd()


#define the data directory where the data is held
data.dir <- "/Users/craigmcdougall/Documents/Project 2/Data/Fletcher"

#list all the files in the data directory: barcodes, genes, matrix
list.files(data.dir)
#have now seen that the files are there

############################# 1. Load libraries
#libraries required
library(dplyr)
library(Seurat)
library(patchwork)

############################ 2. Load the data  

#load raw data using Read10x 
raw.data <- Read10X(data.dir=data.dir)

#initialise a Seurat object with the raw (not normalised) data. This is a count matrix
data <- CreateSeuratObject(counts=raw.data, project = "data", min.cells = 3, min.features = 200) 

######################### 3. Pre-processing quality control
#selection & filtration of cells based on 
  #QC metrics
  #data normalisation
  #scaling
  #detection of highly desireable features

#QC a& selecting cells for analysis
  #- # genes per cell (id low quality cells or multiple reads)
  #- total # molecules in a cell
  #- % of reads that map to mitochondrial genome - calculatemitochondrial QC using PercentagFeatureSet function, genes starting MT as mitochondrial genes

#note Seurat recomments MT- prefix for mitochondria data
#cheung, ho, mayran use mt
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^mt")
#fletcher uses Mt-
#data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^Mt-")

# the number of unique genes and total molecules are calculated durign CreateSeuratObjkect and can be seen in the metadat

#use QC metrics to filter the cells
#visualise the QC metrics - use a violin plot 
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(data, features = c("nFeature_RNA"), ncol = 1)
#make the plots
#plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#then plot
#plot1 + plot2

############################ 4. filter the cells- based on the QC visualisations
#remove unwanted cells from the data set:
#using the QC visualisations, filter the cells
  #filter cells with unique feature counts of a certain range eg > 2,500 or <200 based on plot to exclude bad cells or multiple reads
  #filter cells >5% MT

#use subset function to do this
data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA <4000 & percent.mt < 5)

########################### 5. Normalising the data
#previous step removed unwanted cells from the dataset
#next step is to normalise the dataset  -default global -scaling "LogNormalize" . Normalise values stores in pbmc[["RNA]]@data

data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)

########################## 6. ID identify highly variable features -high cell-cell variation
  #looking for genes which are highly expressed in come cells but lowly in others
  #community has found these often highlight biological relevant signals

#use the FindVariableFeature function

#find the features
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)

#now identify the 10 most highly variable genes
top10 <- head(VariableFeatures(data),10)

#now plot the variable features with and with out lables
#create plots
plot1 <- VariableFeaturePlot(data)
plot2 <- LabelPoints(plot=plot1, points=top10, repel = TRUE)

#plot them
plot1+plot2

###################### 7. Scaling the data
#standard procedure before dimension reduction techniques
#apply a linear transformation (scaling) 
#use the ScaleData function

#shifts the expression of each gene so mean expression across all cells is 0
#scales expression of each gene so that variance across cells is 1 to give equal weight in downstream analysis and prevents heavily expressed genes from dominating
#result stored in data[["RNA]]@scale.data

all.genes <- rownames(data)
data <- ScaleData(data, features = all.genes)

##################### 8. linear dimensional reduction -----PCA on scaled data
#use the previously determined variable features (by default)

data <- RunPCA(data, features = VariableFeatures(object = data))

#for PCA visualisation, VizDimReduction, DimPlot, DimHeatmap

#can look a the PCA
#print(data[["pca"]], dims = 1:5, nfeatures = 5)

#can visualise the top genes associated with the  reduction components
VizDimLoadings(data, dims=1:2, reduction = "pca")

#then a dimensional reduction plot, 2d scatter plot, each point is a cell an position based on cell embedding determined by PCA 
#DimPlot(data, reduction="pca")

#can also DimHeatmap- easy visualisation of primary sources of heterogeneity.
#DimHeatmap(m_8wk, dims=1, cells = 500, balanced=TRUE)
#dims 1 is only first pc

#DimHeatmap(m_8wk, dims=1:15, cells = 500, balanced=TRUE)
#dims 1:15 is pcs 1-15

#################### 9. determine the dimensionality of the dataset

#seurat clusters based on PCA score - overcome background noise
#deciding how many PC to keep and which to drop....use Jackstraw procedure- takes a long time for big datasets!
#data <- JackStraw(data, num.replicate = 100)#this takes time!
#data <- ScoreJackStraw(data, dims = 1:20)

#visualise to compare the distribution of p-values for each PC
#significant PCs solid curve above the dashed line)
#JackStrawPlot(data, dims = 1:20)

#alternate quick and dirty approach is an elbow plot
ElbowPlot(data)
#look for elbow, majoritory of signal is captured in the PCs before elbow.

###################### 10. clustering the cells
#graph based clustering techniques
#driven by distance metric - based on PCs
#K-nearest neighbour, then partition to communities

#first find neighbours
data <- FindNeighbors(data, dims = 1:20)

#then sets the granularity of the clustering
data <- FindClusters(data, resolution = 0.5)

#can look at the cluster IDs for the cells using idents, eg for first 5 cells
#head(Idents(data),5)

#####################  11.Run non-linear dimensional reduction
#Visualise and explore using UMAP or t-SNE: place cells together in similar low dimension space 
#use the same PCs as clustering analysis

data <- RunUMAP(data, dims = 1:20)

#visualise the clusters
DimPlot(data, reduction = "umap")

#possible to save here so that you dont have to run all the heavy computing again!
#saveRDS(data, file = "data.rds")



#######################   12. Finding deferentially expressed features (biomarkers for each)

#find the markers that define clusters via differential expression 
#what makes each cluster different?

#use findAllMarkers, 25% de threshold, only +ve ones
data.markers <-FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#can look to see top 2 DE genes in each cluster
#increase display optiions to show all
options(dplyr.print_max=1e9)
data.markers %>% group_by(cluster) %>% top_n(n=2, wt= avg_logFC)

#######################   13. Visualising the marker expression
#FEATURE PLOT######### looking for cells differentially expressing POMC, but not Pax7 or Pcsk2
#FeaturePlot(data, features = c("Pomc","Pax7", "Pcsk2" ),label = TRUE)

#add extra makers identified in Cheung et al table 2, most significant cort markers, Crh1, Gpc5, Tbx19, Tnnt1, Tnni3, Gm15543, Cplx3, Erg4, Cdh8, Galnt9, Pomc
# FeaturePlot(data, features = c("Crhr1", "Gpc5","Tbx19","Tnnt1","Tnni3","Gm15543","Cplx3", "Egr4", "Cdh8", "Galnt9", "Pomc", "Pax7", "Pcsk2" ),label = TRUE)
FeaturePlot(data, features = c("Crhr1", "Gpc5","Tbx19","Tnnt1"),label = TRUE)
FeaturePlot(data, features = c("Tnni3","Gm15543","Cplx3", "Egr4" ),label = TRUE)
FeaturePlot(data, features = c("Cdh8", "Galnt9", "Pomc", "Pax7", "Pcsk2" ),label = TRUE)

FeaturePlot(data, features = c("Pcsk2"),label = TRUE)

 FeaturePlot(data, features = c( "Pomc", "Pax7", "Pcsk2" ),label = TRUE)
#select corticotroph cluster based on the feature plot.



#can also generate a heatmap for all cells and features, for top 10 features for each cluster
#top10 <- pbmc.markers %>% group_by(cluster)%>%top_n(n=10, wt=avg_logFC)
#DoHeatmap(pbmc, features = top10$gene)+NoLegend()

# can put names on the featureplot 
#enter the cotricoph cell ids, found in visualisation. Note change to correct cluster #s or will missmatch
new.cluster.ids <- c("0", "1", "2", "3", "4", "cort", "6","7")
 
names(new.cluster.ids) <- levels(data)
 
data <- RenameIdents(data, new.cluster.ids)
#then put on the plot to confirm slected the correct one (sanity check)
DimPlot(data, reduction ="umap", label = TRUE, pt.size = 0.5)+NoLegend()
FeaturePlot(data, features = c("Pomc","Pax7", "Pcsk2" ),label = TRUE)

####################### 14. pull out the corticotrophs & recluster
#extract the expression matrix for cort only
#cort.raw.data <- as.matrix(GetAssayData(data, slot = "counts")[, WhichCells(data, ident = "cort")])

#or pull out a  seurat of cort only
cortico <- subset(data, idents = "cort")


#the data has already ben through QC, scaled and normalised

########################## 15. ID identify highly variable features -high cell-cell variation
#looking for genes which are highly expressed in come cells but lowly in others

#find the features
cortico <- FindVariableFeatures(cortico, selection.method = "vst", nfeatures = 2000)

########################## 16. linear dimensional reduction -----PCA on scaled data
#run PCA on corticotrophs only
cortico <- RunPCA(cortico, features = VariableFeatures(object = cortico))

#can visualise the top genes associated with the  reduction components
#VizDimLoadings(cortico, dims=1:2, reduction = "pca")

#can produce a dimensional reduction plot, 2d scatter plot, each point is a cell an postion based on cel embedding determined by pCA 
#DimPlot(cortico, reduction="pca")


######################### 17. determine the dimensionality of the dataset

#seurat cluseters based on PCA score - overcome backround noise
#deciding how many PC to keep and which to drop....
#use Jackstraw procedure- takes a long time for big datasets!
#ortico <- JackStraw(cortico, num.replicate = 100)#this takes time!
#cortico <- ScoreJackStraw(cortico, dims = 1:20)

#visualise to compare the distribution of p-values for each PC
#significant PCs solid curve above the dashed line)
#JackStrawPlot(cortico, dims = 1:20)

#can see p values drop off after pc 10-12
#alternate quick and dirty approach is an elbow plot
ElbowPlot(cortico)
#elbow seen ~ pc9/10, suggests majoritory of signal is captured in the first 10 PCs.



######################### 18. clustering the cells
#first find neighbours
cortico <- FindNeighbors(cortico, dims = 1:20)

#then find the clustersresolution sets the granularity of the clustering
cortico <- FindClusters(cortico, resolution = 0.5)

#######################  19.Run non-linear dimensional reducuction
#Visualise and explore using UMAP
#place cells together in similar low dimension spae 
#use the same PCs as clustering analysis

cortico <- RunUMAP(cortico, dims = 1:20)

#now need to visualise the clusters, note can use label=TRU to lable the clusters too

DimPlot(cortico, reduction = "umap")
#see 14 different clusters.....cheung did t-sne

#it is possible to save here so that you dont have to run all the heavy computing again!
#saveRDS(m_8wk, file = "Cheung_8wk_m.rds")



####################### 20. Finding differentially expressed features (cluster biomarkers)

#find the markers that define clusters via differential expression 
#what makes each cluster different?

#ident.1 identifies + & _ve markers of a cluster
#min.pct is % threshold of differential expression between the two groups

#use findallmarkers, in cluster 1, 25% de threshold 
#cluster1.markers <- FindMarkers(m_8wk, ident.1 = 1,min.pct = 0.25 )
#then look at first 5 markers
#head(cluster1.markers, n=5)

#find the markers for EVERY CLUSTER when compared to remaing cells, only +ve ones
cortico.markers <-FindAllMarkers(cortico, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


#look to see top 2 differentially expressed genes in each cluster
cortico.markers %>% group_by(cluster) %>% top_n(n=2, wt= avg_logFC)

#test.use sets test for differential expression. ROC test gives classification power
#cluster1.markers <- FindMarkers(m_8wk, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

#######################   21. Visualising the marker expression
#FEATURE PLOT######### from nico presentation, Gh, Prl, Pomc, Pax7, Lhb, Tshb
#FeaturePlot(m_8wk, features = c("Gh", "Prl","Pomc","Pax7","Lhb","Tshb", "Pcsk2" ),label = TRUE)
FeaturePlot(cortico, features = c("Pomc"),label = TRUE)
#FeaturePlot(cortico, features = c("Pomc" ),label = TRUE,cols = c( "darkblue","blue", "green", "yellow"))

