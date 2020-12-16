####ReadME
#Script prompts user for a file, then reads the scRNAseq counts in this file, performs qc, transforms and attributes cell type based on thresholds ofor gene expression.


#############Housekeeping
#First clear the environment;
rm(list=ls())
#Then the cache
gc()

###########Call libraries
library(svDialogs)
library(Seurat)
library(dplyr) 
library(patchwork)
library(sctransform)
library(ggplot2)

##########Set working directories

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

###############Transform the data
#normalise, scale & regress out mt
data <- SCTransform(data,vars.to.regress = "percent.mt")#, verbose = FALSE )



###################Thresholding
#Can the cells be identified using thresholds of gene expression?

#Extract the rows for genes of interest, these are canonical genes to aid cell-type identification;
#Create a Seurat object with just the rows with the labels of interest for all columns (all cells) and use this to retrieve the expression matrix for these genes.
canonical <- data[c('Gh', 'Prl','Cga','Tshb','Pomc','Pcsk2','Pax7','Crhr1'), ] 
#retrieve the expression matrix for these markers
canonical <-GetAssayData(canonical)


###########prepare data
#require data in a dataframe, need to convert 'markers' from a sparse matrix to a matrix
#must then transpose the matrix (rows become columns)
#convert to a data frame (for ggplot)
canonical.tr <- data.frame(t(as.matrix(canonical)))
#how many cells
totalCells <- length(rownames(canonical.tr))

##########threshold data to pull out cell types
#plot histogram to identify threshold for each cell type

#begin with hormone producing cell types:

#Somatotrophs
#plot histograms of the canonical marker Gh.
ggplot(data=canonical.tr, aes(Gh))+
  geom_histogram(binwidth = 0.1)
#create an object of which cells to keep based on canonical gene expression levels. 
Somato <- rownames(canonical.tr[canonical.tr$Gh > 7,])
#pull these cells out of the large seurat
Somatotrophs <- subset(data, cells=Somato)
#remove the somatotroph cells from data so cannot be selected twice
data1 <- data[,!colnames(data) %in% colnames(Somatotrophs)]
#check these are split properly
test <- colnames(Somatotrophs) %in% colnames(data1)
table(test)#returns how many true and how many false
#test done. Cells removed. Write to data
data<-data1
#rebuild matrix
canonical <- data[c('Gh', 'Prl','Cga','Tshb','Pomc','Pcsk2','Pax7','Crhr1'), ] 
canonical <-GetAssayData(canonical)
canonical.tr <- data.frame(t(as.matrix(canonical)))


#Lactotrophs
ggplot(data=canonical.tr, aes(Prl))+
  geom_histogram(binwidth = 0.1)
#create an object of which cells to keep based on canonical gene expression levels. 
Lacto <- rownames(canonical.tr[canonical.tr$Prl > 5,])
#pull these cells out of the large seurat
Lactotrophs <- subset(data, cells=Lacto)
data1 <- data[,!colnames(data) %in% colnames(Lactotrophs)]
test <- colnames(Lactotrophs) %in% colnames(data1)
table(test)#returns how many true and how many false
#test done write properly
data <- data1
#rebuild matrix
canonical <- data[c('Gh', 'Prl','Cga','Tshb','Pomc','Pcsk2','Pax7','Crhr1'), ] 
canonical <-GetAssayData(canonical)
canonical.tr <- data.frame(t(as.matrix(canonical)))

#Now for remaining secretory cell types

#Gonadotrophs
ggplot(data=canonical.tr, aes(Cga))+
  geom_histogram(binwidth = 0.1)
#create an object of which cells to keep based on canonical gene expression levels. 
Gonado <- rownames(canonical.tr[canonical.tr$Cga >5,])
#subset the dataset to exclude cells to be filtered out
Gonadotrophs <- subset(data, cells=Gonado)
data1 <- data[,!colnames(data) %in% colnames(Gonadotrophs)]
test <- colnames(Gonadotrophs) %in% colnames(data1)
table(test)#returns how many true and how many false
#test done write properly
data <- data1
#rebuild matrix
canonical <- data[c('Gh', 'Prl','Cga','Tshb','Pomc','Pcsk2','Pax7','Crhr1'), ] 
canonical <-GetAssayData(canonical)
canonical.tr <- data.frame(t(as.matrix(canonical)))


#Thyrotrophs
ggplot(data=canonical.tr, aes(Tshb))+
  geom_histogram(binwidth = 0.1)
#create an object of which cells to keep based on canonical gene expression levels. 
Thyro <- rownames(canonical.tr[canonical.tr$Tshb > 0.5,])
#subset the dataset to exclude cells to be filtered out
Thyrotrophs <- subset(data, cells=Thyro)
data1 <- data[,!colnames(data) %in% colnames(Thyrotrophs)]
test <- colnames(Thyrotrophs) %in% colnames(data1)
table(test)#returns how many true and how many false
#test done write properly
data <- data1
#rebuild matrix
canonical <- data[c('Gh', 'Prl','Cga','Tshb','Pomc','Pcsk2','Pax7','Crhr1'), ] 
canonical <-GetAssayData(canonical)
canonical.tr <- data.frame(t(as.matrix(canonical)))


#Melanotrophs
#Will hve highest Pomc expression (3rd peak) & some Pcsk2
ggplot(data=canonical.tr, aes(Pomc))+
  geom_histogram(binwidth = 0.1)
#create an object of which cells to keep based on canonical gene expression levels. 
Melano <- rownames(canonical.tr[canonical.tr$Pomc > 7 & canonical.tr$Pcsk2 >0,])
#subset the dataset to exclude cells to be filtered out
Melanotrophs <- subset(data, cells=Melano)
data1 <- data[,!colnames(data) %in% colnames(Melanotrophs)]
test <- colnames(Melanotrophs) %in% colnames(data1)
table(test)#returns how many true and how many false
#test done write properly
data <- data1
#rebuild matrix
canonical <- data[c('Gh', 'Prl','Cga','Tshb','Pomc','Pcsk2','Pax7','Crhr1'), ] 
canonical <-GetAssayData(canonical)
canonical.tr <- data.frame(t(as.matrix(canonical)))



#Corticotrophs
#Will hve second highest Pomc expression (2nd peak) and some Crhr1
ggplot(data=canonical.tr, aes(Pomc))+
  geom_histogram(binwidth = 0.1)
#create an object of which cells to keep based on canonical gene expression levels. 
Cort <- rownames(canonical.tr[canonical.tr$Pomc >3 & canonical.tr$Pomc <8.5 ,])
#subset the dataset to exclude cells to be filtered out
Corticotrophs <- subset(data, cells=Cort)
data1 <- data[,!colnames(data) %in% colnames(Corticotrophs)]
test <- colnames(Corticotrophs) %in% colnames(data1)
table(test)#returns how many true and how many false
#test done write properly
data <- data1
#rebuild matrix
canonical <- data[c('Gh', 'Prl','Cga','Tshb','Pomc','Pcsk2','Pax7','Crhr1'), ] 
canonical <-GetAssayData(canonical)
canonical.tr <- data.frame(t(as.matrix(canonical)))


##########Cell proportions
#Somatotrophs
numSomatotrophs <- length(colnames(Somatotrophs))
propSom <- (numSomatotrophs/totalCells)*100

#Lactotrophs
numLactotrophs <- length(colnames(Lactotrophs))
propLac <- (numLactotrophs/totalCells)*100

#Gonadotrophs
numGonadotrophs <- length(colnames(Gonadotrophs))
propGon <- (numGonadotrophs/totalCells)*100

#Thyrotrophs
numThyrotrophs <- length(colnames(Thyrotrophs))
propThy <- (numThyrotrophs/totalCells)*100

#Melanotrophs
numMelanotrophs <- length(colnames(Melanotrophs))
propMel <- (numMelanotrophs/totalCells)*100

#Corticotrophs
numCorticotrophs <- length(colnames(Corticotrophs))
propCor <- (numCorticotrophs/totalCells)*100

#write these to a csv
types <- c("Somatotrophs","Lactotrophs","Gonadotrophs","Thyrotrophs","Melanotrophs","Corticotrophs")
proportions <- c(propSom, propLac,propGon,propThy,propMel,propCor)
props <- data.frame(types,proportions)

#write the proportions to csv file
write.csv(props,"ThresholdCellProportions.csv")

##############Compare threshold corticotrophs to the ones obtained from cluster method

#Need to load in and isolate the corticotrophs from the cluster data
#####Load in the cluster information from the .rds file
data <- readRDS("CorticotrophSubClusters.rds")

#visualise the clusters
DimPlot(data, label=TRUE) + NoLegend()

########Isolate coticotrophs
cortFromCluster <- data#subset(data, idents = "Corticotroph")

cortFromThreshold <- Corticotrophs

#how many cells in each?
CompareCells <- table(Cells(cortFromCluster) %in% Cells(cortFromThreshold))
# %in% operator tells you whether each of the elements of what's on the left of it is contained in what is on the right of it.

#It is true that 528 of the cells obtained from the threshold method are found using the cluster method

#percentage cells same using threshold 528/648
percentSameUsingThresh <- (CompareCells[2]/length(colnames(cortFromCluster)))*100
#81.5% of cluster method cells are obtained useing threshold method.

#sanity check, 0.815 * 648 is
sanity <- 0.815*length(colnames(cortFromCluster))
#see 528 cells by both ways

#120 cells are different. Need to obtain their tags

#find idxs of missing cells (all FALSE)
cellIdx<- colnames(cortFromCluster) %in% colnames(cortFromThreshold)
#sanity check
sanity <- table(cellIdx)
#returns 120 FALSE (missing) and 528 matched as expected.

#pull out the indexes for the 120 missing
missingIdxs <- which(cellIdx == FALSE)

# get the cell names that are missing
missingCells <- subset(cortFromCluster, cells= missingIdxs)

sanity <- length(colnames(missingCells))
#has cut down to 120 as expected 120 missing cells

#which clusters are the missing cells from?
clustersmissing <- missingCells$seurat_clusters
table(clustersmissing)


######examine DE genes in the threshold method clusters

#to do this need to cluster the threshold identified corticotrophs
#identify principal components
cortFromThreshold <- RunPCA(cortFromThreshold,verbose= FALSE) 
ElbowPlot(cortFromThreshold, ndims=30)

#dimensional reduction 
#UMAP dimensional reduction of the identified PCs, change dims!
cortFromThreshold <- RunUMAP(cortFromThreshold, dims=1:30, verbose=FALSE) #note because SCTransform passes 3000 features can run more PCs

#Find nearest neighbours, change dims!
cortFromThreshold <- FindNeighbors(cortFromThreshold, dims=1:30, verbose=FALSE)

#CLuster based on neighbour distances
cortFromThreshold <- FindClusters(cortFromThreshold, verbose=FALSE, resolution = 0.5)

############Visualise the corticotroph sub clusters
#plot the clusters
DimPlot(cortFromThreshold, pt.size = 5, label=TRUE) + NoLegend()
p<-DimPlot(cortFromThreshold, pt.size = 5, label=TRUE) + NoLegend()
ggsave(plot = p, width = 10, height = 10, dpi = 300, filename = "Threshold corticotroph sub-clusters.png")

#########Export the cluster information
saveRDS(cortFromThreshold, file = "../scRNAseq_pipeline/THresholdCorticotrophSubClusters.rds")

####################Differential Expression between clusters

######examine DE genes in the   clusters
#find the markers for EVERY CLUSTER when compared to remaining cells, only +ve ones, 2 fold difference (logfc = 1)
markers <-FindAllMarkers(cortFromThreshold, only.pos = TRUE, min.pct = 0.25, logfc.threshold =1)
#examine top 2 deferentially expressed genes in each corticotroph cluster
df <- markers %>% group_by(cluster) %>% top_n(n=10, wt= avg_logFC)
#plot top 2 of each cluster
head(df, n=30)

#Write tables
write.csv(df,"ThresholdCorticotrophDifferentialyExpressedGenes.csv")





