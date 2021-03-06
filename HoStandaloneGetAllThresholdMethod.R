####ReadME
#Script prompts user for a file, then reads the scRNAseq counts in this file, performs qc, transforms and attributes cell type based on thresholds ofor gene expression.


#############Housekeeping
#First clear the environment;
rm(list=ls())
#clear all plots
dev.off(dev.list()["RStudioGD"])
#Then the cache
gc()

###########Call libraries
library(svDialogs)
library(Seurat)
library(dplyr) 
library(patchwork)
library(sctransform)
library(ggplot2)
library(goseq)
library(GO.db)
library(velocyto.R) 
library(scales)

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
canonical <- data[c('Pou1f1','Gh','Ghrhr', 'Prl','Fshb','Lhb', 'Gnhrhr' , 'Pomc','Pax7','Tbx19','Crhr1'), ] 
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

#Melanotrophs
#POMC, Pax7,Tbx19
ggplot(data=canonical.tr, aes(Pomc))+
  geom_histogram(binwidth = 0.1)
#create an object of which cells to keep based on canonical gene expression levels. 
#Melano <- rownames(canonical.tr[canonical.tr$Pomc > 7 & canonical.tr$Pcsk2 >0,])
Melano <- rownames(canonical.tr[canonical.tr$Pomc >0 & canonical.tr$Pax7>0 & canonical.tr$Tbx19 >0 ,])
#subset the dataset to exclude cells to be filtered out
Melanotrophs <- subset(data, cells=Melano)
data1 <- data[,!colnames(data) %in% colnames(Melanotrophs)]
test <- colnames(Melanotrophs) %in% colnames(data1)
table(test)#returns how many true and how many false
#test done write properly
data <- data1
#rebuild matrix for Ho, note they found no thyrotrophs.....
canonical <- data[c('Pou1f1','Gh','Ghrhr', 'Prl','Fshb','Lhb', 'Gnhrhr' , 'Pomc','Pax7','Tbx19','Crhr1'), ] 
canonical <-GetAssayData(canonical)
canonical.tr <- data.frame(t(as.matrix(canonical)))



#Corticotrophs
#POMC, Tbx19, Crhr1
ggplot(data=canonical.tr, aes(Pomc))+
  geom_histogram(binwidth = 0.1)
#create an object of which cells to keep based on canonical gene expression levels. 
Cort <- rownames(canonical.tr[canonical.tr$Pomc >0 & canonical.tr$Tbx19 >0 & canonical.tr$Crhr1 > 0 ,])

#subset the dataset to exclude cells to be filtered out
Corticotrophs <- subset(data, cells=Cort)
data1 <- data[,!colnames(data) %in% colnames(Corticotrophs)]
test <- colnames(Corticotrophs) %in% colnames(data1)
table(test)#returns how many true and how many false
#test done write properly
data <- data1
#rebuild matrix
canonical <- data[c('Pou1f1','Gh','Ghrhr', 'Prl','Fshb','Lhb', 'Gnhrhr' , 'Pomc','Pax7','Tbx19','Crhr1'), ] 
canonical <-GetAssayData(canonical)
canonical.tr <- data.frame(t(as.matrix(canonical)))


#Somatotrophs
#Pou1f1, Gh, Ghrhr
ggplot(data=canonical.tr, aes(Gh))+
  geom_histogram(binwidth = 0.1)
#create an object of which cells to keep based on canonical gene expression levels. 
Somato <- rownames(canonical.tr[canonical.tr$Pou1f1 > 0 & canonical.tr$Gh > 0 & canonical.tr$Ghrhr > 0,])
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
canonical <- data[c('Pou1f1','Gh','Ghrhr', 'Prl','Fshb','Lhb', 'Gnhrhr' , 'Pomc','Pax7','Tbx19','Crhr1'), ] 
canonical <-GetAssayData(canonical)
canonical.tr <- data.frame(t(as.matrix(canonical)))


#Lactotrophs
#Pou1f1, Prl
ggplot(data=canonical.tr, aes(Prl))+
  geom_histogram(binwidth = 0.1)
#create an object of which cells to keep based on canonical gene expression levels. 
Lacto <- rownames(canonical.tr[canonical.tr$Pou1f1 > 0 & canonical.tr$Prl > 0,])
#pull these cells out of the large seurat
Lactotrophs <- subset(data, cells=Lacto)
data1 <- data[,!colnames(data) %in% colnames(Lactotrophs)]
test <- colnames(Lactotrophs) %in% colnames(data1)
table(test)#returns how many true and how many false
#test done write properly
data <- data1
#rebuild matrix
canonical <- data[c('Pou1f1','Gh','Ghrhr', 'Prl','Fshb','Lhb', 'Gnhrhr' , 'Pomc','Pax7','Tbx19','Crhr1'), ] 
canonical <-GetAssayData(canonical)
canonical.tr <- data.frame(t(as.matrix(canonical)))


#Gonadotrophs
#Fshb, Lhb, Gnhrhr
ggplot(data=canonical.tr, aes(Fshb))+
  geom_histogram(binwidth = 0.1)
#create an object of which cells to keep based on canonical gene expression levels. 
Gonado <- rownames(canonical.tr[canonical.tr$Fshb >0 & canonical.tr$Lhb >0 & canonical.tr$Gnhrhr >0 ,])
#subset the dataset to exclude cells to be filtered out
Gonadotrophs <- subset(data, cells=Gonado)
data1 <- data[,!colnames(data) %in% colnames(Gonadotrophs)]
test <- colnames(Gonadotrophs) %in% colnames(data1)
table(test)#returns how many true and how many false
#test done write properly
data <- data1
#rebuild matrix
canonical <- data[c('Pou1f1','Gh','Ghrhr', 'Prl','Fshb','Lhb', 'Gnhrhr' , 'Pomc','Pax7','Tbx19','Crhr1'), ] 
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



############Visualise the corticotroph sub clusters
#plot the clusters
DimPlot(Corticotrophs, pt.size = 5, label=TRUE) + NoLegend()
p<-DimPlot(Corticotrophs, pt.size = 5, label=TRUE) + NoLegend()
ggsave(plot = p, width = 10, height = 10, dpi = 300, filename = "Original threshold method corticotroph sub-clusters.png")
#see two unexpected clusters, 3 & 4.
#look at how many cells in clusters
dist <- (Idents(Corticotrophs))
extraClusterCellNum <-which(as.numeric(dist)==4,5)
length(extraClusterCellNum)
#it is 52 cells. 

#does differential expression at this point indicate what the unexpected clusters might be?are they corticotrophs?
#find the markers for EVERY CLUSTER when compared to remaining cells, only +ve ones.
markers <-FindAllMarkers(Corticotrophs, only.pos = TRUE, min.pct = 0.25, logfc.threshold =0.25)
#examine top 2 deferentially expressed genes in each corticotroph cluster
idExtraCells<- markers %>% group_by(cluster) %>% top_n(n=2, wt= avg_logFC)
#plot top 2 of each cluster
head(idExtraCells, n=10)

#can see the unexpected clusters 3 & 4 are a mix of somatotrophs(Gh) lactotrophs(Prl) and GOnadotrophs (Lhb) respectively. 
#remove these clusters from the data.
Corticotrophs <- subset(Corticotrophs, idents =c("0","1","2") )  

#replot to check
DimPlot(Corticotrophs, pt.size = 5, label=TRUE) + NoLegend()
p<-DimPlot(Corticotrophs, pt.size = 5, label=TRUE) + NoLegend()
ggsave(plot = p, width = 10, height = 10, dpi = 300, filename = "Threshold corticotroph sub-clusters.png")

#########Export the cluster information
saveRDS(cortFromThreshold, file = "../scRNAseq_pipeline/THresholdCorticotrophSubClusters.rds")


######################################################################################
#                                 ANALYSIS                                           #
#########################examine differential expression##############################



####################Differential Expression between clusters

######examine DE genes in the clusters
#find the markers for EVERY CLUSTER when compared to remaining cells, only +ve ones.
markers <-FindAllMarkers(Corticotrophs, only.pos = TRUE, min.pct = 0.25, logfc.threshold =0.25)
#examine top 10 deferentially expressed genes in each corticotroph cluster
de <- markers %>% group_by(cluster) %>% top_n(n=10, wt= avg_logFC)
#plot top 10 of each cluster
head(de, n=30)

#Write tables
write.csv(de,"CorticotrophDifferentialyExpressedGenes.csv")


######examine 2 fold DE genes in the clusters
#find the markers for EVERY CLUSTER when compared to remaining cells, only +ve ones, 2 fold difference (logfc = 1)
markers <-FindAllMarkers(Corticotrophs, only.pos = TRUE, min.pct = 0.25, logfc.threshold =1)
#examine top 10 deferentially expressed genes in each corticotroph cluster
df <- markers %>% group_by(cluster) %>% top_n(n=10, wt= avg_logFC)
#plot top 10 of each cluster
head(df, n=30)

#Write tables
write.csv(df,"Corticotroph2foldDifferentialyExpressedGenes.csv")


#################GO terms for differentially expressed genes
#pull out the DE genes
de.genes <- de$gene
assayed.genes <- rownames(Corticotrophs)


#construct vector for goseq 
#which genes are there (0/1)
genes=as.integer(assayed.genes%in%de.genes)
#get the names for the genes
names(genes)=assayed.genes
#get a summary to check our 17 2 fold DE genes are there
table(genes)


#Perform GO analysis

#calculate bias (PWF)
#obtain weighting of each gene depending on its length to remove length bias.
#using gene symbols
pwf=nullp(genes,"mm9","geneSymbol")
head(pwf)

#use Wallenius approximation to calculate over and under expressed GO categories amonsgst DE genes.
#goseq will automatically fetch the data describing relationship between the gene IDs and the GO categories
GO.wall=goseq(pwf,"mm9","geneSymbol")
head(GO.wall)

#################GO terms for 2 fold differentially expressed genes
#pull out the 2 fold DE genes
de.genes1 <- df$gene
assayed.genes <- rownames(Corticotrophs)


#construct vector for goseq 
#which genes are there (0/1)
genes=as.integer(assayed.genes%in%de.genes1)
#get the names for the genes
names(genes)=assayed.genes
#get a summary to check our 17 2 fold DE genes are there
table(genes)


#Perform GO analysis

#calculate bias (PWF)
#obtain weighting of each gene depending on its length to remove length bias.
#using gene symbols
pwf=nullp(genes,"mm9","geneSymbol")
head(pwf)

#use Wallenius approximation to calculate over and under expressed GO categories amonsgst DE genes.
#goseq will automatically fetch the data describing relationship between the gene IDs and the GO categories
GO.wall=goseq(pwf,"mm9","geneSymbol")
head(GO.wall)

#################GO terms for 2fold DE for cluster 1

#pull out the genes for cluster 1
#get indexes for xluster1
cluster1idx <- which(df$cluster == 1)

#pull out the desired genes using index
whichgenes <- df[cluster1idx,]
de.genes2 <- whichgenes$gene

assayed.genes <- rownames(Corticotrophs)


#construct vector for goseq 
genes=as.integer(assayed.genes%in%de.genes2)
names(genes)=assayed.genes
head(genes)

#Perform GO analysis

#calculate bias (PWF)
#obtain weighting of each gene depending on its length to remove length bias.
#using gene symbols
pwf=nullp(genes,"mm9","geneSymbol")
head(pwf)

#use Wallenius approximation to calculate over and under expressed GO categories amonsgst DE genes.
#goseq will automatically fetch the data describing relationship between the gene IDs and the GO categories
GO.wall=goseq(pwf,"mm9","geneSymbol")
head(GO.wall)
########Make sense of GO analysis results for the 2 fold DE genes in cluster One

#identify categories significantly enriched/unenriched below p=value cut-off

#overenriched using 0.5 FDR (Benjamini and Hochberg 1995)
enriched.GO=GO.wall$category[p.adjust(GO.wall$over_represented_pvalue,method="BH")<.05]
head(enriched.GO)

#save the information to a txt file
sink(file="Corticotroph cluster 1 2 fold DE genes GO terms.txt")
#Pull information on thes GO terms
for(go in enriched.GO[1:10]){
  print(GOTERM[[go]])
  cat("--------------------------------------\n")
}
sink(file=NULL)
#view the info
for(go in enriched.GO[1:10]){
  print(GOTERM[[go]])
  cat("--------------------------------------\n")
}
#########################Visualise corticotroph DE genes


#######Heatmaps
#Heatmap for the de genes
DoHeatmap(subset(Corticotrophs, downsample = 100), features = de.genes, size = 3) + ggtitle("Corticotroph subclusters; top 10 DE genes " )
p <- DoHeatmap(subset(Corticotrophs, downsample = 100), features = de.genes, size = 3)+ ggtitle("Corticotroph subclusters; top 10 DE genes " )
ggsave(plot = p, width = 10, height = 10, dpi = 300, filename = "DE heatmap.png")

#Heatmap for the 2 fold de genes
DoHeatmap(subset(Corticotrophs, downsample = 100), features = de.genes1, size = 3) + ggtitle("Corticotroph subclusters; 2 fold DE genes " )
p <- DoHeatmap(subset(Corticotrophs, downsample = 100), features = de.genes, size = 3)+ ggtitle("Corticotroph subclusters; 2 fold DE genes " )
ggsave(plot = p, width = 10, height = 10, dpi = 300, filename = "2fold heatmap.png")


#Heatmap for the 2 fold de genes in cluster 1
DoHeatmap(subset(Corticotrophs, downsample = 100), features = de.genes2, size = 3) + ggtitle("Corticotroph subcluster 1; 2 fold DE genes " )
p <- DoHeatmap(subset(Corticotrophs, downsample = 100), features = de.genes2, size = 3) + ggtitle("Corticotroph subcluster 1; 2 fold DE genes " )
ggsave(plot = p, width = 10, height = 10, dpi = 300, filename = "cluster 1 2fold heatmap.png")

#Dotplots
#dotplot for the DE genes
DotPlot(Corticotrophs, features = de.genes) + RotatedAxis()+ ggtitle("Corticotroph subclusters; top 10 DE genes " )
p <- DotPlot(Corticotrophs, features = de.genes) + RotatedAxis()+ ggtitle("Corticotroph subclusters; top 10 DE genes " )
ggsave(plot = p, width = 10, height = 10, dpi = 300, filename = "Corticotroph subclusters; top 10 DE genes dot plot.png")

#dotplot for the 2 fold DE genes
DotPlot(Corticotrophs, features = de.genes1) + RotatedAxis()+ ggtitle("Corticotroph subclusters; 2 fold DE genes " )
p <- DotPlot(Corticotrophs, features = de.genes1) + RotatedAxis()+ ggtitle("Corticotroph subclusters; 2 fold DE genes " )
ggsave(plot = p, width = 10, height = 10, dpi = 300, filename = "Corticotroph subclusters; 2 fold DE genes dot plot.png")

#dotplot for the 2 fold DE genes in cluster 1
DotPlot(Corticotrophs, features = de.genes2) + RotatedAxis()+ ggtitle("Corticotroph subcluster 1; 2 fold DE genes " )
p <- DotPlot(Corticotrophs, features = de.genes1) + RotatedAxis()+ ggtitle("Corticotroph subcluster 1; 2 fold DE genes " )
ggsave(plot = p, width = 10, height = 10, dpi = 300, filename = "Corticotroph subcluster 1; 2 fold DE genes dot plot.png")


