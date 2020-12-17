####ReadME
#Script to study heterogenity, differential gene expression and trajectories in corticotrophs


#############Housekeeping
#First clear the environment;
rm(list=ls())
#clear all plots
dev.off(dev.list()["RStudioGD"])
#Then the cache
gc()

##########Load libraries
library(dplyr) 
library(Seurat)
library(patchwork)
library(sctransform)
library(ggplot2)
library(goseq)
library(GO.db)
library(velocyto.R) 
library(scales)

##########Set seed and working directories
#set seed so plots same each time
set.seed(1234)

#set the working directory to directory wish to read data from and save analysis results to
data.dir <- setwd(here::here())

#####Load in the cluster information from the .rds file
raw<- readRDS("Clusters.rds")


data <- readRDS("NamedClusters.rds")

#visualise the raw clusters
DimPlot(raw, label=TRUE, cols = c("azure3", "aquamarine" ,"blue","blueviolet" , "brown",  "burlywood", "brown1" , "cadetblue","chocolate1", "cornflowerblue","chartreuse", "darkcyan","darkmagenta", "darkslategray", "cyan","darkorange","darkred", "violetred" ,"blue2",  "lawngreen", "magenta", "maroon", "pink")) + NoLegend()+ ggtitle("Pituitary Clusters")
p<-DimPlot(raw, label=TRUE, cols = c("azure3", "aquamarine" ,"blue","blueviolet" , "brown",  "burlywood", "brown1" , "cadetblue","chocolate1", "cornflowerblue","chartreuse", "darkcyan","darkmagenta", "darkslategray", "cyan","darkorange","darkred", "violetred" ,"blue2",  "lawngreen", "magenta", "maroon", "pink")) + NoLegend()+ ggtitle("Pituitary Clusters")
ggsave(plot = p, width = 10, height = 10, dpi = 300, filename = "Pituitary Clusters.png")


#visualise the named clusters
DimPlot(data, label=TRUE) + NoLegend()
p<-DimPlot(data, label=TRUE) + NoLegend()
ggsave(plot = p, width = 10, height = 10, dpi = 300, filename = "Named Pituitary Clusters.png")

########Isolate coticotrophs
cort <- subset(data, idents = "Corticotroph")
DimPlot(cort, label=TRUE) + NoLegend()
#######Re-cluster corticotrophs on their own

#identify principal components
cort <- RunPCA(cort,verbose= FALSE) 
ElbowPlot(cort, ndims=30)

#dimensional reduction 
#UMAP dimensional reduction of the identified PCs, change dims!
cort <- RunUMAP(cort, dims=1:30, verbose=FALSE) #note because SCTransform passes 3000 features can run more PCs

#Find nearest neighbours, change dims!
cort <- FindNeighbors(cort, dims=1:30, verbose=FALSE)

#CLuster based on neighbour distances
cort <- FindClusters(cort, verbose=FALSE, resolution = 0.5)

############Visualise the corticotroph sub clusters
#plot the clusters
DimPlot(cort, pt.size = 5, label=TRUE) + NoLegend()
p<-DimPlot(cort, pt.size = 5, label=TRUE) + NoLegend()
ggsave(plot = p, width = 10, height = 10, dpi = 300, filename = "corticotroph sub-clusters.png")

#########Export the cluster information
saveRDS(cort, file = "../scRNAseq_pipeline/CorticotrophSubClusters.rds")

####################Differential Expression between clusters

######examine DE genes in the clusters
#find the markers for EVERY CLUSTER when compared to remaining cells, only +ve ones.
markers <-FindAllMarkers(cort, only.pos = TRUE, min.pct = 0.25, logfc.threshold =0.25)
#examine top 10 deferentially expressed genes in each corticotroph cluster
de <- markers %>% group_by(cluster) %>% top_n(n=10, wt= avg_logFC)
#plot top 10 of each cluster
head(de, n=30)

#Write tables
write.csv(de,"CorticotrophDifferentialyExpressedGenes.csv")


######examine 2 fold DE genes in the clusters
#find the markers for EVERY CLUSTER when compared to remaining cells, only +ve ones, 2 fold difference (logfc = 1)
markers <-FindAllMarkers(cort, only.pos = TRUE, min.pct = 0.25, logfc.threshold =1)
#examine top 10 deferentially expressed genes in each corticotroph cluster
df <- markers %>% group_by(cluster) %>% top_n(n=10, wt= avg_logFC)
#plot top 10 of each cluster
head(df, n=30)

#Write tables
write.csv(df,"Corticotroph2foldDifferentialyExpressedGenes.csv")


#################GO terms for differentially expressed genes
#pull out the DE genes
de.genes <- de$gene
assayed.genes <- rownames(cort)


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
assayed.genes <- rownames(cort)


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

assayed.genes <- rownames(cort)


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
DoHeatmap(subset(cort, downsample = 100), features = de.genes, size = 3) + ggtitle("Corticotroph subclusters; top 10 DE genes " )
p <- DoHeatmap(subset(cort, downsample = 100), features = de.genes, size = 3)+ ggtitle("Corticotroph subclusters; top 10 DE genes " )
ggsave(plot = p, width = 10, height = 10, dpi = 300, filename = "DE heatmap.png")

#Heatmap for the 2 fold de genes
DoHeatmap(subset(cort, downsample = 100), features = de.genes1, size = 3) + ggtitle("Corticotroph subclusters; 2 fold DE genes " )
p <- DoHeatmap(subset(cort, downsample = 100), features = de.genes, size = 3)+ ggtitle("Corticotroph subclusters; 2 fold DE genes " )
ggsave(plot = p, width = 10, height = 10, dpi = 300, filename = "2fold heatmap.png")


#Heatmap for the 2 fold de genes in cluster 1
DoHeatmap(subset(cort, downsample = 100), features = de.genes2, size = 3) + ggtitle("Corticotroph subcluster 1; 2 fold DE genes " )
p <- DoHeatmap(subset(cort, downsample = 100), features = de.genes2, size = 3) + ggtitle("Corticotroph subcluster 1; 2 fold DE genes " )
ggsave(plot = p, width = 10, height = 10, dpi = 300, filename = "cluster 1 2fold heatmap.png")

#Dotplots
#dotplot for the DE genes
DotPlot(cort, features = de.genes) + RotatedAxis()+ ggtitle("Corticotroph subclusters; top 10 DE genes " )
p <- DotPlot(cort, features = de.genes) + RotatedAxis()+ ggtitle("Corticotroph subclusters; top 10 DE genes " )
ggsave(plot = p, width = 10, height = 10, dpi = 300, filename = "Corticotroph subclusters; top 10 DE genes dot plot.png")

#dotplot for the 2 fold DE genes
DotPlot(cort, features = de.genes1) + RotatedAxis()+ ggtitle("Corticotroph subclusters; 2 fold DE genes " )
p <- DotPlot(cort, features = de.genes1) + RotatedAxis()+ ggtitle("Corticotroph subclusters; 2 fold DE genes " )
ggsave(plot = p, width = 10, height = 10, dpi = 300, filename = "Corticotroph subclusters; 2 fold DE genes dot plot.png")

#dotplot for the 2 fold DE genes in cluster 1
DotPlot(cort, features = de.genes2) + RotatedAxis()+ ggtitle("Corticotroph subcluster 1; 2 fold DE genes " )
p <- DotPlot(cort, features = de.genes1) + RotatedAxis()+ ggtitle("Corticotroph subcluster 1; 2 fold DE genes " )
ggsave(plot = p, width = 10, height = 10, dpi = 300, filename = "Corticotroph subcluster 1; 2 fold DE genes dot plot.png")

########################Perform trajectory analysis

###########Read data

# Load velocyto results from cluster
spliced <- readRDS("VelocytoSpliced.rds")
unspliced <- readRDS("VelocytoUnspliced.rds")

# Match cell names
# Velocyto cell names are saved as something like
# pool1_possorted_genome_bam_HL3XD:ACCGTAAAGTGTTTGCx
# We just want the cell tags in the matrices
colnames(spliced) <- substr(colnames(spliced), 34, 49)
colnames(unspliced) <- substr(colnames(unspliced), 34, 49)

# Velocyto was run only on pool1, so we add -1 at the end
#need to change this in the seurat object
colnames(spliced) <- paste0(colnames(spliced), "-1")
colnames(unspliced) <- paste0(colnames(unspliced), "-1")

# Find common tags between the dataset and velocity results
cell.ids <- intersect(colnames(spliced), Cells(cort))

# Remove data from the other cells
spliced <- spliced[, cell.ids]
unspliced <- unspliced[, cell.ids]

cort<- subset.data.frame(cort, cells = cell.ids)

# Sanity checks - should all return TRUE
all(colnames(spliced) %in% Cells(cort))
all(colnames(unspliced) %in% Cells(cort))
all(Cells(cort) %in% colnames(cort))
all(Cells(cort) %in% colnames(cort))

# Calculate gene correlation
cell.dist <- as.dist(1-armaCor(t(cort[['pca']]@cell.embeddings)))

# Filter genes before velocity calculation
spliced <- filter.genes.by.cluster.expression(spliced, Idents(cort),
                                              min.max.cluster.average = 0.2)
unspliced <- filter.genes.by.cluster.expression(unspliced, Idents(cort), 
                                                min.max.cluster.average = 0.2)

fit.quantile <- 0.02
# Calculate velocity estimates
rvel.cd <- gene.relative.velocity.estimates(spliced, unspliced, deltaT=1, kCells=25,
                                            cell.dist=cell.dist, 
                                            fit.quantile=fit.quantile)

# Get the UMAP coordinates
emb <- cort[["umap"]]@cell.embeddings

# hue_pal returns the ggplot default palette colors
# We generate a named vector with a color for each cell
colors <- hue_pal()(3)[Idents(cort)]
names(colors) = Cells(cort)

#open the plot
png("Corticotroph trajectory anlalysis.png")

#make the plot
p <- show.velocity.on.embedding.cor(emb, rvel.cd, n=50, scale='sqrt', cex=.8,
                               arrow.scale=0.75, show.grid.flow=TRUE, min.grid.cell.mass=0.5,
                               grid.n=40, arrow.lwd=1,
                               do.par=F, cell.colors = colors,n.cores=4)
#save plot
dev.off()

#view the plot
how.velocity.on.embedding.cor(emb, rvel.cd, n=50, scale='sqrt', cex=.8,
                              arrow.scale=0.75, show.grid.flow=TRUE, min.grid.cell.mass=0.5,
                              grid.n=40, arrow.lwd=1,
                              do.par=F, cell.colors = colors,n.cores=4)

