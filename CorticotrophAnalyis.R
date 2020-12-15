####ReadME
#Script to study heterogenity, differential gene expression and trajectories in corticotrophs


#############Housekeeping
#First clear the environment;
rm(list=ls())
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
data <- readRDS("NamedClusters.rds")

#visualise the clusters
DimPlot(data, label=TRUE) + NoLegend()


########Isolate coticotrophs
cort <- subset(data, idents = "Corticotroph")

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


####################Differential Expression between clusters

######examine DE genes in the   clusters
#find the markers for EVERY CLUSTER when compared to remaining cells, only +ve ones, 2 fold difference (logfc = 1)
markers <-FindAllMarkers(cort, only.pos = TRUE, min.pct = 0.25, logfc.threshold =1)
#examine top 2 deferentially expressed genes in each corticotroph cluster
df <- markers %>% group_by(cluster) %>% top_n(n=10, wt= avg_logFC)
#plot top 2 of each cluster
head(df, n=30)

#Write tables
write.csv(df,"CorticotrophDifferentialyExpressedGenes.csv")

#################GO terms for differentially expressed genes

#pull out the genes
de.genes <- df$gene
assayed.genes <- rownames(cort)


#construct vector for goseq 
genes=as.integer(assayed.genes%in%de.genes)
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


########Make sense of GO analysis results

#identify categories significantly enriched/unenriched below p=value cut-off

#overenriched using 0.5 FDR (Benjamini and Hochberg 1995)
enriched.GO=GO.wall$category[p.adjust(GO.wall$over_represented_pvalue,method="BH")<.05]
head(enriched.GO)

#save the information to a txt file
sink(file="Corticotroph cluster GO terms.txt")
#Pull information on thes GO terms
for(go in enriched.GO[1:10]){
  print(GOTERM[[go]])
  cat("--------------------------------------\n")
}
sink(file=NULL)

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
png("Corticotroph trajectory anlalysis")

#make the plot
show.velocity.on.embedding.cor(emb, rvel.cd, n=50, scale='sqrt', cex=.8,
                               arrow.scale=0.75, show.grid.flow=TRUE, min.grid.cell.mass=0.5,
                               grid.n=40, arrow.lwd=1,
                               do.par=F, cell.colors = colors,n.cores=4)
#save plot
dev.off()
