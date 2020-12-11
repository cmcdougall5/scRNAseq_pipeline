#housekeeping
#First clear the environment;
rm(list=ls())
#Then the cache
gc()

#load libraries
library(velocyto.R) 
library(scales)
library(Seurat)

data.dir <- setwd(here::here())

#so clusters are the same each time
set.seed(1234)

###########Read data

# Load velocyto results
spliced <- readRDS("VelocytoSpliced.rds")
unspliced <- readRDS("VelocytoUnspliced.rds")

# Load the data - Craig, you can load your own data, 
# this is simply the Seurat object that I saved onto a rds file
data <- readRDS("CheungCorticotrophs.rds")

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
cell.ids <- intersect(colnames(spliced), Cells(data))

# Remove data from the other cells
spliced <- spliced[, cell.ids]
unspliced <- unspliced[, cell.ids]

data<- subset.data.frame(data, cells = cell.ids)

# Sanity checks - should all return TRUE
all(colnames(spliced) %in% Cells(data))
all(colnames(unspliced) %in% Cells(data))
all(Cells(data) %in% colnames(data))
all(Cells(data) %in% colnames(data))

# Calculate gene correlation
cell.dist <- as.dist(1-armaCor(t(data[['pca']]@cell.embeddings)))

# Filter genes before velocity calculation
spliced <- filter.genes.by.cluster.expression(spliced, Idents(data),
                                              min.max.cluster.average = 0.2)
unspliced <- filter.genes.by.cluster.expression(unspliced, Idents(data), 
                                              min.max.cluster.average = 0.2)

fit.quantile <- 0.02
# Calculate velocity estimates
rvel.cd <- gene.relative.velocity.estimates(spliced, unspliced, deltaT=1, kCells=25,
                                            cell.dist=cell.dist, 
                                            fit.quantile=fit.quantile)

# Get the UMAP coordinates
emb <- data[["umap"]]@cell.embeddings

# hue_pal returns the ggplot default palette colors
# We generate a named vector with a color for each cell
colors <- hue_pal()(3)[Idents(data)]
names(colors) = Cells(data)

show.velocity.on.embedding.cor(emb, rvel.cd, n=50, scale='sqrt', cex=.8,
                               arrow.scale=0.2, show.grid.flow=TRUE, min.grid.cell.mass=0.5,
                               grid.n=40, arrow.lwd=1,
                               do.par=F, cell.colors = colors,n.cores=4)

