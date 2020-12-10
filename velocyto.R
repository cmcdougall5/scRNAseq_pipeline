library(velocyto.R) 
library(scales)

setwd("/lab/Lab/scRNAseqPit/")

# Load velocyto results
spliced <- readRDS("VelocytoSpliced.rds")
unspliced <- readRDS("VelocytoUnspliced.rds")

# Load the data - Craig, you can load your own data, 
# this is simply the Seurat object that I saved onto a rds file
seurat.Pomc <- readRDS("seuratPOMC.rds")

# Match cell names
# Velocyto cell names are saved as something like
# pool1_possorted_genome_bam_HL3XD:ACCGTAAAGTGTTTGCx
# We just want the cell tags
colnames(spliced) <- substr(colnames(spliced), 34, 49)
colnames(unspliced) <- substr(colnames(unspliced), 34, 49)

# Velocyto was run only on pool1, so we add -1 at the end
colnames(spliced) <- paste0(colnames(spliced), "-1")
colnames(unspliced) <- paste0(colnames(unspliced), "-1")

# Find common tags between the dataset and velocity results
cell.ids <- intersect(colnames(spliced), Cells(seurat.Pomc))

# Remove data from the other cells
spliced <- spliced[, cell.ids]
unspliced <- unspliced[, cell.ids]

seurat.Pomc <- subset(seurat.Pomc, cells = cell.ids)

# Sanity checks - should all return TRUE
all(colnames(spliced) %in% Cells(seurat.Pomc))
all(colnames(unspliced) %in% Cells(seurat.Pomc))
all(Cells(seurat.Pomc) %in% colnames(spliced))
all(Cells(seurat.Pomc) %in% colnames(unspliced))

# Calculate gene correlation
cell.dist <- as.dist(1-armaCor(t(seurat.Pomc[['pca']]@cell.embeddings)))

# Filter genes before velocity calculation
spliced <- filter.genes.by.cluster.expression(spliced, Idents(seurat.Pomc),
                                              min.max.cluster.average = 0.2)
unspliced <- filter.genes.by.cluster.expression(unspliced, Idents(seurat.Pomc), 
                                              min.max.cluster.average = 0.2)

fit.quantile <- 0.02
# Calculate velocity estimates
rvel.cd <- gene.relative.velocity.estimates(spliced, unspliced, deltaT=1, kCells=25,
                                            cell.dist=cell.dist, 
                                            fit.quantile=fit.quantile)

# Get the UMAP coordinates
emb <- seurat.Pomc[["umap"]]@cell.embeddings

# hue_pal returns the ggplot default palette colors
# We generate a named vector with a color for each cell
colors <- hue_pal()(3)[Idents(seurat.Pomc)]
names(colors) = Cells(seurat.Pomc)

show.velocity.on.embedding.cor(emb, rvel.cd, n=50, scale='sqrt', cex=.8,
                               arrow.scale=5, show.grid.flow=TRUE, min.grid.cell.mass=0.5,
                               grid.n=40, arrow.lwd=1,
                               do.par=F, cell.colors = colors)

