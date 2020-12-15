####ReadME
#Script to enale cell type attirbution to cluster


#############Housekeeping
#First clear the environment;
rm(list=ls())
#Then the cache
gc()

###########Call libraries
library(dplyr) 
library(Seurat)
library(patchwork)
library(ggplot2)

##########Set seed and working directories
#set seed so plots same each time
set.seed(1234)

#set the working directory to directory wish to read data from and save analysis results to
data.dir <- setwd(here::here())

#####Load in the cluster information from the .rds file
data <- readRDS("Clusters.rds")

#open plot of clusters
png("unamed clusters")
#visualise the clusters- make plot
DimPlot(data, label=TRUE) + NoLegend()
#save plot
dev.off
########Visualise expression
#How well do the clusters represent pituitary cell types?
#use canonical markers to attriubute cell type.
#viloin plot of expression to look for most common cell tpes (Somatotrophs, Lactotrophs,Gonadotrophs,Thyrotrophs Melanotrophs, Corticotrophs. )
#se canonical features
canonical <- c("Gh","Prl","Cga","Tshb","Pomc", "Pcsk2")

#different visualisations
VlnPlot(data, features = canonical, ncol = 2, pt.size=0)
RidgePlot(data, features = canonical, ncol = 3)
FeaturePlot(data, features = canonical)
DotPlot(data, features = canonical) + RotatedAxis()
DoHeatmap(subset(data, downsample = 100), features = canonical, size = 3)

######identify clusters
#Use feature plots to confirm and attibute pituitary secretary cells

#Lactotrophs
FeaturePlot(data, features = c("Pou1f1","Prl"),  ncol = 2, pt.size = 1)&scale_color_viridis_c()
#plot<- FeaturePlot(data, features = c("Prl"),  ncol = 2, pt.size = 1)&scale_color_viridis_c()
#select and assign lactotrophs
#data <- CellSelector(plot = plot, object = data, ident = "Lactotrophs")
#object <- HoverLocator(plot = plot, information = FetchData(data, vars = c("ident")))
#nned to check back when lasso function has been enabled as rectangle selection is not proctical

#Somatotrophs
FeaturePlot(data, features = c("Pou1f1","Gh","Ghrhr"), pt.size = 1)&scale_color_viridis_c()

#Gonadotrophs
FeaturePlot(data, features = c("Fshb", "Lhb", "Cga", "Gnrhr"), pt.size = 1)&scale_color_viridis_c()

#Thyrotophs
FeaturePlot(data, features = c("Pou1f1","Tshb"), pt.size = 1)&scale_color_viridis_c()

#Melaotroph
FeaturePlot(data, features = c("Pomc","Pcsk2","Pax7","Tbx19"), pt.size = 1)&scale_color_viridis_c()

#Corticotroph
FeaturePlot(data, features = c("Pomc","Crhr1","Tbx19","Avpr1b"), pt.size = 1)&scale_color_viridis_c()

#Also identify other cell types

#posteriour cells
FeaturePlot(data, features = c("Nkx2-1","Lhx2"), pt.size = 1)&scale_color_viridis_c()

#Mki67+
FeaturePlot(data, features = c("Mki67","Top2a"), pt.size = 1)&scale_color_viridis_c()

#stem cells
FeaturePlot(data, features = c("Sox2","Fstl1"), pt.size = 1)&scale_color_viridis_c()

#collagen (FSC cells)
FeaturePlot(data, features = c("Col1a1","Dcn"), pt.size = 1)&scale_color_viridis_c()

#White blood cells
FeaturePlot(data, features = c("C1qa","Cd53"), pt.size = 1)&scale_color_viridis_c()

#Red blood cells
FeaturePlot(data, features = c("Hbb-bt","Hba-a1"), pt.size = 1)&scale_color_viridis_c()

#endothelial cells
FeaturePlot(data, features = c("Pecam1","Plvap"), pt.size = 1)&scale_color_viridis_c()

############Rename the clusters

#enter new cluster ids
new.cluster.ids <- c("Somatotroph","Somatotroph","Somatotroph","Somatotroph", "Lactotroph", "Lactotroph", "Lactotroph","Gonadotroph", "Lactotroph","Somatotroph","Stem Cells","Endothelial cells","Mki67+","Somatotroph","Corticotroph","Corticotroph","Melanotroph","Somatotroph","Collagen","WBC","Thyrotroph","RBC","posteriour pituitary")

#make new cluster ids names
names(new.cluster.ids) <- levels(data)

#apply new names to the lusters
data <- RenameIdents(data, new.cluster.ids)

#make a plot of named clusters
png("named clusters")
#visualise the named clusters-make plot
DimPlot(data, reduction ="umap", label = TRUE, pt.size = 0.5)+NoLegend()
#saveplot
dev.off

#export the named clusters
saveRDS(data, file = "../scRNAseq_pipeline/NamedClusters.rds")
