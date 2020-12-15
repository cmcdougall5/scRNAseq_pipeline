####ReadME
#Script to calculate the differentially expressed genes of each cluster.


#############Housekeeping
#First clear the environment;
rm(list=ls())
#Then the cache
gc()

##########Load libraries
library(dplyr) 

##########Set seed and working directories
#set seed so plots same each time
set.seed(1234)

#set the working directory to directory wish to read data from and save analysis results to
data.dir <- setwd(here::here())

#####Load in the cluster information from the .rds file
data <- readRDS("NamedClusters.rds")

#visualise the clusters
DimPlot(data, label=TRUE) + NoLegend()

######obtain all genes (for GO analysis)
Genes <- rownames(data)

######examine DE genes in the   clusters
#find the markers for EVERY CLUSTER when compared to remaining cells, only +ve ones, 2 fold difference (logfc = 1)
markers <-FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold =1)
#examine top 2 deferentially expressed genes in each corticotroph cluster
df <- markers %>% group_by(cluster) %>% top_n(n=2, wt= avg_logFC)

#plot top 2 of each cluster
head(df, n=26)

#Write tables
write.csv(df,"DifferentialyExpressedGenes.csv")
write.csv(Genes,"AllGenes.csv")
