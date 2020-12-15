####ReadME
#SScript to calculate the proportions of each cell type.


#############Housekeeping
#First clear the environment;
rm(list=ls())
#Then the cache
gc()

##########Set seed and working directories
#set seed so plots same each time
set.seed(1234)

#set the working directory to directory wish to read data from and save analysis results to
data.dir <- setwd(here::here())

#####Load in the cluster information from the .rds file
data <- readRDS("NamedClusters.rds")

#visualise the clusters
DimPlot(data, label=TRUE) + NoLegend()


###########Evaluate how many cells in each cluster/total

#Find how many cells per cluster
cellscluster <- table(Idents(data))
#How many cells total and per cluster?
cellstotal <- sum(cellscluster)
#make cells per cluster a datafram  
cellspercluster <- data.frame(cellscluster)

############ Calculate cell proportions

#how to do this in for loop or vectorised?

#Lactotrophs
#Find the lactotroph index
lactoidx <- which(cellspercluster=="Lactotroph")
lactos <- cellspercluster[lactoidx,"Freq"]
proplac <- (lactos/cellstotal)*100

#Somatotrophs
somatoidx <- which(cellspercluster=="Somatotroph")
somatos <- cellspercluster[somatoidx,"Freq"]
propsom <- (somatos/cellstotal)*100

#Gonadotrophs
gonadoidx <- which(cellspercluster=="Gonadotroph")
gonados <- cellspercluster[gonadoidx,"Freq"]
propgon <- (gonados/cellstotal)*100

#Thyrotophs
thyroidx <- which(cellspercluster=="Thyrotroph")
thyros <- cellspercluster[thyroidx,"Freq"]
propthy <- (thyros/cellstotal)*100

#Melaotroph
melanoidx <- which(cellspercluster=="Melanotroph")
melanos <- cellspercluster[melanoidx,"Freq"]
propmel <- (melanos/cellstotal)*100

#Corticotroph
cortidx <- which(cellspercluster=="Corticotroph")
corts <- cellspercluster[cortidx,"Freq"]
propcor <- (corts/cellstotal)*100

#posteriour cells
postidx <- which(cellspercluster=="posteriour pituitary")
posts <- cellspercluster[postidx,"Freq"]
proppost <- (posts/cellstotal)*100

#Mki67+
Mki67idx <- which(cellspercluster=="Mki67+")
Mki67s <- cellspercluster[Mki67idx,"Freq"]
propMki67 <- (Mki67s/cellstotal)*100

#stem cells
stemidx <- which(cellspercluster=="Stem Cells")
stems <- cellspercluster[stemidx,"Freq"]
propstem <- (stems/cellstotal)*100

#collagen (FSC cells)
#Find the corticotroph index
colidx <- which(cellspercluster=="Collagen")
cols <- cellspercluster[colidx,"Freq"]
propcol <- (cols/cellstotal)*100

#White blood cells
#Find the corticotroph index
wbcidx <- which(cellspercluster=="WBC")
wbcs <- cellspercluster[wbcidx,"Freq"]
propwbc <- (wbcs/cellstotal)*100

#Red blood cells
rbcidx <- which(cellspercluster=="RBC")
rbcs <- cellspercluster[rbcidx,"Freq"]
proprbc <- (rbcs/cellstotal)*100

#endothelial cells
endoidx <- which(cellspercluster=="Endothelial cells")
endos <- cellspercluster[endoidx,"Freq"]
propendo <- (endos/cellstotal)*100


#####Write the data to file

#pull out info
ids <- c("Lactotrophs","Somatotrophs","Gonadotrophs","Thyrotophs",
        "Melaotroph","Corticotroph","posteriour cells", "Mki67+",
        "Stem Cells","Collagen","WBC","RBC","Endothelial cells")

proportions <-c(proplac,propsom,propgon,propthy,propmel,propcor,proppost,propMki67,propstem,propcol,propwbc,proprbc,propendo)

#build table
props <- data.frame(Cell_Types =ids, Percent_Proportion =proportions)


#Write table
write.csv2(props,"Cell type proportions.csv")






