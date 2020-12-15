####ReadME
#Script to detect gene ontology, accounts for selection bias (Young MD, Wakefield MJ, Smyth GK, Oshlack A (2010). “Gene ontology analysis for RNA-seq: accounting for selection bias.” Genome Biology, 11, R14. ).
#follow goseq [vignette](http://bioconductor.org/packages/release/bioc/vignettes/goseq/inst/doc/goseq.pdf)
#Gene Ontology branches; Cellular Components, BiologicalProcesses and Molecular Functions.
#############Housekeeping
#First clear the environment;
rm(list=ls())
#Then the cache
gc()

##########Load libraries
library(dplyr) 
library(goseq)
library(GO.db)

##########Set seed and working directories
#set seed so plots same each time
set.seed(1234)

#set the working directory to directory wish to read data from and save analysis results to
data.dir <- setwd(here::here())

###########Prep the gene data

#####Load in the gene tables
de.genes <- read.csv(file="DifferentialyExpressedGenes.csv")
assayed.genes <-  read.csv(file="AllGenes.csv")

#extract only the genes
de.genes <-de.genes$gene
assayed.genes <- assayed.genes$x

#construct vector for goseq 
genes=as.integer(assayed.genes%in%de.genes)
names(genes)=assayed.genes
head(genes)

#tell goseq what genome and ID format were used to summarise the data
#use mm9 build of mouse genome
#check which gene ID code this corresponds to
#supportedOrganisms()[supportedOrganisms()$Genome=="mm9",]
# Ensembl format so the gene ID code to use is "ensGene".
# Entrez format so the gene ID code to use is "knownGene".

##########Perform GO analysis

#fit probability weighting function (PWF) 
#obtain weighting of each gene depending on its length to remove length bias.
#using gene symbols
pwf=nullp(genes,"mm9","geneSymbol")
head(pwf)

##########Option One
#use Wallenius approximation to calculate over and under expressed GO categories amonsgst DE genes.
#goseq will automatically fetch the data describing relationship between the gene IDs and the GO categories
GO.wall=goseq(pwf,"mm9","geneSymbol")
head(GO.wall)

########Option Two
#Use random sampling to build null distribution
#accuracy of the sampling method is limited by the number of samples generated

GO.samp=goseq(pwf,"mm9","geneSymbol",method="Sampling",repcnt=1000)
head(GO.samp)

########Make sense of GO analysis results

#identify categories significantly enriched/unenriched below p=value cut-off

#overenriched using 0.5 FDR (Benjamini and Hochberg 1995)
enriched.GO=GO.wall$category[p.adjust(GO.wall$over_represented_pvalue,method="BH")<.05]
head(enriched.GO)

#save the information to a txt file
sink(file="GO terms.txt")
#Pull information on thes GO terms
for(go in enriched.GO[1:10]){
           print(GOTERM[[go]])
           cat("--------------------------------------\n")
}
sink(file=NULL)


############# KEGG pathway analysis

#Map genes to KEGG pathways
# ENSEMBL<->Entrez and Entrez<->KEGG mappings.

# Get the mapping from ENSEMBL 2 Entrez
en2eg=as.list(org.Mm.egENSEMBL2EG)
# Get the mapping from Entrez 2 KEGG
eg2kegg=as.list(org.Mm.egPATH)
# Define a function which gets all unique KEGG IDs> # associated with a set of Entrez IDs
grepKEGG=function(id,mapkeys){unique(unlist(mapkeys[id],use.names=FALSE))}
# Apply this function to every entry in the mapping from
# ENSEMBL 2 Entrez to combine the two maps
kegg=lapply(en2eg,grepKEGG,eg2kegg)
head(kegg)

#produce the PWF
pwf=nullp(genes,"mm9","geneSymbol")

#perform KEGG analysis
KEGG=goseq(pwf,gene2cat=kegg)
head(KEGG)
