###split up the fletcher data

#Script to split up the Fletcher data into male and females

#From the paper:
#scRNA-seq was done using four lanes on the 10X Genomics Chromium Single
#Cell Controller (10X Genomics, Pleasanton, CA), providing
#two technical replicates for each sex. A total of 7,861 cells
#were recovered: 1,868 and 1,861 for diestrus females, and
#2,047 and 2,085 per lane for males.

#cells all sequenced at the same time (same batch)

#4 data sets numbered 1,2,3,4.
#number tag at end of cell tag eg cell from lane 1: "AAACCTGAGACAGACC-1"


#############Housekeeping
#First clear the environment;
rm(list=ls())
#clear all plots
dev.off(dev.list()["RStudioGD"])
#Then the cache
gc()

##########Load libraries
library(svDialogs)
library(dplyr) 
library(Seurat)
library(patchwork)
library(sctransform)


#load the data in
#set the working directory to directory wish to read data from and save analysis results to
data.dir <- setwd(here::here())
########Load data
#use a dialog for user to input the folder containint the data
user.data <- dlgInput("Enter name of folder containing data to be split by sex", Sys.info()["user"])$res

#load the data
#load raw data using Read10x, this will be relative to your path, name desired folder e.g.  "Cheung Data"
raw.data <- Read10X(user.data)

#initialise a Seurat object with the raw (not normalised) data. This is a count matrix
data <- CreateSeuratObject(counts=raw.data, project = "data", min.cells = 3, min.features = 200) 


#get the identifier off each one
Sexid <- substr(colnames(data), 18,18)

#add this to the Seurat
data[["Sex"]] <- Sexid

#check this worked and sex has been added
head(x=data[[]])

#assume from paper 1&2 femals, 3&4 Males
#now subset based on sex
Female<- subset(data, Sex <3)
Male<- subset(data, Sex >2)


#sanity check all cells theres
sanity <- sum(length(colnames(Female))+length(colnames(Male)))
length(Sexid) == sanity

table(colnames(Female) %in% colnames(Male))
#have split the cell ids up and have checked all accounted for and unique to each sex

#export as two sepeare Seurat objects
saveRDS(Female, file = "FletcherFemales.rds")
saveRDS(Male, file = "FletcherMales.rds")

