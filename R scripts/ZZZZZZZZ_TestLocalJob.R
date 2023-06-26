
memory.limit(4000)

library(tidyverse)
library(plotly)
require(data.table)
library(reshape2)
library(cowplot)
library(lsa)
library(word2vec)
library(text2vec)
library(Matrix)
library(mclust)
library(infotheo)
library(pheatmap)
library(qvalue)
library(scrutiny)

#Set wd to directory this file is in and create subfolders
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
getwd()
dir.create("Labtests_Subsets")

#Read in global file
#data <- fread("swiss_bioref_dataset_export_20220531.csv",sep=";",header=T)
#data <- readRDS("T:/Desktop/Data/FullData.rds")

#Source utils
source(paste0(getwd(),"/R Scripts/Global_Functions.R"))
source(paste0(getwd(),"/R Scripts/Diag_Functions.R"))
source(paste0(getwd(),"/R Scripts/Ichihara.R"))
load("Labtests.RData")

#createDB(Labtests_dict)

# For loop for running createDB() on every LabTest ------------------------

pathways = NULL
for(i in names(Labtests_dict)){
  p <- print(paste0(getwd(),"/Labtests_Subsets/",i))
  pathways <- c(pathways,p)
}
pathways

for(i in pathways[c(2,6)]){
  createDB(i)
}



#saveRDS(cata,paste0(getwd(),"/Labtests_Subsets/","Creatinine","/","Creatinine","Data.rds"))













