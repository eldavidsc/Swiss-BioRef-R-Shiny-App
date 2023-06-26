
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

#Set wd to directory this file is located in and create subfolders
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
dir.create("Labtests_Subsets")

#Read in global file
data <- fread(paste0(getwd(),"/swiss_bioref_dataset_export_20220531.csv"),sep=";",header=T)


#Source utils
source(paste0(getwd(),"/R Scripts/Diag_Functions.R"))
source(paste0(getwd(),"/R Scripts/Global_functions.R"))


#Specify Labtests of interest
Labtests_dict <- c("Hemoglobin"="Hemoglobin A1c/Hemoglobin.total in Blood by IFCC protocol",
                   "Creatinine"="Creatinine [Moles/volume] in Serum or Plasma",
                   "Potassium"="Potassium [Moles/volume] in Serum or Plasma",
                   "Leukocytes"="Leukocytes [#/volume] in Blood by Automated count",
                   "Aspartate"="Aspartate aminotransferase [Enzymatic activity/volume] in Serum or Plasma by With P-5'-P",
                   "Cholesterol"="Cholesterol [Moles/volume] in Serum or Plasma")

#Create Labtests subsets from full data file
create_subsets(Labtests_dict)

# Create dictionary containing pathways of the created Labtests data subsets 
Labtest_pathways <- create_pathways(Labtests_dict)
print(Labtest_pathways)

#Create DB for all specified Labtests in Labtest_pathways
for(i in Labtest_pathways){
  createDB(i)
}






