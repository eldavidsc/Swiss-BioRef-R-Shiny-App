
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

#Check for .csv files in directory
for(i in dir()){
  if(grepl(".csv",i)){
    user_data <- i
  }
}

#Read in user file
if (exists("user_data")){
  print(paste("Reading in user file:",user_data))
  data <- fread(paste0(getwd(),"/",user_data),sep=";",header=T)} else {
    warning("No .csv file detected in the directory. Make sure you copy your file in the BioRef directory in the correct format")
}

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

exists("Labtests_dict")

#Create Labtests subsets from full data file
create_subsets(data,Labtests_dict)

# Create dictionary containing pathways of the created Labtests data subsets 
Labtest_pathways <- create_pathways(Labtests_dict)
print(Labtest_pathways)

#Create DB for all specified Labtests in Labtest_pathways
print("Creating p-table libraries for:")
print(names(Labtests_dict))
for(i in Labtest_pathways){
  createDB(i)
}






