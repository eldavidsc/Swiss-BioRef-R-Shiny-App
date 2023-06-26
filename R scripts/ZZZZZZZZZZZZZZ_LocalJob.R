
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
setwd("T:/Desktop/Data")
memory.limit(4000)
#cata <- readRDS("CreatinineData.rds")
source("T:/Desktop/BioRef/R scripts/Diag_Functions.R")
source("T:/Desktop/R Methods/Ichihara.R")

load("Labtests.RData")

for(i in pathways[c(1,6)]){
  createDB(i)
}