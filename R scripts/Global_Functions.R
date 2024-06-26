
#These functions are used to create the file structure and ptable files 
#used in the R shiny app ShinyApp.R
#----------------------------------------------------------------------------

# Create pathways for Labtest selection -------------------------------------
create_pathways <- function(Labtests_dict){
  #Writes a data file pathway list to be used by createDB()
  pathways = NULL
  for(i in names(Labtests_dict)){
    p <- paste0(getwd(),"/Labtests_Subsets/",i)
    pathways <- c(pathways,p)}
  return(pathways)
}

create_subsets <- function(data,Labtests_dict){
  #Create separate .rds files in different directories from the full dataset 
  for(i in names(Labtests_dict)){
    print(i)
    subset <- data %>% 
      dplyr::filter(Labtest == Labtests_dict[i])
    subset$Age <- round(subset$Age,0)
    if(i == "Hemoglobin"){subset <- subset %>% filter(LabResultValue < 130)}
    if(i == "Creatinine"){subset <- subset %>% filter(LabResultValue < 170)}
    if(i == "Potassium"){subset <- subset %>% filter(LabResultValue < 10)}
    if(i == "Leukocytes"){subset <- subset %>% filter(LabResultValue < 50)}
    if(i == "Aspartate"){subset <- subset %>% filter(LabResultValue < 350)}
    if(i == "Cholesterol"){subset <- subset %>% filter(LabResultValue < 18)}
    dir.create(paste0("Labtests_Subsets/",i))
    saveRDS(subset,paste0(getwd(),"/Labtests_Subsets/",i,"/",i,"Data.rds"))}
  return()
}

# Create ptable Database ----------------------------------------------------

createDB <- function(path){
  #Input the path to a Labtest subset .rds file (only to the folder = the name of the Labtest, it will find the file itself)
  #Function will create folder containing the ptables and one for the bigslice clusters
  #that are needed for the cluster ptables
  
  filename <- list.files(path,pattern=".rds")
  print(paste("Datafile:",filename))
  filepath <- paste0(path,"/",filename)
  print(filepath)
  data <- readRDS(filepath)
  dir.create(paste0(path,"/Ptables"))
  dir.create(paste0(path,"/BigsliceClusters/"))
  
  #Parameters
  #print(head(data,1))
  age <- c(10,20,30,40,50,60,70,80)
  #age <- c(10,20)
  sex <- c("f","m")
  #sex <- c("f")
  cluster_n <- c(400,800)
  #cluster_n <- c(400)
  cluster_type <- c("hierarchical","kmeans")
  #cluster_type <- c("hierarchical")
  test <- c("t-test","wilcox-test")

  
  createDB.app_table <- function(table){
    print(table)
    sex <- table[1]
    age <- as.numeric(table[2])
    age2 <- age+10
    test <- table[3]
    slice <- get_slice_shiny(data,sex,age,age2)
    aa <- paste0(age,"to",age2) 
    name <- paste0(path,"/Ptables/",paste(sex,aa,test,sep="_"),".txt")
    p <- tt_pval_x(slice,test=test)
    p <- p %>% 
      select("Diag","n","mu","sd","pval")
    p <- p[order(p$pval),]
    write.table(p, file=name, sep="\t", col.names=T)
    return()
  }
  
  createDB.app_table_create_clusters <- function(table){
    age <- 20
    age2 <- 80
    print(table)
    #sex <- c("m","f")
    cluster_type <- table[1]
    cluster_n <- as.numeric(table[2])
    slice1 <- get_slice_shiny(data,"f",age,age2)
    slice2 <- get_slice_shiny(data,"m",age,age2)
    slice <- rbind(slice1,slice2)
    name <- paste0(path,"/BigsliceClusters/",paste(cluster_type,cluster_n,sep="_"),".rds")
    clust <- get_clusters(data,cluster_n,cluster_type) 
    saveRDS(clust,file=name)
    return()
  }
  
  createDB.app_table_clust <- function(table){
    sex <- table[1]
    age <- as.numeric(table[2])
    age2 <- age+10
    slice <- get_slice_shiny(data,sex,age,age2)
    clust <- readRDS(table[3])
    test <- table[4]
    clust_name <- str_extract(table[3],"[[:alpha:]]+_[0-9]{3}")
    print(paste(sex,age,clust_name,test))
    aa <- paste0(age,"to",age2) 
    name <- paste0(path,"/Ptables/",paste(sex,aa,test,clust_name,sep="_"),".txt")
    p <- tt_pval_x(slice,clust = clust,test=test)
    p <- p %>% 
      select("Diag","n","mu","sd","pval")
    write.table(p, file=name, sep="\t", col.names=T)
    return()
  }
  
  param <- NULL
  for (s in sex){
    for (a in age){
      for(t in test){
        vec <- c(s,a,t)
        #print(vec)
        param <- rbind(param,vec)}}}
  param <- as.data.frame(param)
  colnames(param) <- c("sex","age","test")
  
  #Run non-cluster stuff
  print("Running ptable calculations...")
  apply(param,1,createDB.app_table) #--------------------------------------------
  
  #Create parameter matrix for bigslice generation
  param_create <- NULL

  for (t in cluster_type){
    for (k in cluster_n){
      vec <- c(t,k)
        #print(vec)
      param_create <- rbind(param_create,vec)}}
  param_create <- as.data.frame(param_create)
  colnames(param_create) <- c("cluster_type","cluster_n")
  
  #Rund bigslice clusters (Being the cluster groups themselves based on age 20-80 for 
  #male and female) grouping based on this cluster is run in an additional step
  print("Running cluster bigslice calculations")
  apply(param_create,1,createDB.app_table_create_clusters) #----------------------
  
  #Read in clusters
  cluster_list <- NULL
  for (i in dir(paste0(path,"/BigsliceClusters"))){
    if(grepl(".rds",i)!=0){
      name <- str_extract(i,"[[:alpha:]]+_[0-9]{3}")
      #name <- list.files()
      filename <- paste0(path,"/BigsliceClusters/",name,".rds")
      cluster_list <- c(cluster_list,filename)
      print(i)}}
  
  print(paste("clusters read in:",cluster_list))
  
  #Create parameter matrix to be run with clustering
  param_clust <- NULL
  for (s in sex){
    for (a in age){
      for (c in cluster_list){
        for(t in test){
          vec <- c(s,a,c,t)
          #print(paste("&&&&&&&&",vec))
          param_clust <- rbind(param_clust,vec)}}}}
  param_clust <- as.data.frame(param_clust)
  colnames(param_clust) <- c("sex","age","cluster","test")
  
  #print(param_clust)
  
  
  print("Running clusters ptables calculations...")
  apply(param_clust,1,createDB.app_table_clust) #----------------------
}


# Create Labtests Pathway Dictionary --------------------------------------

get_LabtestPW_dict <- function(){
  #Create pathways dict. Used to pair the folder/shown name of the Labtest
  #to its location
  Labtests <- c("Hemoglobin","Creatinine" ,"Potassium"  ,"Leukocytes" ,"Aspartate","Cholesterol")
  out <- NULL
  for(i in Labtests){
    out[i] <- paste0(getwd(),"/Labtests_Subsets/",i,"/",i,"Data.rds")
  }
  return(out)
}


























