
# Diag functions ---------------------------------------------------------------

get_slice <- function(data,sex,age_min,age_max,val_min=0,val_max=170){
  #Create dataslice of creatinine dataset based on age, gender and max/min LabResultValue
  d <- data %>% 
    filter(AdministrativeGender== sex) %>% 
    filter(Age %in% (age_min:age_max)) %>% 
    filter(LabResultValue %in% (val_min:val_max))
  return(d)
}
get_slice_shiny <- function(data,sex,age_min,age_max){
  #Create dataslice of creatinine dataset based on age, gender and max/min LabResultValue
  d <- data %>% 
    filter(AdministrativeGender== sex) %>% 
    filter(Age %in% (age_min:age_max)) %>% 
  return(d)
}


# Ex/Include  ------------------------------------------------------------
# Remove significantly different diags from global plot again
#Important note: If working with clusters stemming from tt_pval() use 
#exclude(dissect(tt_pval(x),T))!!! Otherwise include/exclude will only consider
#Combined occurrences of the Diags instead of each individual occurrence

exclude <- function(tes,set,only_length=F){
  #includes all rows containing any of the Diags in vector set
  d <- tes %>% 
    filter(!if_any(paste("Diag",01:05,sep="0"), function(x,set){
      return(x %in% set)}, set))
  if(only_length == T ){
    return(dim(d)[1])}
  return(d)
}

include <- function(tes,set,only_length=F){
  #includes all rows containing any of the Diags in vector set
  d <- tes %>% 
    filter(if_any(paste("Diag",01:05,sep="0"), function(x,set){
      return(x %in% set)}, set))
  if(only_length == T ){
    return(dim(d)[1])}
  return(d)
}

# Get Diags ---------------------------------------------------------------

count_slice_old <- function(x){
  #Count occurrences for each Diagnosis in a given slice
  #warning: count EVERY occurrence. Diagnoses can occur multiple times per row in slice
  out <- x %>% 
    select(LabResultValue,Diag01,Diag02,Diag03,Diag04,Diag05) %>% 
    gather(Stemsfrom,Diag,Diag01,Diag02,Diag03,Diag04,Diag05) %>% 
    group_by(Diag) %>% 
    tally() %>% 
    filter(!row_number() %in% c(1))
  return(out[order(out$n,decreasing=T),])
}

get_diags <- function(tes,only_length=F){
  #Get individual diagnoses in slice
  #additional argument to get length of diagnoses vector instead of explicit diagnoses
  out <- tes %>% 
    select(LabResultValue,Diag01,Diag02,Diag03,Diag04,Diag05) %>% 
    gather(Stemsfrom,Diag,Diag01,Diag02,Diag03,Diag04,Diag05) %>% 
    group_by(Diag) %>% 
    tally()
  out <- out[-1,] #consider that first element is the count of empty diagnoses!
  if(only_length == T){
    return(length(out$Diag))}
  return(out$Diag)
}

c_slice <- function(tes){
  #Count rows with at least one diag occurrence for every diag
  #Diags can occur multiple times per datapoint)
  d <- lapply(get_diags(tes),function(x) return(list(Diag=x,n=include(tes,x,only_length=T))))
  out <- do.call(rbind.data.frame, d)
  return(out[order(out$n,decreasing=T),])
}

slice_table <- function(x){
  #enter slice element for x to get table containing each mention of a Diag with its LabResultValue
  #(contains duplicates as up to 5 diags per one LabResultValue)
  tab <- x %>% 
    select(LabResultValue,Diag01,Diag02,Diag03,Diag04,Diag05) %>% 
    gather(Stemsfrom,Diag,Diag01,Diag02,Diag03,Diag04,Diag05) %>% 
    select(-Stemsfrom) 
  tab <- tab %>% 
    filter(Diag != "")
  tab <- tab[order(tab$Diag),]
  return(tab)
}

slice_matrix <- function(x){
  #Create overwiev of slice in matrix form. 
  #Each Diag as a column with all corresponding LabResultValues as rows
  stable <- slice_table(x)
  max.element <- count_slice_old(x)$n[1] #Get largest element  for matrix creation
  
  #create empty matrix size:(#diags * # max Diag entries)
  M <- matrix(rep(NA,length(unique(stable$Diag))*max.element[[1]]),ncol=length(unique(stable$Diag)))
  colnames(M) <- sort(unique(stable$Diag))
  M <- data.frame(M)
  
  #fill matrix. 
  for (i in seq(1,length(unique(stable$Diag)))){
    print(paste0(i,"/",length(unique(stable$Diag)),": ",sort(unique(stable$Diag))[i]))
    d <- stable %>%
      filter(Diag == sort(unique(stable$Diag))[i])
    M[1:length(d$LabResultValue),i] <- d$LabResultValue
  }
  return(M)
}

# tt_pval -----------------------------------------------------------------
#Get pvalues for t.test between Diag and global distribution without that Diag

app_tt_pval <- function(tag,tes,minimum){
  diag_i <- include(tes,tag)
  if(length(diag_i$LabResultValue) >= minimum){
    globalslice <- exclude(tes,tag)
    ou <- t.test(globalslice$LabResultValue,diag_i$LabResultValue)$p.value
    return(list(Diag=tag,n=length(diag_i$LabResultValue),mu=mean(diag_i$LabResultValue),sd=sd(diag_i$LabResultValue),pval=ou))}
}

tt_pval <- function(tes,alpha=1,HB_correction=F,minimum=5){
  #input slice for tes, slice_table(tes)
  #alpha the upper limit of pvalues to include
  #mimimum of observations for each Diag to be compared to global
  #option for Holm-Bonferri correction 
  out <- lapply(get_diags(tes),app_tt_pval,tes,minimum)
  out <- do.call(rbind.data.frame, out)
  #colnames(out) <- c("Diag","pval")
  if(HB_correction == T){
    out <- out[order(out$pval),]
    out <- out %>% 
      mutate(afraqm = seq(length(out$pval),1)) %>% 
      mutate(corr_alpha =  alpha/afraqm) %>% 
      mutate(HB_significance = pval < alpha/afraqm)
    return(out <- out[out$HB_significance == T,])} #out <- out[out$HB_significance == T,]
  return(out <- out[out$pval < alpha,])
}

#p01 <- tt_pval(m3035,0.01,HB_correction=T)
#p05 <- tt_pval(m3035,0.05,HB_correction=T)
#p0 <- tt_pval(m3035,0.1,HB_correction=T)


plot_shit_old <- function(tes,smatrix,set){
  #provide slice, slice matrix and vector of 7 Diags to plot
  #set <- head(p01[order(p01$pval),],7)$Diag
  ggplot()+
    geom_density(aes(tes$LabResultValue),lty="dashed",lwd=1.2)+
    geom_density(aes(smatrix[,set[1]],colour=paste0(set[1]," n=",sum(!is.na(smatrix[,set[1]])) )))+
    geom_density(aes(smatrix[,set[2]],colour=paste0(set[2]," n=",sum(!is.na(smatrix[,set[2]])) )))+
    geom_density(aes(smatrix[,set[3]],colour=paste0(set[3]," n=",sum(!is.na(smatrix[,set[3]])) )))+
    geom_density(aes(smatrix[,set[4]],colour=paste0(set[4]," n=",sum(!is.na(smatrix[,set[4]])) )))+
    geom_density(aes(smatrix[,set[5]],colour=paste0(set[5]," n=",sum(!is.na(smatrix[,set[5]])) )))+
    geom_density(aes(smatrix[,set[6]],colour=paste0(set[6]," n=",sum(!is.na(smatrix[,set[6]])) )))+
    geom_density(aes(smatrix[,set[7]],colour=paste0(set[7]," n=",sum(!is.na(smatrix[,set[7]])) )))+
    xlab("Creatinine LabresultValue")+
    ggtitle("Significantly different Diags")+
    theme(legend.title = element_blank())
}

plot_shit <- function(tes,p0,density=T){
  #provide slice, slice matrix and vector of 7 Diags to plot
  set <- p0$Diag
  if ( density==F){
    ggplot()+
      geom_histogram(aes(tes$LabResultValue),lty="dashed",lwd=1.2,bins=70)+
      geom_histogram(aes(include(tes,set[1])$LabResultValue,fill=paste0(set[1]," n=",p0$n[1] )),bins=70)+
      geom_histogram(aes(include(tes,set[2])$LabResultValue,fill=paste0(set[2]," n=",p0$n[2] )),bins=70)+
      xlab("Creatinine LabresultValue")+
      xlim(c(0,200))+
      ggtitle(paste("Significantly different Diags. Global n=", dim(tes)[1]))+
      theme(legend.title = element_blank())}
  else{
  ggplot()+
    geom_density(aes(tes$LabResultValue),lty="dashed",lwd=1.2)+
    geom_density(aes(include(tes,set[1])$LabResultValue,colour=paste0(set[1]," n=",p0$n[1] )))+
    geom_density(aes(include(tes,set[2])$LabResultValue,colour=paste0(set[2]," n=",p0$n[2] )))+
    geom_density(aes(include(tes,set[3])$LabResultValue,colour=paste0(set[3]," n=",p0$n[3] )))+
    geom_density(aes(include(tes,set[4])$LabResultValue,colour=paste0(set[4]," n=",p0$n[4] )))+
    geom_density(aes(include(tes,set[5])$LabResultValue,colour=paste0(set[5]," n=",p0$n[5] )))+
    geom_density(aes(include(tes,set[6])$LabResultValue,colour=paste0(set[6]," n=",p0$n[6] )))+
    geom_density(aes(include(tes,set[7])$LabResultValue,colour=paste0(set[7]," n=",p0$n[7] )))+
    xlab("Creatinine LabresultValue")+
    xlim(c(0,200))+
    ggtitle(paste("Significantly different Diags. Global n=", dim(tes)[1]))+
    theme(legend.title = element_blank())
  }
}

plot_shit <- function(tes,p0,density=T){
  #provide slice, slice matrix and vector of 7 Diags to plot
  set <- lapply(p0$Diag,function(x) str_split_1(x,pattern=" "))

  if ( density==F){
    ggplot()+
      geom_histogram(aes(tes$LabResultValue),lty="dashed",lwd=1.2,bins=70)+
      geom_histogram(aes(include(tes,set[1])$LabResultValue,fill=paste0(set[1]," n=",p0$n[1] )),bins=70)+
      geom_histogram(aes(include(tes,set[2])$LabResultValue,fill=paste0(set[2]," n=",p0$n[2] )),bins=70)+
      xlab("Creatinine LabresultValue")+
      xlim(c(0,200))+
      ggtitle(paste("Significantly different Diags. Global n=", dim(tes)[1]))+
      theme(legend.title = element_blank())}
  else{
    ggplot()+
      geom_density(aes(tes$LabResultValue),lty="dashed",lwd=1.2)+
      geom_density(aes(include(tes,set[[1]])$LabResultValue,colour=paste0(set[1]," n=",p0$n[1] )))+
      geom_density(aes(include(tes,set[[2]])$LabResultValue,colour=paste0(set[2]," n=",p0$n[2] )))+
      geom_density(aes(include(tes,set[[3]])$LabResultValue,colour=paste0(set[3]," n=",p0$n[3] )))+
      geom_density(aes(include(tes,set[[4]])$LabResultValue,colour=paste0(set[4]," n=",p0$n[4] )))+
      geom_density(aes(include(tes,set[[5]])$LabResultValue,colour=paste0(set[5]," n=",p0$n[5] )))+
      geom_density(aes(include(tes,set[[6]])$LabResultValue,colour=paste0(set[6]," n=",p0$n[6] )))+
      #geom_density(aes(include(tes,set[[7]])$LabResultValue,colour=paste0(set[7]," n=",p0$n[7] )))+
      #geom_density(aes(include(tes,set[[8]])$LabResultValue,colour=paste0(set[8]," n=",p0$n[8] )))+
      #geom_density(aes(include(tes,set[[9]])$LabResultValue,colour=paste0(set[9]," n=",p0$n[9] )))+
      #geom_density(aes(include(tes,set[[10]])$LabResultValue,colour=paste0(set[10]," n=",p0$n[10] )))+
      
      xlab("Creatinine LabresultValue")+
      xlim(c(0,200))+
      ggtitle(paste("Significantly different Diags. Global n=", dim(tes)[1]))+
      theme(legend.title = element_blank())
  }
}



# Clustering ---------------------------------------------------------------

app_ltodf <- function(x,max_length){
  #function for lapply in ltodf_dont_fill()
  if(length(x)==max_length){
    return(x)}
  else{
    diff <- max_length - length(x)
    x <- sort(x)
    x <- c(x,rep(NA,diff))
    return(x)
  }
}

ltodf_dont_fill <- function(x){
  #do.call(rbind.data.frame, x) repeats the list element when transforming it to a df
  #for all elements that arent the maximum elements length. This function fills the difference with NA
  #instead of repeating
  max_length <- max(sapply(x,length))
  x <- lapply(x,app_ltodf,max_length)
  x <- do.call(rbind.data.frame,x)
  colnames(x) <- seq(1,max_length)
  return(x)
}

get_single <- function(W,first){
  #Get matrix extract from similarity matrix W, 
  #containing all combinations from letters starting with first and second
  out <- as.data.frame(W) %>% 
    select(starts_with(c(first))) %>% 
    mutate(Names = row.names(W)) %>% 
    filter(str_detect(Names,paste0("^",first)) ) %>% 
    select(-Names)
  return(out)
}

get_double <- function(W,first,second){
  #Get matrix extract from similarity matrix W, 
  #containing all combinations from letters starting with first and second
  out <- as.data.frame(W) %>% 
    select(starts_with(c(first,second))) %>% 
    mutate(Names = row.names(W)) %>% 
    filter(str_detect(Names,paste0("^",first)) | str_detect(Names,paste0("^",second))) %>% 
    select(-Names)
  return(out)
}

get_triple <- function(W,first,second,third){
  out <- as.data.frame(W) %>% 
    select(starts_with(c(first,second,third))) %>% 
    mutate(Names = row.names(W)) %>% 
    filter(str_detect(Names,paste0("^",first)) | str_detect(Names,paste0("^",second)) | str_detect(Names,paste0("^",third)) ) %>% 
    select(-Names)
  return(out)
}

fill_N <- function(N,double){
  #Fills new matrix N with the double letter combination matrix based on row and column
  #Used in banding()
  for(i in colnames(double)){
    for(j in colnames(double)){
      N[i,j] <- double[i,j]}}
  return(N)
}

get_diags_W <- function(W){
  #Get vector of Letters of diagnoses occurring in similarity matrix W
  lettrs <- c()
  colnam <- colnames(W)
  for(i in colnam){
    let <- str_split_fixed(i,"",2)[1]
    if(!let %in% lettrs){
      lettrs <- c(lettrs,let)}}
  return(lettrs)
}

banding2 <- function(W){
  N <- matrix(0,ncol=ncol(W),nrow=nrow(W))
  colnames(N) <- colnames(W)
  rownames(N) <- rownames(W)
  N <- as.data.frame(N)
  W_diags <- get_diags_W(W)
  for( i in seq(1,length(W_diags)-1)){
    print(paste("Filling banding:",W_diags[i],W_diags[i+1]))
    N <- fill_N(N,get_double(W,W_diags[i],W_diags[i+1]))}
  return(N)
}

banding3 <- function(W){
  N <- matrix(0,ncol=ncol(W),nrow=nrow(W))
  colnames(N) <- colnames(W)
  rownames(N) <- rownames(W)
  N <- as.data.frame(N)
  W_diags <- get_diags_W(W)
  for( i in seq(1,length(W_diags)-2)){
    print(paste("Filling banding:",W_diags[i],W_diags[i+1],W_diags[i+2]))
    N <- fill_N(N,get_triple(W,W_diags[i],W_diags[i+1],W_diags[i+2]))}
  return(N)
}

get_clusters_old <- function(tes,k,banding=0,seed=1){
  #partition diagnoses of slice into k clusters based on common occurrence
  #Both word2vec and hclust are seed dependent
  set.seed(seed)
  sec <- tes %>% 
    mutate(Diags = paste(Diag01,Diag02,Diag03,Diag04,Diag05)) %>% 
    select(Diags) 
  modl <- as.matrix(word2vec(sec$Diags,type="cbow",min_count = 1))
  print(dim(modl))
  #Calculate cosine similarity
  W <- sim2(modl,method="cosine")
  diag(W) <- 0
  #transform cosine similarity scale to purely positive
  W <- W + 1
  #Order alphabetically
  W <- W[,order(colnames(W))]
  W <- W[order(row.names(W)),]
  #Remove Bullshit
  W <- W[2:dim(W)[1],2:dim(W)[2]]
  #Banding
  if(banding == 1){
    W <- banding2(W)}
  if(banding==2){
    W <- banding3(W)}
  clustrs <- kmeans(W,k)
  print(length(clustrs$cluster))
  toll <- split(names(clustrs$cluster),factor(clustrs$cluster))
  toller <- ltodf_dont_fill(toll)
  return(toller)
}

get_clusters <- function(slice,k,method="hierarchical",banding=0,group_only=F){
  #partition diagnoses of slice into k clusters based on common occurrence
  #Both word2vec and hclust are seed dependent
  set.seed(1)
  sec <- slice %>% 
    mutate(Diags = paste(Diag01,Diag02,Diag03,Diag04,Diag05)) %>% 
    select(Diags) 
  modl <- as.matrix(word2vec(sec$Diags,type="cbow",min_count = 1))
  #Calculate cosine similarity
  W <- sim2(modl,method="cosine")
  diag(W) <- 0
  #Exclude negative similarity
  W[W<0] <- 0
  #Convert to distance
  W <- 1-W
  #Order alphabetically
  W <- W[,order(colnames(W))]
  W <- W[order(row.names(W)),]
  #Remove Bullshit
  W <- W[2:dim(W)[1],2:dim(W)[2]]
  #Banding
  if(banding == 1){
    W <- banding2(W)}
  if(banding==2){
    W <- banding3(W)}
  #Clustering
  if(method == "kmeans"){
    print("kmeans")
    clustrs <- kmeans(W,k)$cluster}
  if(method == "hierarchical"){
    print("hierar")
    W <- as.dist(W)
    clustrs <- hclust(W)
    clustrs <- cutree(clustrs,k=k)}
  if(group_only==T){
    return(clustrs)}
  groups <- split(names(clustrs),factor(clustrs))
  goups <- ltodf_dont_fill(groups)
  return(goups)
}

get_sim <- function(tes,banding=0){
  #partition diagnoses of slice into k clusters based on common occurrence
  #Both word2vec and hclust are seed dependent
  #set.seed(1)
  sec <- tes %>% 
    mutate(Diags = paste(Diag01,Diag02,Diag03,Diag04,Diag05)) %>% 
    select(Diags) 
  modl <- as.matrix(word2vec(sec$Diags,type="cbow",min_count = 1))
  print(dim(modl))
  #Calculate cosine similarity
  W <- sim2(modl,method="cosine")
  diag(W) <- 0
  #Exclude negative similarity
  W[W<0] <- 0
  #Convert to distance
  W <- 1-W
  #Order alphabetically
  W <- W[,order(colnames(W))]
  W <- W[order(row.names(W)),]
  #Remove Bullshit
  W <- W[2:dim(W)[1],2:dim(W)[2]]
  #Banding
  if(banding == 1){
    W <- banding2(W)}
  if(banding==2){
    W <- banding3(W)}
  return(W)
}

get_fucked <- function(W,k,method="kmeans",group_only=F){
  #Make clusters from similarity matrix W made in get_sim()
  if(method == "kmeans"){
    print("kmeans")
    clustrs <- kmeans(W,k)$cluster}
  if(method == "hierarchical"){
    print("hierar")
    W <- as.dist(W)
    clustrs <- hclust(W)
    clustrs <- cutree(clustrs,k=k)}
  if(group_only==T){
    return(clustrs)}
  print(length(clustrs))
  groups <- split(names(clustrs),factor(clustrs))
  goups <- ltodf_dont_fill(groups)
  return(goups)
}

# tt_pval2 ----------------------------------------------------------------
#Tryna make this work with Diag clusters instead of only individual Diags

app_tt_pval <- function(tag,tes,minimum){
  diag_i <- include(tes,tag)
  if(length(diag_i$LabResultValue) >= minimum){
    globalslice <- exclude(tes,tag)
    ou <- t.test(globalslice$LabResultValue,diag_i$LabResultValue)$p.value
    return(list(Diag=paste(na.omit(as.character(as.vector(tag))),collapse=" "),n=length(diag_i$LabResultValue),mu=mean(diag_i$LabResultValue),sd=sd(diag_i$LabResultValue),pval=ou))}
}

tt_pval <- function(tes,alpha=1,HB_correction=F,clust,minimum=5){
  #input slice for tes argument: get_slice()
  #alpha = the upper limit of pvalues to include
  #mimimum of observations for each Diag to be compared to global (t.test has min of 5)
  #option for Holm-Bonferri correction 
  if(!missing(clust)){
      out <- apply(clust,1,app_tt_pval,tes,minimum)
  }else{
    out <- lapply(get_diags(tes),app_tt_pval,tes,minimum)}
  out <- do.call(rbind.data.frame, out)
  if(HB_correction == T){
    out <- out[order(out$pval),]
    out <- out %>% 
      mutate(afraqm = seq(length(out$pval),1)) %>% 
      mutate(corr_alpha =  alpha/afraqm) %>% 
      mutate(HB_significance = pval < alpha/afraqm)
    return(out <- out[out$HB_significance == T,])} 
  return(out <- out[out$pval < alpha,])
}

inspect_clust <- function(x,...){
#plots the distribution of cluster sizes
  len <- apply(x,1,function(x)sum(!is.na(x)))
  hist(len,breaks=seq(0.5,max(len)+0.5,1),xlab="Cluster Group Sizes",
       xlim=c(1,max(len)),...)
  return(len)
}

dissect <- function(p,list=T){
  #Get number of individual Diags from tt_pval object
  #list=T for output of all diags in df form for include/exclude usage
  n <- lapply(p$Diag,function(x)length(str_split_1(x," ")))
  out <- do.call(rbind.data.frame, n)
  out$names <- p$Diag
  colnames(out) <- c("n","Diag")
  if(list==T){
    l <- c()
    for(i in p$Diag){
      l <- c(l,str_split_1(i," "))}
    return(l)}
  return(out)
}

# Seed testing ------------------------------------------------------------

enemei <- function(x,y){
  #Normalized Mutual Information Criterion
  #requires output from get_clusters() with group_only=T to work
  out <- mutinformation(x,y) / sqrt(entropy(x)*entropy(y))
  return(out)
}


# Global Minus ------------------------------------------------------------

vline <- function(x = 0, color = col) {
  list(
    type = "line", 
    y0 = 0, 
    y1 = 1, 
    yref = "paper",
    x0 = x, 
    x1 = x, 
    line = list(color = "#2F4F4F")
  )
}

vline_dot <- function(x = 0, color = col) {
  list(
    type = "line", 
    y0 = 0, 
    y1 = 1, 
    yref = "paper",
    x0 = x, 
    x1 = x, 
    line = list(dash="dot",color = "grey")
  )
}

vline_dot2 <- function(x = 0, color = col) {
  list(
    type = "line", 
    y0 = 0, 
    y1 = 1, 
    yref = "paper",
    x0 = x, 
    x1 = x, 
    line = list(dash="dot",color = "green")
  )
}

globalminus_old <- function(slice,p0,main){
  #Supply slice and tt_pval() object p0 for display of global distribution and 
  #global distribution - all Diags in p0
  minus <- exclude(slice,dissect(p0,T))
  g_ichi <- ichihara(slice$LabResultValue)
  m_ichi <- ichihara(minus$LabResultValue)
  g.lower <- g_ichi$lower.limit.low
  g.higher <- g_ichi$upper.limit.high
  m.lower <- m_ichi$lower.limit.low
  m.higher <- m_ichi$upper.limit.low
  r_ich <- rnorm(dim(minus)[1],m_ichi$mean,m_ichi$sd)
  
  plot_ly(alpha=0.6,nbinsx = max(slice$LabResultValue)) %>% 
    add_histogram(x=slice$LabResultValue,name=paste("Global n=",dim(slice)[1])) %>% 
    add_histogram(x=minus$LabResultValue,name=paste("Minus n=",dim(minus)[1]))%>% 
    add_histogram(x=r_ich,name="Ichihara") %>% 
    layout(barmode="overlay",title=main,shapes=list(vline(g.lower),
                                                    vline(g.higher),
                                                    vline_dot(m.lower),
                                                    vline_dot(m.higher)))
}

globalminus2 <- function(slice,p0,main){
  #Supply slice and tt_pval() object p0 for display of global distribution and 
  #global distribution - all Diags in p0
  decimals <- max(decimal_places(slice[1:20,]$LabResultValue))
  minus <- exclude(slice,dissect(p0,T))
  g_ichi <- ichihara(slice$LabResultValue)
  m_ichi <- ichihara(minus$LabResultValue)
  g.lower <- g_ichi$lower.limit.low
  g.higher <- g_ichi$upper.limit.high
  m.lower <- m_ichi$lower.limit.low
  m.higher <- m_ichi$upper.limit.low
  r_ich <- round(rnorm(dim(minus)[1],m_ichi$mean,m_ichi$sd),decimals)
  
  p <- plot_ly(alpha=0.6,nbinsx = max(slice$LabResultValue)) %>% 
    add_histogram(x=slice$LabResultValue,name=paste("Global n=",dim(slice)[1])) %>% 
    add_histogram(x=minus$LabResultValue,name=paste("Minus n=",dim(minus)[1])) %>% 
    add_histogram(x=r_ich,name="Ichihara") %>% 
    layout(barmode="overlay",title=main,shapes=list(vline(g.lower),
                                                    vline(g.higher),
                                                    vline_dot(m.lower),
                                                    vline_dot(m.higher)))
  
  return(list(plot = p, m.ci.lower= round(m.lower,3), m.ci.higher=round(m.higher,3),
              g.ci.lower=round(g.lower,3),g.ci.higher=round(g.higher,3)))
}

globalminus_shiny <- function(slice,p0,main){
  #Supply slice and tt_pval() object p0 for display of global distribution and 
  #global distribution - all Diags in p0
  decimals <- max(decimal_places(slice[1:20,]$LabResultValue))
  minus <- exclude(slice,dissect(p0,T))
  g_ichi <- ichihara(slice$LabResultValue)
  m_ichi <- ichihara(minus$LabResultValue)
  g.lower <- g_ichi$lower.limit.low
  g.higher <- g_ichi$upper.limit.high
  m.lower <- m_ichi$lower.limit.low
  m.higher <- m_ichi$upper.limit.low
  r_ich <- round(rnorm(dim(minus)[1],m_ichi$mean,m_ichi$sd),decimals)
  
  p <- plot_ly(alpha=0.6,nbinsx = max(slice$LabResultValue)) %>% 
    add_histogram(x=slice$LabResultValue,name=paste("Global n=",dim(slice)[1])) %>% 
    add_histogram(x=minus$LabResultValue,name=paste("Minus n=",dim(minus)[1])) %>% 
    add_histogram(x=r_ich,name="Ichihara") %>% 
    layout(barmode="overlay",title=main,shapes=list(vline(g.lower),
                                                    vline(g.higher),
                                                    vline_dot(m.lower),
                                                    vline_dot(m.higher)))
  
  return(list(plot = p, m.ci.lower= round(m.lower,3), m.ci.higher=round(m.higher,3),
              g.ci.lower=round(g.lower,3),g.ci.higher=round(g.higher,3)))
}


#Cut global minus at ~125 to select non-normally distributed part of minus
#distribution. Get table and look for matches along different age groups

get_minus <- function(slice,clust,cut){
  #Extract minus population. Option for cutting to inspect outliers
  p <- tt_pval(slice,0.05,T,clust)
  minus <- exclude(slice,dissect(p,T))
  if(!missing(cut)){
    minus <- minus[minus$LabResultValue>cut,]}
  return(minus)
}

check <- function(i,bsteibl){
  s <- c()
  for(n in seq(1,dim(bsteibl)[2])){
    s <- c(s,i %in% bsteibl[,n] )}
  return(sum(s))
}

get_matches <- function(toll){
  #Get matches of outliers in outlier table above 125
  uniq <- unique(matrix(toll))
  n <- c()
  for(i in uniq){
    n <- c(n,check(i,toll))}
  summ <- as.data.frame(uniq)
  summ$n <- n -1
  summ$percent <- round(summ$n/(dim(toll)[2]-1),2)
  return(summ[order(summ$n,decreasing = T),])
}


# qval / FDR --------------------------------------------------------------

app_tt_pval <- function(tag,slice,minimum){
  diag_i <- include(slice,tag)
  if(length(diag_i$LabResultValue) >= minimum){
    globalslice <- exclude(slice,tag)
    ou <- t.test(globalslice$LabResultValue,diag_i$LabResultValue)$p.value
    return(list(Diag=paste(na.omit(as.character(as.vector(tag))),collapse=" "),n=length(diag_i$LabResultValue),mu=mean(diag_i$LabResultValue),sd=sd(diag_i$LabResultValue),pval=ou))}
}

tt_pval <- function(slice,alpha=1,MTcorr="none",clust,minimum=5){
  #input slice for tes argument: get_slice()
  #alpha = significance level
  #mimimum of observations for each Diag to be compared to global (t.test has min of 5)
  #option for Holm-Bonferri correction or FDR/qvalue correction
  if(!missing(clust)){
    out <- apply(clust,1,app_tt_pval,slice,minimum)
  }else{
    out <- lapply(get_diags(slice),app_tt_pval,slice,minimum)}
  out <- do.call(rbind.data.frame, out)
  if(MTcorr == "HB"){
    out <- out[order(out$pval),]
    out <- out %>% 
      mutate(afraqm = seq(length(out$pval),1)) %>% 
      mutate(corr_alpha =  alpha/afraqm) %>% 
      mutate(HB_significance = pval < alpha/afraqm)
    return(out <- out[out$HB_significance == T,])} 
  if(MTcorr == "FDR"){
    out <- out[order(out$pval),]
    qobj <- qvalue(out$pval, fdr.level = alpha)
    summary(qobj)
    plot(qobj)
    out <- out %>% 
      mutate(qval=qobj$qvalues, significant =  qobj$significant)
    return(out <- out[out$significant == T,])}
  out <- out[order(out$pval),]
  return(out <- out[out$pval < alpha,])
}


qval_summary <- function(qobj){
  #Capture fucking consol output of summary(qobj)
  wasfrau <- capture.output(summary(qobj))
  head <- as.vector(str_split(wasfrau[9],pattern=" ")[[1]][11:18])[-7]
  head <- c("Method",head)
  
  tab <- head
  for (i in seq(10,11)){
    #print(wasfrau[i])
    new <- read.table(textConnection(wasfrau[i]))
    tab <- rbind(tab,new)
  }
  fdr <- read.table(textConnection(wasfrau[12]))[-1]
  colnames(fdr) <- c("V1","V2","V3","V4","V5","V6","V7","V8")
  tab <- rbind(tab,fdr)
  names(tab) <- tab[1,]
  tab <- tab[-1,]
  return(tab)
}


# Sheiny functions --------------------------------------------------------

diag_compare_hist <- function(slice,tag,nbinsx=max(slice$LabResultValue)){
  diag <- include(slice,tag)
  
  plot_ly(alpha = 0.6,nbinsx=nbinsx) %>% 
    add_histogram(x=slice$LabResultValue,name=paste0("Global n=",dim(slice)[1])) %>% 
    add_histogram(x=diag$LabResultValue,name=paste0(tag," n=",dim(diag)[1])) %>% 
    layout(barmode="overlay")
}



diag_compare_density <- function(slice,tag,bw){
  diag <- include(slice,tag)
  densy <- density(diag$LabResultValue,bw=bw)
  label_global <- paste0("Global n=",dim(slice)[1])
  label_diag <- paste0(tag," n=",dim(diag)[1])
  label_title <- paste0("Relative Location of Global and ",tag," Diagnosis Distributions")
  
  plot_ly(alpha=0.6,x=slice$LabResultValue,type="histogram",name=label_global) %>% 
    add_trace(x=densy$x,y=densy$y,type="scatter",mode="lines",fill="tozeroy",yaxis="y2",name=label_diag) %>% 
    layout(yaxis2=list(overlaying="y",side="right",range=c(0,0.1)),
           title=label_title)
}

diag_compare_density2 <- function(slice,tag,bw){
  decimals <- max(decimal_places(slice[1:50,]$LabResultValue))
  dens_ylim <- if(decimals==0){0.1}else{1}
  diag <- include(slice,tag)
  densy <- density(diag$LabResultValue,bw=bw)
  label_global <- paste0("Global n=",dim(slice)[1])
  label_diag <- paste0(tag," n=",dim(diag)[1])
  label_title <- paste0("Relative Location of Global and ",tag," Diagnosis Distributions")
  
  plot_ly(alpha=0.6,x=slice$LabResultValue,type="histogram",name=label_global) %>% 
    add_trace(x=densy$x,y=densy$y,type="scatter",mode="lines",fill="tozeroy",yaxis="y2",name=label_diag) %>% 
    layout(yaxis2=list(overlaying="y",side="right",range=c(0,dens_ylim)),
           title=label_title)
}


# Globalminus Shiny -------------------------------------------------------

globalminus_shiny <- function(slice,p0,main){
  decimals <- max(decimal_places(slice[1:50,]$LabResultValue))
  minus <- exclude(slice,dissect(p0,T))
  g_ichi <- ichihara(slice$LabResultValue)
  m_ichi <- ichihara(minus$LabResultValue)
  g.lower <- g_ichi$lower.limit.low
  g.higher <- g_ichi$upper.limit.high
  m.lower <- m_ichi$lower.limit.low
  m.higher <- m_ichi$upper.limit.low
  r_ich <- round(rnorm(dim(minus)[1],m_ichi$mean,m_ichi$sd),decimals)
  name_global <- paste0("Global n=",dim(slice)[1])
  name_minus <- paste0("Minus n=",dim(minus)[1])
  name_ichi <- paste0("Healthy patients n=",dim(minus)[1])
  densy <- density(minus$LabResultValue,bw=2)
  
  p <- plot_ly(alpha=0.6,x=slice$LabResultValue,type="histogram",name=name_global,colors="blue") %>% 
    add_histogram(x=minus$LabResultValue,name=paste("Minus n=",dim(minus)[1]),marker = list(color = "#458B74")) %>% 
    add_trace(x=densy$x,y=densy$y,type="scatter",mode="lines",yaxis="y2",name=name_ichi,line=list(color="#CD6600")) %>% 
    layout(yaxis2=list(overlaying="y",side="right",range=c(0,0.1)),
           title=main,barmode="overlay",shapes=list(vline_dot(g.lower),vline_dot(g.higher),vline(m.lower),vline(m.higher)))
  
  
  return(list(plot = p, m.ci.lower= round(m.lower,3), m.ci.higher=round(m.higher,3),
              g.ci.lower=round(g.lower,3),g.ci.higher=round(g.higher,3)))
}

globalminus_shiny2 <- function(slice,p0,main){
  decimals <- max(decimal_places(slice[1:50,]$LabResultValue))
  minus <- exclude(slice,dissect(p0,T))
  g_ichi <- ichihara(slice$LabResultValue)
  m_ichi <- ichihara(minus$LabResultValue)
  g.lower <- g_ichi$lower.limit.low
  g.higher <- g_ichi$upper.limit.high
  m.lower <- m_ichi$lower.limit.low
  m.higher <- m_ichi$upper.limit.low
  r_ich <- round(rnorm(dim(minus)[1],m_ichi$mean,m_ichi$sd),decimals)
  name_global <- paste0("Global n=",dim(slice)[1])
  name_minus <- paste0("Minus n=",dim(minus)[1])
  name_ichi <- paste0("Healthy patients n=",dim(minus)[1])
  #densy <- density(minus$LabResultValue,bw=2)
  densy <- density(r_ich,bw=2)
  max_x <- m_ichi$mean*2
  
  
  p <- plot_ly(alpha=0.6,x=slice$LabResultValue,type="histogram",nbinsx=200,name=name_global,colors="blue") %>% 
    add_histogram(x=minus$LabResultValue,name=paste("Minus n=",dim(minus)[1]),marker = list(color = "#458B74")) %>% 
    add_trace(x=densy$x,y=densy$y,type="scatter",mode="lines",yaxis="y2",name=name_ichi,line=list(color="#CD6600")) %>% 
    layout(yaxis2=list(overlaying="y",side="right",range=c(0,0.1)),
           title=main,barmode="overlay",shapes=list(vline_dot(g.lower),vline_dot(g.higher),vline(m.lower),vline(m.higher)),
           xaxis=list(range=c(0,max_x)))
  
  return(list(plot = p, m.ci.lower= round(m.lower,3), m.ci.higher=round(m.higher,3),
              g.ci.lower=round(g.lower,3),g.ci.higher=round(g.higher,3)))
}


globalminus_shiny2.1 <- function(slice,p0,min_x,max_x,main){
  decimals <- max(decimal_places(slice[1:50,]$LabResultValue))
  dens_ylim <- if(decimals==0){0.1}else{1}
  minus <- exclude(slice,dissect(p0,T))
  g_ichi <- ichihara(slice$LabResultValue)
  m_ichi <- ichihara(minus$LabResultValue)
  g.lower <- g_ichi$lower.limit.low
  g.higher <- g_ichi$upper.limit.high
  m.lower <- m_ichi$lower.limit.low
  m.higher <- m_ichi$upper.limit.low
  r_ich <- round(rnorm(dim(minus)[1],m_ichi$mean,m_ichi$sd),decimals)
  name_global <- paste0("Global n=",dim(slice)[1])
  name_minus <- paste0("Minus n=",dim(minus)[1])
  name_ichi <- paste0("Healthy patients n=",dim(minus)[1])
  #densy <- density(minus$LabResultValue,bw=2)
  densy <- density(r_ich,bw=2)

  p <- plot_ly(alpha=0.6,x=slice$LabResultValue,type="histogram",name=name_global,colors="blue") %>% 
    add_histogram(x=minus$LabResultValue,name=paste("Minus n=",dim(minus)[1]),marker = list(color = "#458B74")) %>% 
    add_trace(x=densy$x,y=densy$y,type="scatter",mode="lines",yaxis="y2",name=name_ichi,line=list(color="#CD6600")) %>% 
    layout(yaxis2=list(overlaying="y",side="right",range=c(0,dens_ylim)),
           title=main,barmode="overlay",shapes=list(vline_dot(g.lower),vline_dot(g.higher),vline(m.lower),vline(m.higher)),
           xaxis=list(range=c(min_x,max_x)))
  
  return(list(plot = p, m.ci.lower= round(m.lower,3), m.ci.higher=round(m.higher,3),
              g.ci.lower=round(g.lower,3),g.ci.higher=round(g.higher,3)))
}



# Makearrow ---------------------------------------------------------------

make_arrow <- function(ci.l.g,ci.l.m,ci.u.m,ci.u.g, unit){
  low_lim <- min(ci.l.g,ci.l.m)
  upper_lim <- max(ci.u.g,ci.u.m)
  
  color1 <- if(ci.l.g < ci.l.m){"#1C86EE"}else{"#FF3030"}
  color2 <- if(ci.u.g > ci.u.m){"#1C86EE"}else{"red"}
  value1 <- round(ci.l.m - ci.l.g,2)
  value2 <- round(ci.u.m - ci.u.g,2)
  plusmin1 <- if(value1 > 0){"+"}else{}
  plusmin2 <- if(value2 > 0){"+"}else{}
  label1 <- paste(plusmin1,value1)
  label2 <-  paste(plusmin2,value2)
  
  i <- ggplot() + geom_point() + xlim(c(low_lim,upper_lim)) + ylim(1,3) +
    xlab(unit)+
    ggtitle("Change in Global vs. Minus distribution Reference Value Estimates")+
    theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
  
  i <- i + geom_segment(aes(x = ci.l.g, y = 2, xend = ci.l.m, yend = 2),
                        lineend = "round",
                        linejoin = "mitre",
                        linewidth = 1.5,
                        arrow = arrow(length = unit(0.5, "cm")),colour=color1)
  i <- i + geom_segment(aes(x = ci.u.g, y = 2, xend = ci.u.m, yend = 2),
                        lineend = "round",
                        linejoin = "mitre",
                        linewidth = 1.5,
                        arrow = arrow(length = unit(0.5, "cm")),colour=color2)
  i <- i + annotate(geom="text",x=(ci.l.g+ci.l.m)/2,y=2.65,size=7,label=label1,color=color1)
  i <- i + annotate(geom="text",x=(ci.u.g+ci.u.m)/2,y=2.65,size=7,label=label2,color=color2)
  i <- i + annotate(geom="text",x=(ci.l.g+ci.l.m)/2,y=1.5,size=4,label="Lower")
  i <- i + annotate(geom="text",x=(ci.u.g+ci.u.m)/2,y=1.5,size=4,label="Upper")
  
  return(i)
}

#make_arrow(2,3,6,5,"Geil")




