
########################## Uses multiple database in BioRef folder
memory.limit(4000)
set.seed(1)

library(shiny)
library(tidyverse)
library(plotly)
library(scrutiny)
library(qvalue)
library(shinyWidgets)


#Set wd to directory this file is in 
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

source(paste0(getwd(),"/R Scripts/Diag_Functions.R"))
source(paste0(getwd(),"/R Scripts/Global_Functions.R"))
source(paste0(getwd(),"/R Scripts/Ichihara.R"))




### Input Variables
sex_show <- c("Male","Female")
agesel_dict <- c("10-20"=10,"20-30"=20,"30-40"=30,"40-50"=40,"50-60"=50,"60-70"=60,
                 "70-80"=70,"80-90"=80)
core_test <- c("t-test","wilcox-test")
clustering_type <- c("Hierarchical","Kmeans")
number_cluster <- c(400,800)
significance_level_show <- c(0.05,0.01)
Labtests <- c("Hemoglobin","Creatinine" ,"Potassium"  ,"Leukocytes" ,"Aspartate","Cholesterol")
LabtestsPW_dict <- get_LabtestPW_dict()


### UI ------------------------------------------------------------------------------
cluster_tabs <- tabsetPanel(
  id = "clustering_tabs",
  tabPanel("Single",
           radioButtons("sex","Biological sex",sex_show),
           selectInput("age","Select age range",selected=30,agesel_dict),
           radioButtons("test","Select test",core_test),
           
           selectInput("sign","Significance Level",significance_level_show)),
  tabPanel("Clusters", 
           radioButtons("sex2","Biological sex",sex_show),
           selectInput("age2","Select age range",selected=30,agesel_dict),
           radioButtons("test2","Select test",core_test),
           selectInput("sign2","Significance Level",significance_level_show),
           radioButtons("clust_type","Clustering Type",clustering_type),
           selectInput("numclust","Number of clusters",number_cluster))
)

ui <- fluidPage(
  titlePanel("Swiss BioRef: Reference Interval Estimation"),
  wellPanel(textOutput("try")),
  sidebarLayout(
    sidebarPanel(selectInput("Labtests","Analyte",selected = "Creatinine",Labtests),
                 cluster_tabs),
    mainPanel(dropdownButton(
      tags$h3("Histogram Options"),
      sliderInput("xminmax","Select range to be plotted",min=0,max=200,value=c(0,170)),
      circle = TRUE,
      status = "danger", 
      icon = icon("gear"), width = "300px",
      tooltip = tooltipOptions(title = "Click to see inputs !")),
      plotlyOutput("globalminus",width="850px"),
              div(style = "height: 120px; width: 650px;",
                  # Plot output
                  plotOutput("arrows")),)),

  wellPanel(textOutput("ci")),
  #######################################################################################
  tabsetPanel(
    tabPanel("Browse Significant Diagnoses",
             sidebarLayout(sidebarPanel(textOutput("n_diag_slice"),
                                        textOutput("n_diag_sign"),
                                        textOutput("n_diag_ratio")),
                           mainPanel(dataTableOutput("ptable")))),
    tabPanel("Inspect Single Diagnose",
             sidebarLayout(sidebarPanel(uiOutput("diag_inspect"),
                                        sliderInput("bandwidth","Density Bandwidth",0.1,5,2)),
                           mainPanel(plotlyOutput("diag_compare_density_plot")))),
    tabPanel("Multiple Testing Summary",
             navlistPanel(
               id = "tabset",
               "False Discovery Rate",
               tabPanel("FDR Plots",
                        h3(textOutput("FDR_pi0")),
                        plotOutput("FDR_plot")),
               tabPanel("Rejection Overview",
                        textOutput("tab2_nav1"),
                        dataTableOutput("FDR_summary")))))
)
### Sörvör --------------------------------------------------------------------------
server <- function(input, output) {
  
  #################### Data readIn
  pathtodata <- reactive({LabtestsPW_dict[[input$Labtests]]})
  data <- reactive({readRDS(pathtodata())})

  #################### Basic stuff
  #Sex
  sex <- reactive({if(input$sex == "Male" | input$sex2 == "Male"){"m"}
    else if(input$sex == "Female" | input$sex2 == "Female"){"f"} })
  observeEvent(input$sex,{updateRadioButtons(inputId = "sex2",selected=input$sex)}) #Change clust tabset duplicate to same selection
  observeEvent(input$sex2,{updateRadioButtons(inputId = "sex",selected=input$sex2)})#Change noclust tabset duplicate to same selection
  #Age
  age <- age2 <- reactive({as.numeric(input$age)})
  observeEvent(input$age,{updateSelectInput(inputId = "age2", selected=input$age)})
  observeEvent(input$age2,{updateSelectInput(inputId = "age", selected=input$age2)})
  age2 <- reactive({age()+10})
  #sign
  sign <- sign2 <- reactive({input$sign})
  observeEvent(input$sign,{updateSelectInput(inputId = "sign2", selected=input$sign)})
  observeEvent(input$sign2,{updateSelectInput(inputId = "sign", selected=input$sign2)})
  ss <- reactiveVal()
  #test
  test <- test2 <- reactive({input$test})
  observeEvent(input$test, {updateSelectInput(inputId = "test2", selected=input$test)})
  observeEvent(input$test2, {updateSelectInput(inputId = "test", selected=input$test2)})
  
  #Cluster 
  v <- reactiveValues(clustering_active = F)
  observeEvent(input$clustering_tabs,{
    if(input$clustering_tabs == "Single"){
      v$clustering_active <- F}else{
        v$clustering_active <- T}})
  cluster_type <- reactive({input$clust_type})
  cluster_n <- reactive({input$numclust}) 
  
  #################### File Selection: choose which raw table is read in
  name <- reactive({if(v$clustering_active){paste(sex(),paste0(age(),"to",age2()),test(),cluster_type(),cluster_n(),sep="_")
  }else{paste(sex(),paste0(age(),"to",age2()),test(),sep="_")}})
  path <- reactive({paste0(getwd(),"/Labtests_Subsets/",input$Labtests,"/Ptables/")})
  filename <- reactive(paste0(name(),".txt"))
  pathway <- reactive({paste0(path(),filename())})
  
  #################### FDR Object
  ptable_raw <- reactive({read.table(pathway())})
  FDR_qobj <- reactive({qvalue(ptable_raw()$pval,fdr.level=sign())})
  FDR_summary <- reactive({qval_summary(FDR_qobj())})
  FDR_plot <- reactive(plot(FDR_qobj()))
  FDR_pi0 <- reactive(FDR_qobj()$pi0)
  
  #################### FDR Correction
  ptable_FDR_full <- reactive({if(sign()==0.05){
    ptable_raw() %>% mutate(test = FDR_qobj()$pvalues, qval= FDR_qobj()$qvalues,significant = qval < 0.05) %>% filter(significant==TRUE)}
    else if(sign() == 0.01){
      ptable_raw() %>% mutate(test = FDR_qobj()$pvalues, qval= FDR_qobj()$qvalues,significant = qval < 0.01) %>% filter(significant==TRUE)}})
  ptable_FDR <- reactive({ptable_FDR_full() %>% select(Diag,n,mu,sd,pval)})
  sig_diag_list <- reactive({ ptable_FDR()[order(ptable_FDR()$Diag),] %>% select(Diag) })
  
  #################### FDR output
  output$FDR_summary <- renderDataTable(FDR_summary(),options=list(autoWidth=T,dom="t"))
  output$ptable <- renderDataTable(ptable_FDR(),options = list(scrollX=T,autoWidth=T,iDisplayLength = 10))
  output$FDR_plot <- renderPlot(FDR_plot()) ###
  output$FDR_pi0 <- renderText({paste("Estimated proportion of true null hypothesis pi0:",round(FDR_pi0(),3))})
  
  #################### Slice & Glomalminus
  slice <- reactive({get_slice_shiny(data(),sex(),age(),age2())}) #### <- tise
  globalminus <- reactive({globalminus_shiny(slice(),ptable_FDR(),input$xminmax[1],input$xminmax[2],name())})
  ci <- reactive({paste("Reference Interval Estimate:",globalminus()$m.ci.lower,"-",globalminus()$m.ci.higher,slice()$LabResultUnit[1])})
  g.ci.lower <- reactive({globalminus()$g.ci.lower})
  g.ci.higher <- reactive({globalminus()$g.ci.higher})
  m.ci.lower <- reactive({globalminus()$m.ci.lower})
  m.ci.higher <- reactive({globalminus()$m.ci.higher})
  
  #################### Diag Table stuff (Tab 1)
  n_diag_slice <- reactive({length(unique(ptable_raw()$Diag))})
  n_diag_sign <- reactive({dim(sig_diag_list())[1]})
  n_diag_ratio <- reactive({round((n_diag_sign()/n_diag_slice())*100,2)})
  output$n_diag_slice <- renderText({paste0("Total diagnoses in selected slice: ",n_diag_slice())})
  output$n_diag_sign <- renderText({paste0("Thereof significantly different from global: ",n_diag_sign())})
  output$n_diag_ratio <- renderText({paste0("Percentage: ",n_diag_ratio(),"%")})
  
  #################### Diag inspect stuff (Tab 2)
  output$diag_inspect <- renderUI({selectInput("diag_inspect","Select Diag for inspection",sig_diag_list())})
  #output$diag_inspect <- renderUI({pickerInput("diag_inspect","Select Diag for inspection",sig_diag_list())})
  
  output$mueter <- renderText(input$diag_inspect)
  diag_compare_density_plot <- reactive({diag_compare_density(slice(),input$diag_inspect,input$bandwidth)}) 
  output$diag_compare_density_plot <- renderPlotly(diag_compare_density_plot())
  
  #################### Text output & Stuff
  output$tab2_nav1 <- renderText({"Number of significant Diagnoses per significance level "})
  output$globalminus <- renderPlotly(globalminus()$plot)
  output$ci <- renderText(ci())
  arrowplot <- reactive({make_arrow(globalminus()$g.ci.lower,globalminus()$m.ci.lower,
                                    globalminus()$m.ci.higher,globalminus()$g.ci.higher,slice()$LabResultUnit[1])})
  output$arrows <- renderPlot({arrowplot()}, height = 120,width=650)
  
  #### try
  output$try <- renderText(pathway())
}


# Run the application 
shinyApp(ui = ui, server = server)













