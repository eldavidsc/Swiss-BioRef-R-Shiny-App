# Quickguide Swiss BioRef Shiny App Setup
Master's thesis project of David Schaer

The interactive BioRef App aims at providing refined reference intervals for improved patient diagnosis by healthcare professionals. 
By applying the Minus Distribution Method developed in the scope of this project, this app infers upper and lower reference limits from raw routine blood analyte data, representing the range of values a "healthy" patients blood test results should lie in.

In order to use the BioRef Shiny App move your dataset into this folder as a .csv file.
By sourcing the file Initiate_Session.R your dataset will be read and several .rds subset files will be created, one for each of the labtest analytes Aspartate, Creatinine, Hemoglobin, Cholesterol, Potassium and Leukocytes. Libraries containing statistical calculation results will be created based on these subsets for each of the analytes.
When the calculations are complete you may discard the original dataset and the R shiny app is ready for use.
Start an interactive session by running ShinyApp.R
