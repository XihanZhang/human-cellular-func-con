# Load packages
# install.packages("WGCNA")
# install.packages("BiocManager")
# BiocManager::install("GO.db")
# BiocManager::install("impute")
# BiocManager::install("preprocessCore")
# rm(list=ls())
library(tidyverse)
library(WGCNA)


# set up the base and data dir
base_dir   = '/Volumes/GoogleDrive/My Drive/Gradient_Shift_Cellular_Basis/Milgram'
data_path  = paste(base_dir, '/para/data/ahba', sep = '')

# loop through files derived from different abagen parameters
for ( para in c('0.1', '0.3', '0.5', 'ProbeMax', 'NormZscore','NormZscore0.3') ){
  # read the expression matrix and annotation info
  fname      = paste0(data_path, '/ahba_group_samples_expression_', para,'.csv')
  microdata  = read_csv(fname)
  fname      = paste0(data_path, '/sample_info_vertex_reannot_mapped_', para,'.csv')
  sampleInfo  = read_csv(fname)
  fname      = paste0(data_path, '/ahba_expression_EntrezID_GeneName_Chromosome_full_', para,'.csv')
  AnnotationInfo  = read_csv(fname)
  fname      = paste(base_dir, 'microarray','normalized_microarray_donor9861','Ontology.csv', sep='/')
  ontology  = read_csv(fname)

  # convert the expression to matrix
  microdata_mat = microdata %>% as.data.frame() %>% data.matrix() %>% t()

  # assemble them into Robject
  donorData = NULL
  donorData$micro    = microdata_mat
  donorData$samples  = sampleInfo
  donorData$ontology = ontology
  donorData$probes   = AnnotationInfo

  save(x=donorData, file=paste0(data_path, '/donorData_obj_abagen_', para,'.Rdata'))
}