library(tidyverse)
library(cifti)
#install.packages("cifti") 
#78: USA (MI) [https]
#  personal library: ‘~/R/x86_64-pc-linux-gnu-library/4.0’
library(gifti)
library(XML)
library(DescTools)
library(Cairo)

# This script makes parcel-wise gene expression averages for each ROI in the 
# 200/300/400/1000 7/17 Network parcellations of Scheaffer et al (2018).
# In our work, we used 400 and 7 network solution.

# The aggregation was done across all donors and within each donor, but you need
# within donor aggregation for cell type deconvolution to control for the batch effect.
# Abagen toolbox did the sample and gene normalization within each donor, 
# so it suits our purpose well in this case.

# The number of donor contributing to and the dominant donor within each parcel is reported in csv files,
# which help you to check the donor dominance issue

# The schaefer atlas can be downloaded from here:
# https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/HCP/fslr32k/cifti

# This script is adapted from Anderson et al. NatCom 2020:
# https://github.com/HolmesLab/2020_NatComm_interneurons_cortical_function_schizophrenia/blob/master/scripts/05_foci_schaeffer_parcel_mapping.R


# project directory
base_dir  = '/Volumes/GoogleDrive/My Drive/Gradient_Shift_Cellular_Basis/Milgram/para'

# this is for testing out different parameter combo
#for ( para in c('0.1', '0.3', 'ProbeMax', 'NormZscore') ){
# this is for the final parameter combo:
for ( para in c('NormZscore0.3') ){
  # read sample information
  load(file=paste0(base_dir, '/data/ahba/donorData_obj_abagen_',para,'.Rdata'), verbose=T)
  load(file=paste0(base_dir, '/data/ahba/ahba_ctx_norm_abagen_',para,'.Rdata'), verbose=T)
  
  # read sample-to-vertex projection info
  sample_info = read_csv(paste0(base_dir, '/data/ahba/sample_info_vertex_reannot_mapped_',para,'.csv'))
  # remove samples greater than 4mm away from nearest vertex
  sample_dat  = sample_info[which(abs(sample_info$mm_to_surf) < 4),]
  
  # cortical sample by normalized gene expression data frame
  reg_micro_scale = as.data.frame(t(ctx_data_scale))
  reg_samples = sample_dat
  
  # stitch the left/right verticies together to match Schaeffer parcel cifti format
  reg_samples$bihemi_vertex = reg_samples$vertex + 1 # cifti indices index at 0, R indexes at 1
  right_ctx_idxs = intersect(grep('right', reg_samples$structure_name), which(reg_samples$top_level == 'CTX'))
  reg_samples$bihemi_vertex[right_ctx_idxs] = reg_samples$bihemi_vertex[right_ctx_idxs] + 32492
  
  # Read Schaeffer parcel info (for multiple parcel #'s, network assignments)
  donor_arr =  c('9861','10021','12876','14380','15496','15697')
  donor_specific_expression = NULL
  
  for (parcels in c('100','200','300','400','500','600','700','800','900','1000')){
    for (net in c('7','17')){
      
      # schaeffer parcellation by vertex
      schaeffer  = paste0(base_dir, '/data/Schaefer/Schaefer2018_',parcels,'Parcels_',net,'Networks_order.dscalar.nii')
      schaef_cii = read_cifti(schaeffer, drop_data = TRUE, trans_data = TRUE)
      
      # corresponding parcel labels
      schaeffer_labels = read_csv(col_names=F, paste0(base_dir, '/data/Schaefer/Schaefer2018_',parcels,'Parcels_',net,'Networks_order_info.txt'))
      schaef_labels    = schaeffer_labels$X1[grep('Network', schaeffer_labels$X1)]
      
      head(schaef_cii$data)
      length(schaef_cii$data)
      sort(unique(schaef_cii$data))
      
      # calculate the average gene expression within each parcel
      schaeffer_mat = matrix(NA, ncol=as.numeric(parcels), nrow=ncol(reg_micro_scale))
      for (donor in donor_arr){donor_specific_expression[[donor]] = schaeffer_mat}
      
      # loop over each parcel, find matching samples, calculate average expression
      for (idx in 1:length(schaef_labels)){
        write(idx,'')
        # schaeffer indices
        parcel_idxs   = which(schaef_cii$data == idx)
        # foci within this parcel
        match_idxs    = which(reg_samples$bihemi_vertex %in% parcel_idxs)
        # data for this parcel
        match_samples = reg_samples[match_idxs,]
        
        # average expression of every gene, across samples in this parcel
        schaeffer_expr = colMeans(reg_micro_scale[match_idxs,])
        
        # plug in values to the pre-allocated matrix
        schaeffer_mat[,idx] = schaeffer_expr
        
        # do the same parcel-wise averaging, but separately for each donor
        for (donor in donor_arr){
          donor_idxs     = which(as.character(reg_samples$brain) == donor)
          donor_reg_idxs = intersect(match_idxs, donor_idxs)
          expr           = colMeans(reg_micro_scale[donor_reg_idxs,])
          # average expressino of every gene, across samples in this parcel
          donor_specific_expression[[donor]][,idx] = expr
        }
      }
      # add column/row names
      schaef_out = as.data.frame(schaeffer_mat)
      rownames(schaef_out) = colnames(reg_micro_scale)
      colnames(schaef_out) = schaef_labels
      
      # write the averaged expression matrices for each donor
      for (donor in donor_arr){
        donor_specific_expression[[donor]] = as.data.frame(donor_specific_expression[[donor]])
        rownames(donor_specific_expression[[donor]]) = colnames(reg_micro_scale)
        colnames(donor_specific_expression[[donor]]) = schaef_labels
        
        donor_specific_expression[[donor]]$gene = colnames(reg_micro_scale)
        out_path = paste0(base_dir, '/data/ahba/schaeffer_donor',as.character(donor),'_p', parcels,'_',net,'Net_expr_mat_abagen_',para,'.csv')
        write_csv(donor_specific_expression[[donor]], path=out_path)
      }
      
      schaef_out$gene = colnames(reg_micro_scale)
      write_csv(schaef_out, path=paste0(base_dir, '/data/ahba/schaeffer',parcels,'_',net,'Net_expr_mat_abagen_',para,'.csv'))
    }
  }
  
  
  # calculate sample coverage for each parcel: 7 network solution
  for (parcels in c('100','200','300','400','500','600','700','800','900','1000')){
    net = '7'
    # schaeffer parcellation by vertex
    schaeffer  = paste0(base_dir, '/data/Schaefer/Schaefer2018_',parcels,'Parcels_',net,'Networks_order.dscalar.nii')
    schaef_cii = read_cifti(schaeffer, drop_data = TRUE, trans_data = TRUE)
    
    # corresponding parcel labels
    schaeffer_labels = read_csv(col_names=F, paste0(base_dir, '/data/Schaefer/Schaefer2018_',parcels,'Parcels_',net,'Networks_order_info.txt'))
    schaef_labels    = schaeffer_labels$X1[grep('Network', schaeffer_labels$X1)]
    
    head(schaef_cii$data)
    length(schaef_cii$data)
    sort(unique(schaef_cii$data))
    
    # for each parcel, calc the # of donors with a sample present, and calc how over-represented any individual donor is
    parcel_max_brain = matrix(NA, ncol=as.numeric(parcels), nrow=1)
    parcel_brain_compisition = matrix(NA, ncol=as.numeric(parcels), nrow=1)
    for (idx in 1:length(schaef_labels)){
      write(idx,'')
      parcel_idxs   = which(schaef_cii$data == idx) # schaeffer indices
      match_idxs    = which(reg_samples$bihemi_vertex %in% parcel_idxs) # foci within this parcel
      match_samples = reg_samples[match_idxs,] # data for this parcel
      
      # number of donors contributing data to this parcel
      num_brains_in_parcel = length(unique(match_samples$brain))
      
      # donor with maximal representation in parcel
      max_brain = names(rev(sort(table(match_samples$brain))))[1]
      max_brain_perc = length(which(match_samples$brain == max_brain))/length(match_samples$brain)
      
      parcel_brain_compisition[,idx] = num_brains_in_parcel
      parcel_max_brain[,idx] = max_brain_perc
    }
    colnames(parcel_brain_compisition) = schaef_labels
    colnames(parcel_max_brain) = schaef_labels
    parcel_brain_compisition = as.data.frame(parcel_brain_compisition)
    parcel_max_brain = as.data.frame(parcel_max_brain)
    
    write_csv(parcel_brain_compisition, path=paste0(base_dir, '/data/ahba/schaeffer',parcels,'_',net,'_numDonorsInEachParcel_abagen_',para,'.csv'))
    write_csv(parcel_max_brain, path=paste0(base_dir, '/data/ahba/schaeffer',parcels,'_',net,'_maxDonorPresenceInParcel_abagen_',para,'.csv'))
  }

  # calculate sample coverage for each parcel: 7 network solutions
  for (parcels in c('100','200','300','400','500','600','700','800','900','1000')){
    net = '17'
    # schaeffer parcellation by vertex
    schaeffer  = paste0(base_dir, '/data/Schaefer/Schaefer2018_',parcels,'Parcels_',net,'Networks_order.dscalar.nii')
    schaef_cii = read_cifti(schaeffer, drop_data = TRUE, trans_data = TRUE)
    
    # corresponding parcel labels
    schaeffer_labels = read_csv(col_names=F, paste0(base_dir, '/data/Schaefer/Schaefer2018_',parcels,'Parcels_',net,'Networks_order_info.txt'))
    schaef_labels    = schaeffer_labels$X1[grep('Network', schaeffer_labels$X1)]
    
    head(schaef_cii$data)
    length(schaef_cii$data)
    sort(unique(schaef_cii$data))
    
    # for each parcel, calc the # of donors with a sample present, and calc how over-represented any individual donor is
    parcel_max_brain = matrix(NA, ncol=as.numeric(parcels), nrow=1)
    parcel_brain_compisition = matrix(NA, ncol=as.numeric(parcels), nrow=1)
    for (idx in 1:length(schaef_labels)){
      write(idx,'')
      parcel_idxs   = which(schaef_cii$data == idx) # schaeffer indices
      match_idxs    = which(reg_samples$bihemi_vertex %in% parcel_idxs) # foci within this parcel
      match_samples = reg_samples[match_idxs,] # data for this parcel
      
      # number of donors contributing data to this parcel
      num_brains_in_parcel = length(unique(match_samples$brain))
      
      # donor with maximal representation in parcel
      max_brain = names(rev(sort(table(match_samples$brain))))[1]
      max_brain_perc = length(which(match_samples$brain == max_brain))/length(match_samples$brain)
      
      parcel_brain_compisition[,idx] = num_brains_in_parcel
      parcel_max_brain[,idx] = max_brain_perc
    }
    colnames(parcel_brain_compisition) = schaef_labels
    colnames(parcel_max_brain) = schaef_labels
    parcel_brain_compisition = as.data.frame(parcel_brain_compisition)
    parcel_max_brain = as.data.frame(parcel_max_brain)
    
    # updated schaefer atlas for the 17 resolution
    write_csv(parcel_brain_compisition, path=paste0(base_dir, '/data/ahba/schaeffer',parcels,'_',net,'_numDonorsInEachParcel_abagen_',para,'.csv'))
    write_csv(parcel_max_brain, path=paste0(base_dir, '/data/ahba/schaeffer',parcels,'_',net,'_maxDonorPresenceInParcel_abagen_',para,'.csv'))
  }
}