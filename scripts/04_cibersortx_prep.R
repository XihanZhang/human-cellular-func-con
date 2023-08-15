# rm(list=ls())

library(SingleCellExperiment)
library(Seurat)
library(mclust)
library(dplyr)
library(GEOquery)
library(EWCE)
library(tidyverse)
#BiocManager::install("biomaRt")
library(biomaRt)
#install.packages("bigmemory")
#library(bigmemory)

# This script pull out common genes from Lake data and AHBA data, and convert the format required by cibersortX.
# Write single-cell expression matrix, and AHBA mixture files.
# This script is adapted from Anderson et al. NatCom2020
# https://github.com/HolmesLab/2020_NatComm_interneurons_cortical_function_schizophrenia/blob/master/scripts/06b_cibersortx_prep.R

# loop through files derived from different abagen parameters
for ( para in c('0.1', '0.3', '0.5', 'ProbeMax', 'NormZscore','NormZscore0.3') ){
    # set up directories
    base_dir = '/Volumes/GoogleDrive/My Drive/Gradient_Shift_Cellular_Basis/Milgram'
    umi_dir  = file.path(base_dir, 'data/singlecell/')
    out_dir = '/Volumes/GoogleDrive/My Drive/Gradient_Shift_Cellular_Basis/Milgram/para'
    umi_out_dir = file.path(out_dir, 'data/singlecell/')
    
    region  = 'VisualCortex'
    in_file = paste0(umi_dir, 'lake2018/GSE97930_',region,'_snDrop-seq_UMI_Count_Matrix_08-01-2017_seurat_processed.Rdata')
    in_file
    load(verbose=T, file = in_file)
    
    #Loading objects:
    #  lake_vis_seuset
    #  lake_vis_gene_dict # 32693
    lake_vis_seuset_data = lake_vis_seuset@data # 32693 * 19368
    lake_vis_seuset_ident = lake_vis_seuset@ident # 19368
    rm(lake_vis_seuset) # save the memory
    lake_vis_seuset_data_colnames = colnames(lake_vis_seuset_data)
    
    # load AHBA data
    load(file=paste0(out_dir, '/data/ahba/donorData_obj_abagen_',para,'.Rdata'), verbose=T)
    # grab the probe information
    ahba_probes = donorData$probes
    rm(donorData) # save the memory
    
    #load(file=paste0(base_dir, '/data/ahba/ahba_data_object.Rdata'), verbose=T)
    #donor_nums = names(ahba_data)
    donor_nums = c("9861","10021","12876","14380","15496","15697")
    #rm(ahba_data)
    
    # clean the dupicates within lake and ahba
    # what merge does: The rows in the two data frames that match on the specified columns are extracted, 
    # and joined together. If there is more than one match, 
    # all possible matches contribute one row each.
    # https://stackoverflow.com/questions/24150765/why-does-merge-result-in-more-rows-than-original-data
    # Thus before merge, the duplicates from the 2 data frames needs to be cleaned
    # check if ahba probes has removed duplicates:
    unique_rows = length(unique(ahba_probes$entrez_id)) # 15459
    all_rows = nrow(ahba_probes) # 15633
    unique_rows = length(unique(lake_vis_gene_dict$entrez_id)) # 21978
    all_rows = nrow(lake_vis_gene_dict) # 32693
    # Thus in this way, merge will generate many more rows to all possible match among duplicates
    # rename the ahba gene name array
    # get rid of duplicates in ahba
    keep_idxs2     = which(!duplicated(ahba_probes$entrez_id)) # 15459
    ahba_probes = ahba_probes[keep_idxs2,] # 15459
    # get rid of duplicates in lake
    keep_idxs2     = which(!duplicated(lake_vis_gene_dict$entrez_id)) # 21978
    lake_vis_gene_dict = lake_vis_gene_dict[keep_idxs2,] # 21978
    lake_vis_seuset_data = lake_vis_seuset_data[keep_idxs2,]
    
    # merge the lake with ahba
    lake_vis_ahba_combo = merge(x=ahba_probes, y=lake_vis_gene_dict, by.x='entrez_id', by.y='entrez_id') # 14088
    
    # ENTREZ IDs that match between the two datasets
    id_keep = intersect(lake_vis_gene_dict$entrez_id, ahba_probes$entrez_id) #14088
    
    # prep SINGLE CELL SIGNATURE FILE
    lake_vis_keep_idxs = which(lake_vis_gene_dict$entrez_id %in% id_keep) # 14088 rows
    lake_vis_micro = 10^as.matrix(lake_vis_seuset_data[lake_vis_keep_idxs,]) # 14088 genes * 19368 samples
    #lake_vis_micro_orig = as.matrix(lake_vis_seuset_data[lake_vis_keep_idxs,])
    #lake_vis_micro_orig = big.matrix(lake_vis_seuset_data[lake_vis_keep_idxs,], col.names=lake_vis_seuset_data_colnames)
    
    rm(lake_vis_seuset_data) # save the memory
    
    # replace the lake gene names with the ones that have been aligned to the AHBA data
    lake_vis_subset_df = data.frame(genes=rownames(lake_vis_micro), lakeidx=1:nrow(lake_vis_micro))
    lake_vis_probe_info = merge(x=lake_vis_subset_df, y=lake_vis_ahba_combo, by.x='genes', by.y='lake_orig_Symbol')
    lake_vis_probe_info = lake_vis_probe_info[order(lake_vis_probe_info$lakeidx),]
    
    # re-arrange the row according to the order of entrez_id
    lake_vis_micro = lake_vis_micro[order(lake_vis_probe_info$entrez_id),] # data matrix must come before the probe table!
    lake_vis_probe_info = lake_vis_probe_info[order(lake_vis_probe_info$entrez_id),] # data matrix must come before the probe table!
    
    # write data
    save(lake_vis_micro, lake_vis_probe_info, lake_vis_seuset_ident, lake_vis_seuset_data_colnames, donor_nums, file=paste0(umi_out_dir, 'lake2018/lake_vis_ahba_match_new_',para,'.Rdata'))
    
    load(file=paste0(umi_out_dir, 'lake2018/lake_vis_ahba_match_new_',para,'.Rdata'), verbose=T)
    #load(file=paste0(umi_dir, 'lake2018/lake_vis_ahba_match.Rdata'), verbose=T)
    load(file=paste0(out_dir, '/data/ahba/donorData_obj_abagen_',para,'.Rdata'), verbose=T)
    
    # find index of genes that are missing in each subject
    missing_idx = c()
    for (donor in donor_nums){
        #write(donor, '')
        
        donor_ctx_idxs = intersect(which(donorData$samples$brain == donor), which(donorData$samples$top_level == 'CTX'))
        
        ahba_probes = donorData$probes
        # get rid of duplicates in ahba
        keep_idxs2  = which(!duplicated(ahba_probes$entrez_id)) # 15459
        ahba_probes = ahba_probes[keep_idxs2,] # 15459
        # subset the matrix for this subject
        ahba_ctx_dat    = 2^donorData$micro[keep_idxs2,donor_ctx_idxs]
        # subset genes common with lake data
        id_keep = intersect(lake_vis_probe_info$entrez_id, ahba_probes$entrez_id) # 13694
        keep_idxs = which(ahba_probes$entrez_id %in% id_keep) # 13694
        ahba_ctx_probes = ahba_probes[keep_idxs,]
        ahba_ctx_dat = ahba_ctx_dat[keep_idxs,]
        # re-arrange the order according to entrez_id
        ahba_ctx_dat = ahba_ctx_dat[order(ahba_ctx_probes$entrez_id), ] # data matrix must come before the probe table!
        ahba_ctx_probes = ahba_ctx_probes[order(ahba_ctx_probes$entrez_id), ] # order of probe table will change!
        
        ahba_ctx_out = cbind(ahba_ctx_probes$gene_symbol, ahba_ctx_dat)
        nan_ids = which(is.na(ahba_ctx_out), arr.ind=TRUE)
        missing_row_num = dim(nan_ids)[1]/length(donor_ctx_idxs)
        missing_rows = nan_ids[1:missing_row_num]
        missing_idx = c(missing_idx,missing_rows)
    }
    
    missing_idx = unique(missing_idx)
    missing_idx = missing_idx[!is.na(missing_idx)] # 20
    
    # remove these missing genes from all genes
    if (length(missing_idx)>0){
        lake_vis_probe_info = lake_vis_probe_info[-missing_idx,] # 14088->14068
        lake_vis_micro = lake_vis_micro[-missing_idx,] # 14088*19368 -> 14068*19368
    }
    
    # CIBERSORT formatting
    vis_dexmat = cbind(lake_vis_probe_info$gene_symbol, lake_vis_micro) #14068*19369
    
    cell_ident = unlist(lapply(lake_vis_seuset_data_colnames, function(x) strsplit(x,'_')[[1]][[1]]))
    
    colnames(vis_dexmat) = c('!Sample_title', cell_ident)
    sample_file = paste0(out_dir, '/data/cibersortX/Lake_',region,'_ahba_matched_sc_sample_orig_celllabels_new_',para,'.txt')
    write.table(x=vis_dexmat, file=sample_file, quote=FALSE, row.names = FALSE, sep='\t')
    
    cell_ident = as.character(lake_vis_seuset_ident)
    cell_ident[grep('Ex6',cell_ident)] = 'Ex6'
    cell_ident[grep('Ex3',cell_ident)] = 'Ex3'
    cell_ident[grep('Ex5',cell_ident)] = 'Ex5'
    cell_ident[grep('In1',cell_ident)] = 'In1'
    cell_ident[grep('In4',cell_ident)] = 'In4'
    cell_ident[grep('In6',cell_ident)] = 'PVALB'
    cell_ident[grep('In7|In8',cell_ident)] = 'SST'
    
    colnames(vis_dexmat) = c('!Sample_title', cell_ident)
    
    # write
    sample_file = paste0(out_dir, '/data/cibersortX/Lake_',region,'_ahba_matched_sc_sample_new_',para,'.txt')
    write.table(x=vis_dexmat, file=sample_file, quote=FALSE, row.names = FALSE, sep='\t')
    
    
    # CIBERSORT phenotyping data
    row = 1
    pheno_mat = matrix(2, nrow=length(unique(unique(lake_vis_seuset_ident))), ncol=length(lake_vis_seuset_ident))
    for (cell in unique(lake_vis_seuset_ident)){
        pheno_mat[row,lake_vis_seuset_ident == cell] = 1
        row = row + 1
    }
    uniq_cell_arr = unique(lake_vis_seuset_ident)
    out_pheno_mat = cbind(as.character(unique(lake_vis_seuset_ident)), pheno_mat)
    pheno_file = paste0(out_dir, '/data/cibersortX/Lake_',region,'_ahba_matched_sc_pheno_new_',para,'.txt')
    write.table(x=out_pheno_mat, file=pheno_file, quote=FALSE, row.names = FALSE, sep='\t', col.names=F)
    
    
    load(file=paste0(out_dir, '/data/ahba/donorData_obj_abagen_',para,'.Rdata'), verbose=T)
    
    # subset AHBA data and WRITE MIXTURE FILE
    for (donor in donor_nums){
        write(donor, '')
        
        # subset cortical samples of this donor
        donor_ctx_idxs = intersect(which(donorData$samples$brain == donor), which(donorData$samples$top_level == 'CTX'))
        
        ahba_probes = donorData$probes
        # get rid of duplicates in ahba
        keep_idxs2  = which(!duplicated(ahba_probes$entrez_id)) # 15459
        ahba_probes = ahba_probes[keep_idxs2,] # 15459
        # subset the matrix for this subject
        ahba_ctx_dat    = 2^donorData$micro[keep_idxs2,donor_ctx_idxs]
        # subset genes common with lake data
        id_keep = intersect(lake_vis_probe_info$entrez_id, ahba_probes$entrez_id) # 13694
        keep_idxs = which(ahba_probes$entrez_id %in% id_keep) # 13694
        ahba_ctx_probes = ahba_probes[keep_idxs,]
        ahba_ctx_dat = ahba_ctx_dat[keep_idxs,]
        # re-arrange the order according to entrez_id
        ahba_ctx_dat = ahba_ctx_dat[order(ahba_ctx_probes$entrez_id), ] # data matrix must come before the probe table!
        ahba_ctx_probes = ahba_ctx_probes[order(ahba_ctx_probes$entrez_id), ] # order of probe table will change!
        
        ahba_ctx_samples = donorData$samples[donor_ctx_idxs,]
        
        ahba_ctx_out = cbind(ahba_ctx_probes$gene_symbol, ahba_ctx_dat)
        colnames(ahba_ctx_out) = c('!Sample_title', ahba_ctx_samples$well_id)
        
        ciber_dir = '/Volumes/GoogleDrive/My Drive/Gradient_Shift_Cellular_Basis/Milgram/para/data/cibersortX'
        mixture_file = paste0(ciber_dir, '/', region,'_ahba_',donor,'_mixture_new_',para,'.txt')
        write(mixture_file,'')
        write.table(x=ahba_ctx_out, file=mixture_file, quote=FALSE, row.names = FALSE, sep='\t')
    }
    
    
    rm(lake_vis_micro)
    rm(lake_vis_probe_info)
    
    
    # set up directories
    base_dir = '/Volumes/GoogleDrive/My Drive/Gradient_Shift_Cellular_Basis/Milgram'
    umi_dir  = file.path(base_dir, 'data/singlecell/')
    out_dir = '/Volumes/GoogleDrive/My Drive/Gradient_Shift_Cellular_Basis/Milgram/para'
    umi_out_dir = file.path(out_dir, 'data/singlecell/')
    
    region = 'FrontalCortex'
    load(verbose=T, file = paste0(umi_dir, 'lake2018/GSE97930_',region,'_snDrop-seq_UMI_Count_Matrix_08-01-2017_seurat_processed.Rdata'))
    
    # Loading objects:
    #  lake_dfc_seuset
    #  lake_dfc_gene_dict
    # lake_dfc_gene_dict =lake_gene_dict
    lake_dfc_seuset_data = lake_dfc_seuset@data # 24654 * 10319
    lake_dfc_seuset_ident = lake_dfc_seuset@ident # 10319
    rm(lake_dfc_seuset) # save the memory
    lake_dfc_seuset_data_colnames = colnames(lake_dfc_seuset_data)
    
    # grab the probe information
    load(file=paste0(out_dir, '/data/ahba/donorData_obj_abagen_',para,'.Rdata'), verbose=T)
    ahba_probes = donorData$probes
    rm(donorData) # save the memory
    
    # check if ahba probes has removed duplicates:
    unique_rows = length(unique(ahba_probes$entrez_id)) # 15459
    all_rows = nrow(ahba_probes) # 15633
    unique_rows = length(unique(lake_dfc_gene_dict$entrez_id)) # 18352
    all_rows = nrow(lake_dfc_gene_dict) # 24654
    # Thus in this way, merge will generate many more rows to all possible match among duplicates
    # rename the ahba gene name array
    # get rid of duplicates in ahba
    keep_idxs2     = which(!duplicated(ahba_probes$entrez_id)) # 15459
    ahba_probes = ahba_probes[keep_idxs2,] # 15459
    # get rid of duplicates in lake
    keep_idxs2     = which(!duplicated(lake_dfc_gene_dict$entrez_id)) # 18352
    lake_dfc_gene_dict = lake_dfc_gene_dict[keep_idxs2,] # 18352
    lake_dfc_seuset_data = lake_dfc_seuset_data[keep_idxs2,] # 18352 * 10319
    
    # ENTREZ IDs that match between the two datasets
    id_keep = intersect(lake_dfc_gene_dict$entrez_id, ahba_probes$entrez_id) #13694
    
    # rename the ahba gene name array
    #ahba_probes$ahba_hg37_GeneNames = ahba_probes$hg37_GeneNames
    lake_dfc_ahba_combo = merge(x=ahba_probes, lake_dfc_gene_dict, by.x='entrez_id', by.y='entrez_id') #13694
    
    # prep SINGLE CELL SIGNATURE FILE
    lake_dfc_keep_idxs = which(lake_dfc_gene_dict$entrez_id %in% id_keep) # 13694
    lake_dfc_micro = 10^as.matrix(lake_dfc_seuset_data[lake_dfc_keep_idxs,]) # 13694 * 10319
    #lake_dfc_micro_orig = as.matrix(lake_dfc_seuset@data[lake_dfc_keep_idxs,])
    #lake_dfc_scale = as.matrix(lake_dfc_seuset@scale.data[lake_dfc_keep_idxs,])
    rm(lake_dfc_seuset_data) # save the memory
    
    # replace the lake gene names with the ones that have been aligned to the AHBA data
    lake_dfc_subset_df = data.frame(genes=rownames(lake_dfc_micro), lakeidx=1:nrow(lake_dfc_micro)) # 13694
    lake_dfc_probe_info = merge(x=lake_dfc_subset_df, y=lake_dfc_ahba_combo, by.x='genes', by.y='lake_orig_Symbol') # 13694
    lake_dfc_probe_info = lake_dfc_probe_info[order(lake_dfc_probe_info$lakeidx),] # 13694
    
    # re-arrange the row according to the order of entrez_id
    lake_dfc_micro = lake_dfc_micro[order(lake_dfc_probe_info$entrez_id),] # data matrix must come before the probe table!
    lake_dfc_probe_info = lake_dfc_probe_info[order(lake_dfc_probe_info$entrez_id),] # data matrix must come before the probe table!
    
    # write data
    save(lake_dfc_micro, lake_dfc_probe_info, lake_dfc_seuset_data_colnames, lake_dfc_seuset_ident, donor_nums, file=paste0(umi_out_dir, 'lake2018/lake_dfc_ahba_match_new_',para,'.Rdata'))
    
    load(file=paste0(umi_out_dir, 'lake2018/lake_dfc_ahba_match_new_',para,'.Rdata'), verbose=T)
    load(file=paste0(out_dir, '/data/ahba/donorData_obj_abagen_',para,'.Rdata'), verbose=T)
    
    # find index of genes that are missing in each subject
    missing_idx = c()
    for (donor in donor_nums){
        #write(donor, '')
        
        donor_ctx_idxs = intersect(which(donorData$samples$brain == donor), which(donorData$samples$top_level == 'CTX'))
        
        ahba_probes = donorData$probes
        # get rid of duplicates in ahba
        keep_idxs2  = which(!duplicated(ahba_probes$entrez_id)) # 15459
        ahba_probes = ahba_probes[keep_idxs2,] # 15459
        
        ahba_ctx_dat    = 2^donorData$micro[keep_idxs2,donor_ctx_idxs]
        
        id_keep = intersect(lake_dfc_probe_info$entrez_id, ahba_probes$entrez_id) # 13694
        keep_idxs = which(ahba_probes$entrez_id %in% id_keep) # 13694
        
        ahba_ctx_probes = ahba_probes[keep_idxs,]
        ahba_ctx_dat = ahba_ctx_dat[keep_idxs,]
        
        # re-arrange the order according to entrez_id
        ahba_ctx_dat = ahba_ctx_dat[order(ahba_ctx_probes$entrez_id), ] # data matrix must come before the probe table!
        ahba_ctx_probes = ahba_ctx_probes[order(ahba_ctx_probes$entrez_id), ] # order of probe table will change!
        
        ahba_ctx_out = cbind(ahba_ctx_probes$gene_symbol, ahba_ctx_dat)
        nan_ids = which(is.na(ahba_ctx_out), arr.ind=TRUE)
        missing_row_num = dim(nan_ids)[1]/length(donor_ctx_idxs)
        missing_rows = nan_ids[1:missing_row_num]
        missing_idx = c(missing_idx,missing_rows)
    }
    
    missing_idx = unique(missing_idx)
    missing_idx = missing_idx[!is.na(missing_idx)] # 19
    
    # remove these missing genes from all genes
    if (length(missing_idx)>0){
        lake_dfc_probe_info = lake_dfc_probe_info[-missing_idx,] # 13694->13675
        lake_dfc_micro = lake_dfc_micro[-missing_idx,] # 13694->13675
    }
    
    # CIBERSORT formatting
    #lake_vis_micro = 10^lake_vis_micro_orig
    #lake_vis_micro = round(lake_vis_micro,6)
    dfc_dexmat = cbind(lake_dfc_probe_info$gene_symbol, lake_dfc_micro) # 13675 10320
    #cell_ident = as.character(lake_dfc_seuset@ident)
    cell_ident = unlist(lapply(lake_dfc_seuset_data_colnames, function(x) strsplit(x,'_')[[1]][[1]]))
    
    colnames(dfc_dexmat) = c('!Sample_title', cell_ident)
    sample_file = paste0(out_dir, '/data/cibersortX/Lake_',region,'_ahba_matched_sc_sample_orig_celllabels_names_new_',para,'.txt')
    write.table(x=dfc_dexmat, file=sample_file, quote=FALSE, row.names = FALSE, sep='\t')
    
    cell_ident = as.character(lake_dfc_seuset_ident)
    cell_ident[grep('Ex6',cell_ident)] = 'Ex6'
    cell_ident[grep('Ex3',cell_ident)] = 'Ex3'
    cell_ident[grep('Ex5',cell_ident)] = 'Ex5'
    cell_ident[grep('In1',cell_ident)] = 'In1'
    cell_ident[grep('In4',cell_ident)] = 'In4'
    cell_ident[grep('In6',cell_ident)] = 'PVALB'
    cell_ident[grep('In7|In8',cell_ident)] = 'SST'
    
    colnames(dfc_dexmat) = c('!Sample_title', cell_ident)
    
    # write
    sample_file = paste0(out_dir, '/data/cibersortX/Lake_',region,'_ahba_matched_sc_sample_new_',para,'.txt')
    write.table(x=dfc_dexmat, file=sample_file, quote=FALSE, row.names = FALSE, sep='\t')
    
    # CIBERSORT phenotyping data
    row = 1
    pheno_mat = matrix(2, nrow=length(unique(unique(lake_dfc_seuset_ident))), ncol=length(lake_dfc_seuset_ident))
    for (cell in unique(lake_dfc_seuset_ident)){
        pheno_mat[row,lake_dfc_seuset_ident == cell] = 1
        row = row + 1
    }
    uniq_cell_arr = unique(lake_dfc_seuset_ident)
    out_pheno_mat = cbind(as.character(unique(lake_dfc_seuset_ident)), pheno_mat)
    pheno_file = paste0(base_dir, '/data/cibersortX/Lake_',region,'_ahba_matched_sc_pheno_new_',para,'.txt')
    write.table(x=out_pheno_mat, file=pheno_file, quote=FALSE, row.names = FALSE, sep='\t', col.names=F)
    
    load(file=paste0(out_dir, '/data/ahba/donorData_obj_abagen_',para,'.Rdata'), verbose=T)
    
    # subset AHBA data and WRITE MIXTURE FILE
    for (donor in donor_nums){
        write(donor, '')
        
        donor_ctx_idxs = intersect(which(donorData$samples$brain == donor), which(donorData$samples$top_level == 'CTX'))
        
        ahba_probes = donorData$probes
        # get rid of duplicates in ahba
        keep_idxs2  = which(!duplicated(ahba_probes$entrez_id)) # 15459
        ahba_probes = ahba_probes[keep_idxs2,] # 15459
        
        ahba_ctx_dat    = 2^donorData$micro[keep_idxs2,donor_ctx_idxs]
        
        id_keep = intersect(lake_dfc_probe_info$entrez_id, ahba_probes$entrez_id) # 13694
        keep_idxs = which(ahba_probes$entrez_id %in% id_keep) # 13694
        
        ahba_ctx_probes = ahba_probes[keep_idxs,]
        ahba_ctx_dat = ahba_ctx_dat[keep_idxs,] # 13675
        ahba_ctx_samples = donorData$samples[donor_ctx_idxs,]
        
        # re-arrange the order according to entrez_id
        ahba_ctx_dat = ahba_ctx_dat[order(ahba_ctx_probes$entrez_id), ] # data matrix must come before the probe table!
        ahba_ctx_probes = ahba_ctx_probes[order(ahba_ctx_probes$entrez_id), ] # order of probe table will change!
        
        ahba_ctx_out = cbind(ahba_ctx_probes$gene_symbol, ahba_ctx_dat)
        
        colnames(ahba_ctx_out) = c('!Sample_title', ahba_ctx_samples$well_id)
        
        ciber_dir = '/Volumes/GoogleDrive/My Drive/Gradient_Shift_Cellular_Basis/Milgram/para/data/cibersortX'
        mixture_file = paste0(ciber_dir, '/', region,'_ahba_',donor,'_mixture_new_',para,'.txt')
        write(mixture_file,'')
        write.table(x=ahba_ctx_out, file=mixture_file, quote=FALSE, row.names = FALSE, sep='\t')
    }
}