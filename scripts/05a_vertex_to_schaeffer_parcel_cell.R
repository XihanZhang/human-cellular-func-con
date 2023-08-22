library(tidyverse)
library(gplots)
library(RColorBrewer)
library(Cairo)
library(cifti)
library(cowplot)

# set up directories
base_dir  = '/Volumes/GoogleDrive/My Drive/Gradient_Shift_Cellular_Basis/Milgram/para'
old_base_dir  = '/Volumes/GoogleDrive/My Drive/Gradient_Shift_Cellular_Basis/Milgram'
ciber_dir = paste0(base_dir, '/data/cibersortX')
load(file='/Volumes/GoogleDrive/My Drive/Gradient_Shift_Cellular_Basis/Milgram/data/ahba/ahba_data_object.Rdata', verbose=T)

# function to parcel vertex-registered sample-level cell type fractions according to schaefer atlas
parcel_average = function(ciber_dat, ref_region, parcels, net, cell_types, name_string=''){
    donor_specific_expression = NULL
    
    # schaeffer parcellation by vertex
    schaeffer  = paste0(base_dir, '/data/Schaefer/Schaefer2018_',parcels,'Parcels_',net,'Networks_order.dscalar.nii')
    schaef_cii = read_cifti(schaeffer, drop_data = TRUE, trans_data = TRUE)
    # output the fsa6 label to Enigma toolbox for later plotting in python adn matlab
    write.table(schaef_cii$data,paste0(old_base_dir, '/software/ENIGMA-1.1.3/matlab/shared/parcellations/schaefer_',parcels,'_fsa6.csv'), row.names = FALSE,col.names=FALSE)
    write.table(schaef_cii$data,paste0(old_base_dir, '/software/ENIGMA-1.1.3/enigmatoolbox/datasets/parcellations/schaefer_',parcels,'_fsa6.csv'), row.names = FALSE,col.names=FALSE)
    
    
    # corresponding parcel labels
    schaeffer_labels = read_csv(col_names=F, paste0(base_dir, '/data/Schaefer/Schaefer2018_',parcels,'Parcels_',net,'Networks_order_info.txt'))
    schaef_labels    = schaeffer_labels$X1[grep('Network', schaeffer_labels$X1)]
    write.table(schaef_labels,paste0(old_base_dir, '/data/Gradients_Margulies2016/schaefer_',parcels,'_',net,'Net_labels.csv'), row.names = FALSE,col.names=FALSE)
    
    # calculate the average cell fractions within each parcel
    schaeffer_mat = matrix(NA, ncol=as.numeric(parcels), nrow=length(cell_types))
    for (donor in donor_arr){donor_specific_expression[[donor]] = schaeffer_mat}
    
    reg_cell_types = paste0(ref_region, '_', cell_types)
    
    for (idx in 1:length(schaef_labels)){
        write(idx,'')
        parcel_idxs   = which(schaef_cii$data == idx) # schaeffer indices
        match_idxs    = which(ciber_dat$bihemi_vertex %in% parcel_idxs) # foci within this parcel
        match_samples = ciber_dat[match_idxs,] # data for this parcel
        
        schaeffer_expr = colMeans(ciber_dat[match_idxs, reg_cell_types], na.rm=T) # average expressino of every gene, across samples in this parcel
        schaeffer_mat[,idx] = schaeffer_expr # plug in values to the pre-allocated matrix
        
        for (donor in donor_arr){
            donor_idxs     = which(as.character(ciber_dat$brain) == donor)
            donor_reg_idxs = intersect(match_idxs, donor_idxs)
            expr = colMeans(ciber_dat[donor_reg_idxs,reg_cell_types]) # average expressino of every gene, across samples in this parcel
            donor_specific_expression[[donor]][,idx] = expr
        }
    }
    schaef_out = as.data.frame(schaeffer_mat)
    rownames(schaef_out) = reg_cell_types
    colnames(schaef_out) = schaef_labels
    
    for (donor in donor_arr){
        donor_specific_expression[[donor]] = as.data.frame(donor_specific_expression[[donor]])
        rownames(donor_specific_expression[[donor]]) = reg_cell_types
        colnames(donor_specific_expression[[donor]]) = schaef_labels
        
        donor_specific_expression[[donor]]$gene = reg_cell_types
        write_csv(donor_specific_expression[[donor]], path=paste0(base_dir, '/data/ahba/schaeffer_LAKE_',ref_region,'_donor',as.character(donor),'_p', parcels,'_',net,'Net_expr_mat',name_string,'_new_',para,'.csv'))
    }
    
    schaef_out$gene = reg_cell_types
    write_csv(schaef_out, path=paste0(base_dir, '/data/ahba/schaeffer_LAKE_',ref_region,'_',parcels,'_',net,'Net_expr_mat',name_string,'_new_',para,'.csv'))
    
    return(list(schaef_out, donor_specific_expression))
}

######### Frontal ############
region = 'FrontalCortex'
# for different parameter combinations
#for ( para in c('0.1','0.3','0.5','ProbeMax', 'NormZscore', 'NormZscore0.3') ){
# for the finalized parameter combo
for ( para in c('NormZscore0.3') ){
    for (parcels in c('100','200','300','400','500','600','700','800','900','1000')){
        parcel_num = strtoi(parcels)
        for (net in c('7','17')){
            # read cell type expression of "signature genes"
            gep_file = paste0(ciber_dir, '/CIBERSORTx_Lake_',region,'_ahba_matched_sc_SignatureMatrix_new_',para,'.txt')
            gep = read.table(gep_file, header=T)
            
            # create heatmap of gene signatures
            file_out = paste0(paste0(base_dir, '/figures/',region,'_gep_signature_matrix_new_',para,'.pdf'))
            CairoPDF(file_out, height=5, width=4)
            heatmap.2(as.matrix(gep[,2:ncol(gep)]), col=brewer.pal(11,"RdBu"), key=FALSE, dendrogram='col', scale="row", trace="none", labRow=FALSE)
            dev.off()
            file_out
            
            # average expression of signature genes in each cell class
            sig_cormat = cor(gep[,2:ncol(gep)], use = 'pairwise.complete')
            cell_text = t(as.matrix(sig_cormat))
            cell_text = round(cell_text,2)
            # color pallette for the correlation plot
            my_palette = colorRampPalette(c('#67001F','#B2182B','#D6604D','#F4A582','#FDDBC7','#F7F7F7','#D1E5F0','#92C5DE','#4393C3','#2166AC','#053061'))(n = 299)
            my_palette = rev(my_palette)
            
            # plot
            pdf(file = paste0(base_dir, '/figures/lake_',region,'_corr_of_signature_mat_new_',para,'.pdf'), width=10, height=10)
            heatmap.2(as.matrix(sig_cormat), cellnote=cell_text, trace="none", notecol='black',col=my_palette,breaks=seq(-1, 1, length.out=300),
                      distfun=function(x) dist(x, method="euclidean"),dendrogram="col",density.info="none",
                      hclustfun=function(x) hclust(x, method="ward.D2"),Rowv=T, symm=F,symkey=F,symbreaks=T)
            dev.off()
            
            # read AHBA informatioN
            load(file=paste0(base_dir, '/data/ahba/ahba_ctx_norm_abagen_',para,'.Rdata'), verbose=T)
            donor_arr = c("9861","10021","12876","14380","15496","15697")
            
            # read sample-to-vertex projection info
            sample_info = read_csv(paste0(base_dir, '/data/ahba/sample_info_vertex_reannot_mapped_',para,'.csv'))
            sample_dat  = sample_info[which(abs(sample_info$mm_to_surf) < 4),]
            
            write(region,'')
            all_ciber = NULL
            for (donor in donor_arr){
                write(donor, '')
                
                path_to_ciber = paste0(ciber_dir, '/',region,'_ahba_',donor,'_fraction_new_',para,'.txt')
                ciber_dat = read_delim(path_to_ciber, delim='\t') #read.table(path_to_ciber, header=T)
                all_ciber = rbind(all_ciber, ciber_dat)
                
                ciber_field = paste0(region, '_CiberFrac')
                ahba_data[[donor]][[ciber_field]] = merge(x=ciber_dat, y=sample_info, by.y='well_id', by.x='Mixture')
            }
            write_csv(x=all_ciber, path=paste0(ciber_dir, '/Cibersortx_AHBAsamples_Lake_',region,'_new_',para,'.csv'))
            
            # Read cell fractions derived from Visual Cortex cell signatures
            ciber = read_csv(paste0(ciber_dir, '/Cibersortx_AHBAsamples_Lake_',region,'_new_',para,'.csv'))
        
            colnames(ciber) = paste0('DFC_',colnames(ciber))
            
            # merge cibersort data with AHBA sample dat
            cell_types = colnames(ciber)[2:19]
            sample_dat_ciber = merge(x=sample_dat, y=ciber, by.x='well_id', by.y='DFC_Mixture', all.x=T)
            
            # Collapse individual cell fractions into schaeffer parcels
            # stitch the left/right verticies together to match Schaeffer parcel cifti format
            sample_dat_ciber$bihemi_vertex = sample_dat_ciber$vertex + 1 # cifti indices index at 0, R indexes at 1
            right_ctx_idxs = intersect(grep('right', sample_dat_ciber$structure_name), which(sample_dat_ciber$top_level == 'CTX'))
            sample_dat_ciber$bihemi_vertex[right_ctx_idxs] = sample_dat_ciber$bihemi_vertex[right_ctx_idxs] + 32492
            
            # average parcel-wise fraction of cell types (DFC)
            cell_types = gsub('DFC_', '', cell_types)
            
            # do it
            DFC_schaef_out = parcel_average(ciber_dat=sample_dat_ciber, ref_region='DFC', parcels=parcels, net=net, cell_types=cell_types)
            DFC_donor_specific_expression = DFC_schaef_out[[2]]
            DFC_schaef_out = DFC_schaef_out[[1]]
            
            # DFC - create averaged plots
            lake_dfc_dat = read_csv(paste0(base_dir, '/data/ahba/schaeffer_LAKE_DFC_',parcels,'_',net,'Net_expr_mat_new_',para,'.csv'))
            
            # Spatial correlation of DFC cell types
            dfc_sc_spatial_cors = cor(t(lake_dfc_dat[,1:parcel_num]), use = 'pairwise.complete')
            colnames(dfc_sc_spatial_cors) = lake_dfc_dat$gene
            rownames(dfc_sc_spatial_cors) = lake_dfc_dat$gene
            
            
            # get ready for plotting
            my_palette = colorRampPalette(c('#67001F','#B2182B','#D6604D','#F4A582','#FDDBC7','#F7F7F7','#D1E5F0','#92C5DE','#4393C3','#2166AC','#053061'))(n = 299)
            my_palette = rev(my_palette)
            cell_text  = t(as.matrix(dfc_sc_spatial_cors))
            cell_text  = round(cell_text,2)
            
            # plot corrmat!
            pdf(file = paste0(base_dir,'/figures/PaperPlot_lake_dfc_spatial_cors_',para,'.pdf'), width=10, height=10)
            heatmap.2(as.matrix(dfc_sc_spatial_cors), cellnote=cell_text, trace="none", notecol='black',col=my_palette,
                      distfun=function(x) dist(x, method="euclidean"),dendrogram="col",density.info="none",
                      hclustfun=function(x) hclust(x, method="ward.D2"),Rowv=T, symm = TRUE)
            dev.off()
        }
    }
}


######### Visual ############
region = 'VisualCortex'
# for different parameter combinations
#for ( para in c('0.1','0.3','0.5','ProbeMax', 'NormZscore', 'NormZscore0.3') ){
# for the finalized parameter combo
for ( para in c('NormZscore0.3') ){
    for (parcels in c('100','200','300','400','500','600','700','800','900','1000')){
        parcel_num = strtoi(parcels)
        for (net in c('7','17')){
            # read cell type expression of "signature genes"
            gep_file = paste0(ciber_dir, '/CIBERSORTx_Lake_',region,'_ahba_matched_sc_SignatureMatrix_new_',para,'.txt')
            gep = read.table(gep_file, header=T)
            
            # create heatmap of gene signatures
            file_out = paste0(paste0(base_dir, '/figures/',region,'_gep_signature_matrix_new_',para,'.pdf'))
            CairoPDF(file_out, height=5, width=4)
            heatmap.2(as.matrix(gep[,2:ncol(gep)]), col=brewer.pal(11,"RdBu"), key=FALSE, dendrogram='col', scale="row", trace="none", labRow=FALSE)
            dev.off()
            file_out
            
            # average expression of signature genes in each cell class
            sig_cormat = cor(gep[,2:ncol(gep)], use = 'pairwise.complete')
            cell_text = t(as.matrix(sig_cormat))
            cell_text = round(cell_text,2)
            # color pallette for the correlation plot
            my_palette = colorRampPalette(c('#67001F','#B2182B','#D6604D','#F4A582','#FDDBC7','#F7F7F7','#D1E5F0','#92C5DE','#4393C3','#2166AC','#053061'))(n = 299)
            my_palette = rev(my_palette)
            
            # plot
            pdf(file = paste0(base_dir, '/figures/lake_',region,'_corr_of_signature_mat_new_',para,'.pdf'), width=10, height=10)
            heatmap.2(as.matrix(sig_cormat), cellnote=cell_text, trace="none", notecol='black',col=my_palette,breaks=seq(-1, 1, length.out=300),
                      distfun=function(x) dist(x, method="euclidean"),dendrogram="col",density.info="none",
                      hclustfun=function(x) hclust(x, method="ward.D2"),Rowv=T, symm=F,symkey=F,symbreaks=T)
            dev.off()
            
            
            # read AHBA informatioN
            load(file=paste0(base_dir, '/data/ahba/ahba_ctx_norm_abagen_',para,'.Rdata'), verbose=T)
            donor_arr = c("9861","10021","12876","14380","15496","15697")
            
            # read sample-to-vertex projection info
            sample_info = read_csv(paste0(base_dir, '/data/ahba/sample_info_vertex_reannot_mapped_',para,'.csv'))
            sample_dat  = sample_info[which(abs(sample_info$mm_to_surf) < 4),]
            
            write(region,'')
            all_ciber = NULL
            for (donor in donor_arr){
                write(donor, '')
                
                path_to_ciber = paste0(ciber_dir, '/',region,'_ahba_',donor,'_fraction_new_',para,'.txt')
                ciber_dat = read_delim(path_to_ciber, delim='\t') #read.table(path_to_ciber, header=T)
                all_ciber = rbind(all_ciber, ciber_dat)
                
                ciber_field = paste0(region, '_CiberFrac')
                ahba_data[[donor]][[ciber_field]] = merge(x=ciber_dat, y=sample_info, by.y='well_id', by.x='Mixture')
            }
            write_csv(x=all_ciber, path=paste0(ciber_dir, '/Cibersortx_AHBAsamples_Lake_',region,'_new_',para,'.csv'))
            
            # Read cell fractions derived from Visual Cortex cell signatures
            ciber = read_csv(paste0(ciber_dir, '/Cibersortx_AHBAsamples_Lake_',region,'_new_',para,'.csv'))
            colnames(ciber) = paste0('VIS_',colnames(ciber))
            
            # merge cibersort data with AHBA sample dat
            cell_types = colnames(ciber)[2:19]
            sample_dat_ciber = merge(x=sample_dat, y=ciber, by.x='well_id', by.y='VIS_Mixture', all.x=T)
            
            # Collapse individual cell proportions into schaeffer parcels
            # stitch the left/right verticies together to match Schaeffer parcel cifti format
            sample_dat_ciber$bihemi_vertex = sample_dat_ciber$vertex + 1 # cifti indices index at 0, R indexes at 1
            right_ctx_idxs = intersect(grep('right', sample_dat_ciber$structure_name), which(sample_dat_ciber$top_level == 'CTX'))
            sample_dat_ciber$bihemi_vertex[right_ctx_idxs] = sample_dat_ciber$bihemi_vertex[right_ctx_idxs] + 32492
            
            # average parcel-wise fraction of cell types (VIS)
            cell_types = gsub('VIS_', '', cell_types)
            
            # do it
            VIS_schaef_out = parcel_average(ciber_dat=sample_dat_ciber, ref_region='VIS', parcels=parcels, net=net, cell_types=cell_types)
            VIS_donor_specific_expression = VIS_schaef_out[[2]]
            VIS_schaef_out = VIS_schaef_out[[1]]
            
            # VIS - create averaged plots
            lake_vis_dat = read_csv(paste0(base_dir, '/data/ahba/schaeffer_LAKE_VIS_',parcels,'_',net,'Net_expr_mat_new_',para,'.csv'))
            
            
            # Spatial correlation of DFC cell types
            vis_sc_spatial_cors = cor(t(lake_vis_dat[,1:parcel_num]), use = 'pairwise.complete')
            colnames(vis_sc_spatial_cors) = lake_vis_dat$gene
            rownames(vis_sc_spatial_cors) = lake_vis_dat$gene
            
            
            # get ready for plotting
            my_palette = colorRampPalette(c('#67001F','#B2182B','#D6604D','#F4A582','#FDDBC7','#F7F7F7','#D1E5F0','#92C5DE','#4393C3','#2166AC','#053061'))(n = 299)
            my_palette = rev(my_palette)
            cell_text  = t(as.matrix(vis_sc_spatial_cors))
            cell_text  = round(cell_text,2)
            
            # plot corrmat!
            pdf(file = paste0(base_dir,'/figures/PaperPlot_lake_vis_spatial_cors_',para,'.pdf'), width=10, height=10)
            heatmap.2(as.matrix(vis_sc_spatial_cors), cellnote=cell_text, trace="none", notecol='black',col=my_palette,
                      distfun=function(x) dist(x, method="euclidean"),dendrogram="col",density.info="none",
                      hclustfun=function(x) hclust(x, method="ward.D2"),Rowv=T, symm = TRUE)
            dev.off()
        }
    }
}
