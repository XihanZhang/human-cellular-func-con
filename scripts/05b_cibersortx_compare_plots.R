library(tidyverse)
library(gplots)
library(RColorBrewer)
library(Cairo)
library(cifti)
library(cowplot)

# function to plot parcel-wise cortical correlations for viewing with HCP wb_view
plot_matlab = function(values, out_path, parcel_num, net_num){
    
    base_dir         = '/Volumes/GoogleDrive/My Drive/Gradient_Shift_Cellular_Basis/Milgram'
    dscalar_template = paste0(base_dir, '/data/Schaefer/Schaefer2018_',parcel_num,'Parcels_',net_num,'Networks_order.dscalar.nii')
    parcel_info_file = paste0(base_dir, '/data/Schaefer/Schaefer2018_',parcel_num,'Parcels_',net_num,'Networks_order_info.txt')
    write_val_file   = paste0('/Volumes/GoogleDrive/My Drive/Gradient_Shift_Cellular_Basis/Milgram/data/ahba/tmp/tmp.txt')
    
    write_delim(delim=' ', x=as.data.frame(values), path=write_val_file,col_names=F)
    save_path = out_path
    print(save_path)
    
    matfunc = 'plotVolOnSurface'
    cmd = paste0('/Applications/MATLAB_R2017b.app/bin/matlab -nosplash -nodesktop -r "cd(\'/Volumes/GoogleDrive/My Drive/Gradient_Shift_Cellular_Basis/Milgram/scripts\');',
                 matfunc, '(\'', dscalar_template, '\',\'', parcel_info_file, '\',\'',
                 write_val_file, '\',\'', save_path, '\'); exit;"')
    system(cmd)
}

# set up directories
base_dir  = '/Volumes/GoogleDrive/My Drive/Gradient_Shift_Cellular_Basis/Milgram/para'
old_base_dir  = '/Volumes/GoogleDrive/My Drive/Gradient_Shift_Cellular_Basis/Milgram'
ciber_dir = paste0(base_dir, '/data/cibersortX')
region = 'FrontalCortex'
load(file='/Volumes/GoogleDrive/My Drive/Gradient_Shift_Cellular_Basis/Milgram/data/ahba/ahba_data_object.Rdata', verbose=T)

for ( para in c('0.1','0.3','0.5','ProbeMax', 'NormZscore', 'NormZscore0.3','original') ){
    
    # 1. Cell-type signature matrix & its correlation matrix
    region = 'VisualCortex'
    # cell-type signature matrix
    gep_file = paste0(ciber_dir, '/CIBERSORTx_Lake_',region,'_ahba_matched_sc_SignatureMatrix_new_',para,'.txt')
    gep = read.table(gep_file, header=T)
    
    # average expression of signature genes in each cell class
    sig_cormat = cor(gep[,2:ncol(gep)], use = 'pairwise.complete')
    cell_text = t(as.matrix(sig_cormat))
    cell_text = round(cell_text,2)
    # color pallette for the correlation plot
    my_palette = colorRampPalette(c('#67001F','#B2182B','#D6604D','#F4A582','#FDDBC7','#F7F7F7','#D1E5F0','#92C5DE','#4393C3','#2166AC','#053061'))(n = 299)
    my_palette = rev(my_palette)
    
    # plot
    pdf(file = paste0(base_dir, '/figures/lake_',region,'_corr_of_signature_mat_NoCluster_new_',para,'.pdf'), width=10, height=10)
    heatmap.2(as.matrix(sig_cormat), cellnote=cell_text,, trace="none", notecol='black',col=my_palette,Colv = NA, Rowv = NA,breaks=seq(-1, 1, length.out=300))
    dev.off()
    
    
    region = 'FrontalCortex'
    # cell-type signature matrix
    gep_file = paste0(ciber_dir, '/CIBERSORTx_Lake_',region,'_ahba_matched_sc_SignatureMatrix_new_',para,'.txt')
    gep = read.table(gep_file, header=T)
    
    # average expression of signature genes in each cell class
    sig_cormat = cor(gep[,2:ncol(gep)], use = 'pairwise.complete')
    cell_text = t(as.matrix(sig_cormat))
    cell_text = round(cell_text,2)
    # color pallette for the correlation plot
    my_palette = colorRampPalette(c('#67001F','#B2182B','#D6604D','#F4A582','#FDDBC7','#F7F7F7','#D1E5F0','#92C5DE','#4393C3','#2166AC','#053061'))(n = 299)
    my_palette = rev(my_palette)
    
    # plot
    pdf(file = paste0(base_dir, '/figures/lake_',region,'_corr_of_signature_mat_NoCluster_new_',para,'.pdf'), width=10, height=10)
    heatmap.2(as.matrix(sig_cormat), cellnote=cell_text,, trace="none", notecol='black',col=my_palette,Colv = NA, Rowv = NA,breaks=seq(-1, 1, length.out=300))
    dev.off()
    
    
    # 2. Cell-type signature matrix 
    # DFC - create averaged plots
    nparcel = '400'
    net_num = '7'
    lake_dfc_dat = read_csv(paste0(base_dir, '/data/ahba/schaeffer_LAKE_DFC_',nparcel,'_',net_num,'Net_expr_mat_new_',para,'.csv'))
    #for (cell in lake_dfc_dat$gene){
        #write(cell,'')
        #cell_proportion = as.numeric(lake_dfc_dat[lake_dfc_dat$gene == cell,1:as.numeric(nparcel)])
        #cell_proportion[which(cell_proportion == 0)] = .001
        #plot_matlab(values=as.numeric(cell_proportion), out_path=paste0(base_dir, '/figures/surface_plots/PaperFig_Lake_DFC_',cell,'_schaef',nparcel,'_net',net_num,'_new_',para,'.dscalar.nii'), parcel_num=nparcel, net_num=net_num)
    #}
    #sst_proportion   = as.numeric(lake_dfc_dat[lake_dfc_dat$gene == 'DFC_SST',1:as.numeric(nparcel)])
    #pvalb_proportion = as.numeric(lake_dfc_dat[lake_dfc_dat$gene == 'DFC_PVALB',1:as.numeric(nparcel)])
    #sst_minus_pvalb  = scale(sst_proportion) - scale(pvalb_proportion)
    #sst_minus_pvalb[which(sst_minus_pvalb == 0)] = .001
    #plot_matlab(values=as.numeric(sst_minus_pvalb), out_path=paste0(base_dir, '/figures/surface_plots/PaperFig_schaeffer_LAKE_DFC_SST-PVALB_schaef_',nparcel,'_net',net_num,'_new_',para,'.dscalar.nii'), parcel_num=nparcel, net_num=net_num)
    
    # Spatial correlation of DFC cell types
    dfc_sc_spatial_cors = cor(t(lake_dfc_dat[,1:400]), use = 'pairwise.complete')
    colnames(dfc_sc_spatial_cors) = lake_dfc_dat$gene
    rownames(dfc_sc_spatial_cors) = lake_dfc_dat$gene
    
    
    # get ready for plotting
    my_palette = colorRampPalette(c('#67001F','#B2182B','#D6604D','#F4A582','#FDDBC7','#F7F7F7','#D1E5F0','#92C5DE','#4393C3','#2166AC','#053061'))(n = 299)
    my_palette = rev(my_palette)
    cell_text  = t(as.matrix(dfc_sc_spatial_cors))
    cell_text  = round(cell_text,2)
    
    # plot corrmat!
    pdf(file = paste0(base_dir,'/figures/PaperPlot_lake_dfc_spatial_cors_noclustering_',para,'.pdf'), width=10, height=10)
    heatmap.2(as.matrix(dfc_sc_spatial_cors), cellnote=cell_text,, trace="none", notecol='black',col=my_palette,Colv = NA, Rowv = NA)
    dev.off()
}







