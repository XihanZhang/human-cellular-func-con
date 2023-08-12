# Load packages
library(tidyverse)
library(WGCNA)

# this script extract and aggregate necessary information from tables downloaded from AHBA (http://human.brain-map.org/)
# this script was adapted from Anderson et al 2020 NatCom 
# (https://github.com/HolmesLab/2020_NatComm_interneurons_cortical_function_schizophrenia/blob/master/scripts/01_preprocess_ahba.R)

# this function will map AHBA samples to their overarching regional category (e.g. mPFC >> CORTEX; mediodorsal_thal >> THAL)
find_top_level = function(ref_df, struct_id, ontology, top_level) {

    ontology_row  = which(ontology$id == struct_id )

    # print feedback if more than one match, this should never be the case
    if (length(ontology_row) > 1){
        write("ERROR")
    }

    ontology_info = ontology[ontology_row,]
    splits    = strsplit(as.character(ontology_info$structure_id_path), '/')[[1]]
    reg.match = ref_df$id[which(ref_df$id %in% splits)]
    if (length(reg.match) == 1){
        region_cat   = as.character(ref_df$top_level[which(ref_df$id %in% splits)])
        region_name  = ref_df[ref_df$id == reg.match,]$name
        return(c(reg.match, region_cat, region_name))
    } else {
        return(c(NA,NA,NA))
    }
}


# Modify these filepaths for your local directory structure
# Yale cluster milgram dir:
#base_dir = '/gpfs/milgram/project/holmes/xz555/gradient_shift'
# laptop dir:
base_dir = '/Volumes/GoogleDrive/My Drive/Gradient_Shift_Cellular_Basis/Milgram'

# path to downloaded raw AHBA microarray data
data_path  = paste(base_dir, '/data/ahba', sep = '') 
donor_nums = c('9861', '10021', '12876', '14380', '15496', '15697')


# Map ontology structure names to more general sample categories
# e.g. left motor gyrus --> CTX
top_level = read_csv(paste0(base_dir, '/ref_files/top_level_categories.csv'))


# Read data from each subject
sampleInfo = NULL
for ( donor in donor_nums ) {
    file = paste('donor', donor, sep='')
    write(paste('Reading and collapsing data for donor: ', donor, sep = ''),'')

    # read sample Information
    saname    = paste(data_path, file, 'SampleAnnot.csv', sep='/')
    samp_info = read_csv(saname)
    samp_info$brain = donor

    # concatenate sample info into one overall matrix
    sampleInfo = rbind(sampleInfo, samp_info)

    # Gene Ontology Info
    oname    = paste(data_path, file, 'Ontology.csv', sep='/')
    ont_data = read_csv(oname)
    ontology = ont_data

    # identify samples in the regions we want to analyze
    out         = do.call(rbind, lapply(samp_info$structure_id, find_top_level, ref_df=top_level, ontology=ontology))
    reg_info    = as.data.frame(out)
    reg_info$V1 = as.numeric(as.character(reg_info$V1))
    colnames(reg_info) = c('reg_num', 'top_level', 'region_clean')

}

# map each ahba subregion to its overall category
out         = do.call(rbind, lapply(sampleInfo$structure_id, find_top_level, ref_df=top_level, ontology=ontology))
reg_info    = as.data.frame(out)
reg_info$V1 = as.numeric(as.character(reg_info$V1))
colnames(reg_info) = c('reg_num', 'top_level', 'region_clean')


# add new sample super-category info to original sample information dataframe
sampleInfo    = cbind(sampleInfo, reg_info)

# output the csv table
write_csv(sampleInfo, path=paste0(base_dir, '/para/data/ahba/ahba_sampleInfo.csv'))