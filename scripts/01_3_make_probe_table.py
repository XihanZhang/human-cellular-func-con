#!/bin/python

import pandas as pd
import abagen
import os.path

# set up dirs
base_dir = '/Volumes/GoogleDrive/My Drive/Gradient_Shift_Cellular_Basis/Milgram'
data_dir = os.path.join(base_dir,'para','data','ahba')
micro_path = os.path.join(base_dir,'microarray','normalized_microarray_donor9861','Probes.csv')

# loop through gene_entrezid derived from different abagen parameters
paras = ['0.1', '0.3', '0.5', 'ProbeMax', 'NormZscore']

for para in paras:
    # load the converted entrez ID
    gene_symbols_entrezID_name_chromosome = pd.read_csv(os.path.join(data_dir,f'ahba_expression_EntrezID_GeneName_Chromosome_{para}.csv'))
    gene_symbols_entrezID_name_chromosome = gene_symbols_entrezID_name_chromosome.rename(columns={'Unnamed: 0': 'gene_symbol'})
    # find the NA entries
    null_index = gene_symbols_entrezID_name_chromosome[gene_symbols_entrezID_name_chromosome['entrez_id'].isnull()].index.tolist() # 515 missing
    NA_gene_symbols = gene_symbols_entrezID_name_chromosome['gene_symbol'][null_index].tolist()

    # load full list of gene symbols
    expression = pd.read_csv(os.path.join(data_dir,f'ahba_group_samples_expression_{para}.csv'))
    gene_symbols = expression.columns.tolist()

    # load the probe info provided by ahba
    probes = pd.read_csv(micro_path)

    # drop duplicated rows (as one gene has more than one probes)
    probes_sub = probes.drop_duplicates(subset=['entrez_id'], keep='first')
    # drop the useless columns
    probes_sub = probes_sub.drop(['probe_id','probe_name'], axis=1)
    # keep probes by gene_symbols
    probes_sub = probes_sub.loc[probes_sub["gene_symbol"].isin(gene_symbols)]
    # find the entrez ID of NA
    probes_sub_has_NA_entrezID = probes_sub.loc[probes_sub["gene_symbol"].isin(NA_gene_symbols)] # 343 filling
    probes_sub_has_NA_entrezID.reset_index(inplace=True)
    probes_sub_has_NA_entrezID

    for i in range(len(probes_sub_has_NA_entrezID)):
        [this_gene_index] = gene_symbols_entrezID_name_chromosome.index[gene_symbols_entrezID_name_chromosome['gene_symbol'] == probes_sub_has_NA_entrezID["gene_symbol"][i]].tolist()
        gene_symbols_entrezID_name_chromosome['entrez_id'][this_gene_index] = probes_sub_has_NA_entrezID["entrez_id"][i]
        gene_symbols_entrezID_name_chromosome['gene_name'][this_gene_index] = probes_sub_has_NA_entrezID["gene_name"][i]
        gene_symbols_entrezID_name_chromosome['chromosome'][this_gene_index] = probes_sub_has_NA_entrezID["chromosome"][i]

    # output the gene symbol and entrez id pairs
    gene_symbols_entrezID_name_chromosome.to_csv(os.path.join(data_dir,f'ahba_expression_EntrezID_GeneName_Chromosome_full_{para}.csv'),index=False)
