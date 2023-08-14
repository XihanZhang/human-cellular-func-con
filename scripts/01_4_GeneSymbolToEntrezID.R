# this script is to generate entez ID for the gene symbol
# install
#BiocManager::install("org.Hs.eg.db")

# load the package
library('org.Hs.eg.db')
library(dplyr)

# see the available identifiers that can be used in this package
#columns(org.Hs.eg.db)
#keytypes(org.Hs.eg.db)
#head(keys(org.Hs.eg.db, keytype="SYMBOL"))
#head(keys(org.Hs.eg.db, keytype="GENENAME"))

# set dirs
base_dir  = '/Volumes/GoogleDrive/My Drive/Gradient_Shift_Cellular_Basis/Milgram/para'
data_dir = paste0(base_dir, '/data/ahba')

# loop through the genesymbol files derived from different abagen parameters
for ( para in c('0.1', '0.3', '0.5', 'ProbeMax', 'NormZscore', 'NormZscore0.3') ){
  # load the list of gene symbols
  gene_symbol_file = paste0(data_dir, '/ahba_expression_gene_names_', para, '.csv')
  symbols <- read.csv(gene_symbol_file,header=F,stringsAsFactors=F)
  symbols <- unlist(symbols)
  
  # 1. use mapIds method to obtain Entrez IDs
  gene_entrez_ids <- mapIds(org.Hs.eg.db, symbols, 'ENTREZID', 'SYMBOL')
  
  # 2. use mapIds method to obtain gene name
  gene_name <- mapIds(org.Hs.eg.db, symbols, 'GENENAME', 'SYMBOL')
  
  # merge the Entrez IDs and gene name into one table
  full_table <- data.frame(gene_entrez_ids)
  colnames(full_table) <- 'entrez_id'
  full_table$gene_name = gene_name
  
  # 3. make table of corresponding chromosomes
  # (1) create a table of gene symbol-chromosome map
  x <- org.Hs.egCHR
  mapped_genes <- mappedkeys(x)
  xx <- as.list(x[mapped_genes])
  
  table_entrezID_chromosome <- data.frame(names(xx))
  colnames(table_entrezID_chromosome) <- 'entrez_id'
  
  table_chromosome <- list()
  for(id in 1:length(table_entrezID_chromosome$entrez_id)) {
    this = table_entrezID_chromosome$entrez_id[id]
    table_chromosome <- c(table_chromosome,xx[[this]][1])
  }
  
  table_entrezID_chromosome$chromosome = table_chromosome
  
  # (2) fill in the chromosome by entrez id
  for(id in 1:length(full_table$entrez_id)) {
    this = full_table$entrez_id[id]
    if (is.na(this)) {
      full_table$chromosome[id] <-'NA'
    } else {
      full_table$chromosome[id] <- unlist(table_entrezID_chromosome$chromosome[table_entrezID_chromosome$entrez_id==this])
    }
  }
  
  chromosome = full_table$chromosome
  chromosome = unlist(chromosome)
  full_table$chromosome = chromosome
  
  # output the Full table
  gene_entrez_file = paste0(data_dir, '/ahba_expression_EntrezID_GeneName_Chromosome_', para, '.csv')
  write.csv(full_table, gene_entrez_file, row.names = TRUE)
  
}









