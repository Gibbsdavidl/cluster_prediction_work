###
### This file is for finding the intersection of gene symbols across all data sets we're considering

library(readr)
library(stringr)

## EDIT THESE LOCATIONS##

### WHERE THE DATA DIR LIVES ###
data_dir <- 'D:/Work/cluster_prediction_data/'

### WHERE TO WRITE THE NEW TABLES ###
out_dir <- paste0(data_dir, 'data/')


### The Pan-Cancer EB++ (named EBpp) batch corrected data set #########################################################
#
# read in the table of gene symbols per row
ebpp_genes <- readr::read_csv(paste0(data_dir,'data/EBpp/ebpp_genes.txt'))
gene_symbols <- unlist(sapply(
  sapply(ebpp_genes, function(x) str_split(x,pattern='\\|')),
  function(y) y[[1]]
))


### RSubRead ##########################################################################################################
#
###
Sys.setenv(VROOM_CONNECTION_SIZE=500072)
subread <- read_delim(paste0(data_dir,'data/Rsubread/GSM1536837_01_27_15_TCGA_20.Illumina.tumor_Rsubread_TPM.tsv.gz'), delim = '\t')

#extract the genes
subread_genes <- as.vector(subread[,1])[[1]]


### Xena Kallisto TPM #################################################################################################
#
###
gene_map_tpm <- read_csv(paste0(data_dir, 'data/xena_kallisto_tpm/gene_id_map.csv'))


### Xena RSEM FPKM #################################################################################################
#
###
gene_map_fpkm <- read_csv(paste0(data_dir, 'data/xena_rsem_fpkm/gene_id_map.csv'))
#gene_dup <- which( (!duplicated(gene_map$gene)) & (!is.na(gene_map$gene)) )


collected_genes <- intersect(gene_symbols, subread_genes)
collected_genes <- intersect(collected_genes, gene_map_tpm$gene)
collected_genes <- intersect(collected_genes, gene_map_fpkm$gene)

gene_df <- data.frame(Genes=collected_genes)

write.table(gene_df, paste0(out_dir,'available_genes.csv'), sep = ',', quote=F, row.names = F)
