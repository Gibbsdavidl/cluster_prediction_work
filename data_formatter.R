
library(readr)
library(stringr)

## EDIT THESE LOCATIONS##

### WHERE THE DATA DIR LIVES ###
data_dir <- 'F:/Work/cluster_prediction_data/'
#data_dir <- 'C:/Users/dgibbs/ISB_Work/cluster_prediction_data/'

### WHERE TO WRITE THE NEW TABLES ###
out_dir <- paste0(data_dir, 'results/')

### GENES TO INCLUDE AS A FILE ###
gene_table <- read.csv(paste0(data_dir,'formatted/random_genes_table.csv'), header=T, stringsAsFactors = F)
all_genes <- unique(gene_table$Gene)

### smaller size for testing ###
subset_subtypes <- 200

## END EDITING


### The immune subtypes ###
#   Use: ClusterModel1 as the subtype
###
subtypes <- readr::read_delim(paste0(data_dir,'data/five_signature_mclust_ensemble_results.tsv.gz'), delim = '\t')
subtypes$AliquotBarcode <- str_replace_all(string = subtypes$AliquotBarcode, pattern = '\\.', replacement = '-')
subtypes$SampleBarcode2 <- str_sub(subtypes$SampleBarcode, start = 1, end = 15)

if (subset_subtypes > 0){
  random_barcodes <- sample(subtypes$AliquotBarcode, size = subset_subtypes, )
  subtypes <- subtypes[subtypes$AliquotBarcode %in% random_barcodes,]
}

### The Pan-Cancer EB++ (named EBpp) batch corrected data set #########################################################
#
###
ebpp <- read_delim(paste0(data_dir,'data/EBpp/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv'), delim = '\t')

# drop the gene ID column
ebpp$gene_id <- NULL
ebpp <- as.data.frame(ebpp)

# read in the table of gene symbols per row
ebpp_genes <- readr::read_csv(paste0(data_dir,'data/EBpp/ebpp_genes.txt'))
gene_symbols <- unlist(sapply(
  sapply(ebpp_genes, function(x) str_split(x,pattern='\\|')),
  function(y) y[[1]]
))

gene_dup <- which(!duplicated(gene_symbols))
ebpp <- ebpp[gene_dup,]
rownames(ebpp) <- gene_symbols[gene_dup]

# subset data
# and get it into order
ebpp_sub <- ebpp[all_genes,subtypes$AliquotBarcode]

# we transpose it and add in the label
ebpp_sub_t <- as.data.frame(t(as.data.frame(ebpp_sub)))
ebpp_sub_t[is.na(ebpp_sub_t)] <- 0.0

# add in the labels and barcodes
ebpp_sub_t[['Label']] <- subtypes$ClusterModel1
ebpp_sub_t[['Barcode']] <- subtypes$AliquotBarcode

# and write it out for training
write.csv(ebpp_sub_t, paste0(data_dir,'formatted/EBPP_subset.csv'), row.names = FALSE, quote = F)



### RSubRead ##########################################################################################################
#
###
Sys.setenv(VROOM_CONNECTION_SIZE=500072)
subread <- read_delim(paste0(data_dir,'data/Rsubread/GSM1536837_01_27_15_TCGA_20.Illumina.tumor_Rsubread_TPM.tsv.gz'), delim = '\t')

#extract the genes
subread_genes <- as.vector(subread[,1])[[1]]

# remove the gene ID column
subread <- subread[, -c(1)]
subread <- as.data.frame(subread)

gene_dup <- which(!duplicated(subread_genes))
subread <- subread[gene_dup,]
rownames(subread) <- subread_genes[gene_dup]

# subset data
# and get it into order
present_barcodes <- intersect(subtypes$AliquotBarcode, colnames(subread))
subtypes_subread <- subtypes[subtypes$AliquotBarcode %in% present_barcodes, ]
subread_sub <- subread[all_genes,subtypes_subread$AliquotBarcode]

# we transpose it and add in the label
subread_sub_t <- as.data.frame(t(as.data.frame(subread_sub)))
subread_sub_t[is.na(subread_sub_t)] <- 0.0

# add in the labels and barcodes
subread_sub_t[['Label']] <- subtypes_subread$ClusterModel1
subread_sub_t[['Barcode']] <- subtypes_subread$AliquotBarcode

# and write it out for training
write.csv(subread_sub_t, paste0(data_dir,'formatted/Rsubread_subset.csv'), row.names = FALSE, quote = F)



### Xena Kallisto TPM #################################################################################################
#
###
xenak <- read_delim(paste0(data_dir,'data/xena_kallisto_tpm/tcga_Kallisto_tpm.tsv.gz'), delim = '\t')
#gene_map <- read_delim(paste0(data_dir, 'data/annot/ensembl_GRCh37_p13_gene_ID_map.txt.gz'), delim ='\t')
gene_map <- read_csv(paste0(data_dir, 'data/xena_kallisto_tpm/gene_id_map.csv'))

# mapping the gene names
xenak_genes <- xenak[,1]

# remove the gene ID column
xenak <- xenak[, -c(1)]
xenak <- as.data.frame(xenak)

gene_dup <- which( (!duplicated(gene_map$gene)) & (!is.na(gene_map$gene)) )
xenak <- xenak[gene_dup,]
rownames(xenak) <- gene_map$gene[gene_dup]

# subset data and get it into order
present_barcodes <- intersect(subtypes$SampleBarcode2, colnames(xenak))
subtypes_xenak <- subtypes[subtypes$SampleBarcode2 %in% present_barcodes, ]
xenak_sub <- xenak[all_genes,subtypes_xenak$SampleBarcode2]

# we transpose it and replace NAs with 0s
xenak_sub_t <- as.data.frame(t(as.data.frame(xenak_sub)))
xenak_sub_t[is.na(xenak_sub_t)] <- 0.0

# add in the labels and barcodes
xenak_sub_t[['Label']] <- subtypes_xenak$ClusterModel1
xenak_sub_t[['Barcode']] <- subtypes_xenak$AliquotBarcode

# and write it out for training
write.csv(xenak_sub_t, paste0(data_dir,'formatted/xena_kallisto_subset.csv'), row.names = FALSE, quote = F)


### Xena RSEM FPKM #################################################################################################
#
###
xenak <- read_delim(paste0(data_dir,'data/xena_rsem_fpkm/tcga_RSEM_gene_fpkm.tsv.gz'), delim = '\t')
#gene_map <- read_delim(paste0(data_dir, 'data/annot/ensembl_GRCh37_p13_gene_ID_map.txt.gz'), delim ='\t')
gene_map <- read_csv(paste0(data_dir, 'data/xena_rsem_fpkm/gene_id_map.csv'))

# mapping the gene names
xenak_genes <- xenak[,1]

# remove the gene ID column
xenak <- xenak[, -c(1)]
xenak <- as.data.frame(xenak)

gene_dup <- which( (!duplicated(gene_map$gene)) & (!is.na(gene_map$gene)) )
xenak <- xenak[gene_dup,]
rownames(xenak) <- gene_map$gene[gene_dup]

# subset data and get it into order
present_barcodes <- intersect(subtypes$SampleBarcode2, colnames(xenak))
subtypes_xenak <- subtypes[subtypes$SampleBarcode2 %in% present_barcodes, ]
xenak_sub <- xenak[all_genes,subtypes_xenak$SampleBarcode2]

# we transpose it and replace NAs with 0s
xenak_sub_t <- as.data.frame(t(as.data.frame(xenak_sub)))
xenak_sub_t[is.na(xenak_sub_t)] <- 0.0

# add in the labels and barcodes
xenak_sub_t[['Label']] <- subtypes_xenak$ClusterModel1
xenak_sub_t[['Barcode']] <- subtypes_xenak$AliquotBarcode

# and write it out for training
write.csv(xenak_sub_t, paste0(data_dir,'formatted/xena_rsem_fpkm_subset.csv'), row.names = FALSE, quote = F)


### Xena RSEM TPM #################################################################################################
#
###
