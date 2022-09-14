

library(readr)
library(stringr)

### WHERE THE DATA DIR LIVES ###
data_dir <- 'D:/Work/cluster_prediction_data/'
#data_dir <- 'C:/Users/dgibbs/ISB_Work/cluster_prediction_data/'

### WHERE TO WRITE THE NEW TABLES ###
out_dir <- paste0(data_dir, 'formatted/')

# read in the table of gene symbols per row
available_genes <- readr::read_csv(paste0(data_dir,'data/available_genes.csv'))

# first we'll select some genes:
# might need to remove '?s'
random_genes <- sample(available_genes$Genes,size=100)

# then create some gene signatures:
random_sigs <- list()
for (i in 1:5){
  random_sigs[[i]] <- sample(available_genes$Genes,size=25)
  random_genes <- c(random_genes, random_sigs[[i]])
}

gene_table1 <- data.frame(Type='Feature', Gene=random_genes)
gene_table2 <- data.frame(Type='Set1', Gene=random_sigs[[1]])
gene_table3 <- data.frame(Type='Set2', Gene=random_sigs[[2]])
gene_table4 <- data.frame(Type='Set3', Gene=random_sigs[[3]])

gene_table <- rbind(gene_table1, gene_table2, gene_table3, gene_table4)

idx <- which( (gene_table$Gene == '?') )

if (length(idx) > 0) {
  gene_table <- gene_table[-idx,]
}

write.csv(gene_table, paste0(out_dir,'random_genes_table.csv'), row.names = FALSE, quote = F)
