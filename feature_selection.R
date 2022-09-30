
library(readr)
library(stringr)

### Selecting features ###

### WHERE THE DATA DIR LIVES ###
data_dir <- 'C:/Users/dgibbs/ISB_Work/cluster_prediction_data/'

### WHERE TO WRITE THE NEW TABLES ###
out_dir <- paste0(data_dir, 'results/')

### The immune subtypes ###
#   Use: ClusterModel1 as the subtype
###
subtypes <- readr::read_delim(paste0(data_dir,'data/five_signature_mclust_ensemble_results.tsv.gz'), delim = '\t')
subtypes$AliquotBarcode <- str_replace_all(string = subtypes$AliquotBarcode, pattern = '\\.', replacement = '-')
subtypes$SampleBarcode2 <- str_sub(subtypes$SampleBarcode, start = 1, end = 15)


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
ebpp_sub <- ebpp[,subtypes$AliquotBarcode]
ebpp_sub <- ebpp_sub[2:20502,]

# we transpose it and add in the label
ebpp_sub_t <- as.data.frame(t(as.data.frame(ebpp_sub)))
#ebpp_sub_t[is.na(ebpp_sub_t)] <- 0.0

rm(ebpp)

data <- list(subtypes, ebpp_sub_t)

save(data, file = paste0(data_dir,'feature_selection_data.rda'))


# clear memory
rm(list=ls())

### WHERE THE DATA DIR LIVES ###
data_dir <- 'C:/Users/dgibbs/ISB_Work/cluster_prediction_data/'

### WHERE TO WRITE THE NEW TABLES ###
out_dir <- paste0(data_dir, 'results/')

# load the data
load(paste0(data_dir,'feature_selection_data.rda'))

kt_fun <- function(x,subtypes,clusterval){
  try({
    y <- ifelse(subtypes$ClusterModel1 == clusterval,1, 0)
    res0 <- kruskal.test(x,y)
    return(res0$statistic)
  })
  return(0.0)
}

meddiff <- function(x,subtypes,clusterval) {
  y <- ifelse(subtypes$ClusterModel1 == clusterval,1, 0)
  return(median(as.numeric(x[y==1]),na.rm=TRUE) - median(as.numeric(x[y==0]),na.rm=TRUE))
}

c1meddiff <- as.numeric(apply(ebpp, 2, function(a) meddiff(a,subtypes,1)))
c2meddiff <- as.numeric(apply(ebpp, 2, function(a) meddiff(a,subtypes,2)))
c3meddiff <- as.numeric(apply(ebpp, 2, function(a) meddiff(a,subtypes,3)))
c4meddiff <- as.numeric(apply(ebpp, 2, function(a) meddiff(a,subtypes,4)))
c5meddiff <- as.numeric(apply(ebpp, 2, function(a) meddiff(a,subtypes,5)))
c6meddiff <- as.numeric(apply(ebpp, 2, function(a) meddiff(a,subtypes,6)))

meddiff_table <- data.frame(C1diff=c1meddiff, C2diff=c2meddiff, C3diff=c3meddiff, C4diff=c4meddiff, C5diff=c5meddiff, C6diff=c6meddiff)

write.csv(meddiff_table, file=paste0(data_dir, 'meddiff_table.csv'))

c1kt <- as.numeric(apply(ebpp, 2, function(a) kt_fun(a,subtypes,1)))
c2kt <- as.numeric(apply(ebpp, 2, function(a) kt_fun(a,subtypes,2)))
c3kt <- as.numeric(apply(ebpp, 2, function(a) kt_fun(a,subtypes,3)))
c4kt <- as.numeric(apply(ebpp, 2, function(a) kt_fun(a,subtypes,4)))
c5kt <- as.numeric(apply(ebpp, 2, function(a) kt_fun(a,subtypes,5)))
c6kt <- as.numeric(apply(ebpp, 2, function(a) kt_fun(a,subtypes,6)))

stat_table <- data.frame(C1stat=c1kt, C2stat=c2kt, C3stat=c3kt, C4stat=c4kt, C5stat=c5kt, C6stat=c6kt)

write.csv(stat_table, file=paste0(data_dir, 'stat_table.csv'))


stat_table <- read.csv(paste0(data_dir,'stat_table.csv'), header=T)

meds_table <- read.csv(paste0(data_dir,'meddiff_table.csv'))

c1_col <- ifelse( (abs(stat_table$C1stat) > 500) & (abs(meds_table$C1diff) > 2000), "red", "blue")
plot(y=stat_table$C1stat, x=meds_table$C1diff, col=c1_col)

c2_col <- ifelse( (abs(stat_table$C2stat) > 500) & (abs(meds_table$C2diff) > 2000), "red", "blue")
plot(y=stat_table$C2stat, x=meds_table$C2diff, col=c2_col)

