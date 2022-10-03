library(readr)
library(stringr)

### Selecting features ###

### WHERE THE DATA DIR LIVES ###
#data_dir <- 'C:/Users/dgibbs/ISB_Work/cluster_prediction_data/'
data_dir <- 'F:/Work/cluster_prediction_data/'

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
symbs <- gene_symbols[gene_dup]
symbs <- symbs[2:20502]

ebpp <- ebpp[gene_dup,]
rownames(ebpp) <- gene_symbols[gene_dup]

# subset data
# and get it into order
ebpp_sub <- ebpp[,subtypes$AliquotBarcode]
ebpp_sub <- ebpp_sub[2:20502,]

# we transpose it and add in the label
ebpp_ranked <- apply(ebpp_sub, 2, data.table::frank)

rank_norm <- function(x) {
  norx <- (x - min(x)) / (max(x) - min(x))
  return(norx)
}

ebpp_ranknorm <- t(apply(ebpp_ranked, 2, rank_norm))
colnames(ebpp_ranknorm) <- symbs

rm(ebpp, ebpp_ranked)
gc()

y <- ifelse(subtypes$ClusterModel1 == 1,1, 0)

idx <- which(subtypes$ClusterModel1 == 1)

table(y)


kl_fun <- function(vecx,subtypes,clusterval){
  try({
    x <- vecx[subtypes$ClusterModel1 == clusterval]
    y <- vecx[subtypes$ClusterModel1 != clusterval]
    res0 <- KL.divergence(x, y, k = 9, algorithm="kd_tree")
    return(mean(res0,na.rm=TRUE))
  })
  return(0.0)
}


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


samp_up <- function(vecx,subtypes,clusterval){
  try({
    x <- vecx[subtypes$ClusterModel1 == clusterval]
    y <- vecx[subtypes$ClusterModel1 != clusterval]
    xs <- sample(x,10000,TRUE)
    ys <- sample(y,10000,TRUE)
    return(sum(xs > ys)/10000)
  })
}

samp_dn <- function(vecx,subtypes,clusterval){
  try({
    x <- vecx[subtypes$ClusterModel1 == clusterval]
    y <- vecx[subtypes$ClusterModel1 != clusterval]
    xs <- sample(x,10000,TRUE)
    ys <- sample(y,10000,TRUE)
    return(sum(xs < ys)/10000)
  })
}


c1me <- as.numeric(apply(ebpp_ranknorm, 2, function(a) meddiff(a,subtypes,1)))
c1kt <- as.numeric(apply(ebpp_ranknorm, 2, function(a) kt_fun(a,subtypes,1)))
c1KL <- as.numeric(apply(ebpp_ranknorm, 2, function(a) kl_fun(a,subtypes,1)))
c1up <- as.numeric(apply(ebpp_ranknorm, 2, function(a) samp_up(a,subtypes,1)))
c1dn <- as.numeric(apply(ebpp_ranknorm, 2, function(a) samp_dn(a,subtypes,1)))


