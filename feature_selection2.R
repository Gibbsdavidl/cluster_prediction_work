library(readr)
library(stringr)
library(FNN)
library(ggplot2)

### Selecting features ###

### WHERE THE DATA DIR LIVES ###
data_dir <- 'C:/Users/dgibbs/ISB_Work/cluster_prediction_data/'
#data_dir <- 'F:/Work/cluster_prediction_data/'

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


kl_fun <- function(vecx,subtypes,clusterval){
  try({
    x <- vecx[subtypes$ClusterModel1 == clusterval]
    y <- vecx[subtypes$ClusterModel1 != clusterval]
    res0 <- FNN::KL.divergence(x, y, k = 9, algorithm="kd_tree")
    xreturn(mean(res0,na.rm=TRUE))
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


make_c_table <- function(X, subtypes, clusterval, fileout) {
  c1me <- as.numeric(apply(ebpp_ranknorm, 2, function(a) meddiff(a,subtypes,clusterval)))
  c1kt <- as.numeric(apply(ebpp_ranknorm, 2, function(a) kt_fun(a,subtypes,clusterval)))
  c1KL <- as.numeric(apply(ebpp_ranknorm, 2, function(a) kl_fun(a,subtypes,clusterval)))
  c1up <- as.numeric(apply(ebpp_ranknorm, 2, function(a) samp_up(a,subtypes,clusterval)))
  c1dn <- as.numeric(apply(ebpp_ranknorm, 2, function(a) samp_dn(a,subtypes,clusterval)))
  df <- data.frame(MD=c1me, KT=c1kt, KL=c1KL, UP=c1up, DN=c1dn)
  write.csv(df, fileout)
  return(df)
}

fileout <- 'C:/Users/dgibbs/ISB_Work/cluster_prediction_data/C1_feature_sel_vals.csv'
df1 <- make_c_table(ebpp_ranknorm, subtypes, 1, fileout)

fileout <- 'C:/Users/dgibbs/ISB_Work/cluster_prediction_data/C2_feature_sel_vals.csv'
df2 <- make_c_table(ebpp_ranknorm, subtypes, 2, fileout)

fileout <- 'C:/Users/dgibbs/ISB_Work/cluster_prediction_data/C3_feature_sel_vals.csv'
df3 <- make_c_table(ebpp_ranknorm, subtypes, 3, fileout)

fileout <- 'C:/Users/dgibbs/ISB_Work/cluster_prediction_data/C4_feature_sel_vals.csv'
df4 <- make_c_table(ebpp_ranknorm, subtypes, 4, fileout)

fileout <- 'C:/Users/dgibbs/ISB_Work/cluster_prediction_data/C5_feature_sel_vals.csv'
df5 <- make_c_table(ebpp_ranknorm, subtypes, 5, fileout)

fileout <- 'C:/Users/dgibbs/ISB_Work/cluster_prediction_data/C6_feature_sel_vals.csv'
df6 <- make_c_table(ebpp_ranknorm, subtypes, 6, fileout)


plot(y=df1$KL, x=df1$UP, xlab='Probability for greater than', ylab='KL divergence')

plot(y=df2$KL, x=df2$UP, xlab='Probability for greater than', ylab='KL divergence')

plot(y=df3$KL, x=df3$UP, xlab='Probability for greater than', ylab='KL divergence')

plot(y=df4$KL, x=df4$UP, xlab='Probability for greater than', ylab='KL divergence')

plot(y=df5$KL, x=df5$UP, xlab='Probability for greater than', ylab='KL divergence')

plot(y=df6$KL, x=df6$UP, xlab='Probability for greater than', ylab='KL divergence')

#plot(y=c1KL, x=c1up, xlab='Probability for greater than', ylab='KL divergence')
#plot(y=c1KL, x=c1me, xlab='median difference between groups', ylab='KL divergence')
#plot(y=c1KL, x=c1kt, xlab='Kruskal stat', ylab='KL divergence')
#plot(x=c1up, y=c1kt, ylab='Kruskal stat', xlab='prob up')
#plot(x=c1up, y=c1dn, ylab='prob dn', xlab='prob up')
#df2[(df2$UP < 0.5) & (df2$DN < 0.4),]
#df2 <- df[abs(df$MD) > 0, ]
#plot(x=df2$UP, y=df2$DN, ylab='prob dn', xlab='prob up')


df2[(!is.na(df2$KL)) & (!is.infinite(df2$KL)) & (df2$KL > 0.5),]

df3[(!is.na(df3$KL)) & (!is.infinite(df3$KL)) & (df3$KL > 0.5),]