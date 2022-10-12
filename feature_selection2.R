library(readr)
library(stringr)
library(FNN)
library(ggplot2)

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


df1 <- read.csv(paste0(data_dir,'/features/C1_feature_sel_vals.csv'))
df2 <- read.csv(paste0(data_dir,'/features/C2_feature_sel_vals.csv'))
df3 <- read.csv(paste0(data_dir,'/features/C3_feature_sel_vals.csv'))
df4 <- read.csv(paste0(data_dir,'/features/C4_feature_sel_vals.csv'))
df5 <- read.csv(paste0(data_dir,'/features/C5_feature_sel_vals.csv'))
df6 <- read.csv(paste0(data_dir,'/features/C6_feature_sel_vals.csv'))


df1x <- df1[is.finite(rowSums(df1)),]
df1_up <- df1x[df1x$UP > 0.5,]
df1_dn <- df1x[df1x$UP < 0.5,]
idx <- order(df1_up$KL, decreasing = T)
jdx <- order(df1_dn$KL, decreasing = T)

df1_up_X <- df1_up[idx[1:20],'X']
df1_dn_X <- df1_dn[jdx[1:20],'X']

df1x$color <-  "grey"
df1x$color[df1x$X %in% df1_up_X] <- "red"
df1x$color[df1x$X %in% df1_dn_X] <- "blue"

qplot(y=df1x$KL, x=df1x$UP,
      xlab='Probability for greater than', ylab='KL divergence',
      col = df1x$color)


df2x <- df2[is.finite(rowSums(df2)),]
df2_up <- df2x[df2x$UP > 0.5,]
df2_dn <- df2x[df2x$UP < 0.5,]
idx <- order(df2_up$KL, decreasing = T)
jdx <- order(df2_dn$KL, decreasing = T)

df2_up_X <- df2_up[idx[1:20],'X']
df2_dn_X <- df2_dn[jdx[1:20],'X']

df2x$color <-  "grey"
df2x$color[df2x$X %in% df2_up_X] <- "red"
df2x$color[df2x$X %in% df2_dn_X] <- "blue"

qplot(y=df2x$KL, x=df2x$UP,
      xlab='Probability for greater than', ylab='KL divergence',
      col = df2x$color)



df3x <- df3[is.finite(rowSums(df3)),]
df3_up <- df3x[df3x$UP > 0.5,]
df3_dn <- df3x[df3x$UP < 0.5,]
idx <- order(df3_up$KL, decreasing = T)
jdx <- order(df3_dn$KL, decreasing = T)

df3_up_X <- df3_up[idx[1:20],'X']
df3_dn_X <- df3_dn[jdx[1:20],'X']

df3x$color <-  "grey"
df3x$color[df3x$X %in% df3_up_X] <- "red"
df3x$color[df3x$X %in% df3_dn_X] <- "blue"

qplot(y=df3x$KL, x=df3x$UP,
      xlab='Probability for greater than', ylab='KL divergence',
      col = df3x$color)


df4x <- df4[is.finite(rowSums(df4)),]
df4_up <- df4x[df4x$UP > 0.5,]
df4_dn <- df4x[df4x$UP < 0.5,]
idx <- order(df4_up$KL, decreasing = T)
jdx <- order(df4_dn$KL, decreasing = T)

df4_up_X <- df4_up[idx[1:20],'X']
df4_dn_X <- df4_dn[jdx[1:20],'X']

df4x$color <-  "grey"
df4x$color[df4x$X %in% df4_up_X] <- "red"
df4x$color[df4x$X %in% df4_dn_X] <- "blue"

qplot(y=df4x$KL, x=df4x$UP,
      xlab='Probability for greater than', ylab='KL divergence',
      col = df4x$color)



df5x <- df5[is.finite(rowSums(df5)),]
df5_up <- df5x[df5x$UP > 0.5,]
df5_dn <- df5x[df5x$UP < 0.5,]
idx <- order(df5_up$KL, decreasing = T)
jdx <- order(df5_dn$KL, decreasing = T)

df5_up_X <- df5_up[idx[1:20],'X']
df5_dn_X <- df5_dn[jdx[1:20],'X']

df5x$color <-  "grey"
df5x$color[df5x$X %in% df5_up_X] <- "red"
df5x$color[df5x$X %in% df5_dn_X] <- "blue"

qplot(y=df5x$KL, x=df5x$UP,
      xlab='Probability for greater than', ylab='KL divergence',
      col = df5x$color)


df6x <- df6[is.finite(rowSums(df6)),]
df6_up <- df6x[df6x$UP > 0.5,]
df6_dn <- df6x[df6x$UP < 0.5,]
idx <- order(df6_up$KL, decreasing = T)
jdx <- order(df6_dn$KL, decreasing = T)

df6_up_X <- df6_up[idx[1:20],'X']
df6_dn_X <- df6_dn[jdx[1:20],'X']

df6x$color <-  "grey"
df6x$color[df6x$X %in% df6_up_X] <- "red"
df6x$color[df6x$X %in% df6_dn_X] <- "blue"

qplot(y=df6x$KL, x=df6x$UP,
      xlab='Probability for greater than', ylab='KL divergence',
      col = df6x$color)


sf1up <- data.frame(Type='feature', Symbol=symbs[df1_up_X], X=df1_up_X, Cluster='C1_up')
sf1dn <- data.frame(Type='feature', Symbol=symbs[df1_dn_X], X=df1_dn_X, Cluster='C1_dn')

sf2up <- data.frame(Type='feature', Symbol=symbs[df2_up_X], X=df2_up_X, Cluster='C2_up')
sf2dn <- data.frame(Type='feature', Symbol=symbs[df2_dn_X], X=df2_dn_X, Cluster='C2_dn')

sf3up <- data.frame(Type='feature', Symbol=symbs[df3_up_X], X=df3_up_X, Cluster='C3_up')
sf3dn <- data.frame(Type='feature', Symbol=symbs[df3_dn_X], X=df3_dn_X, Cluster='C3_dn')

sf4up <- data.frame(Type='feature', Symbol=symbs[df4_up_X], X=df4_up_X, Cluster='C4_up')
sf4dn <- data.frame(Type='feature', Symbol=symbs[df4_dn_X], X=df4_dn_X, Cluster='C4_dn')

sf5up <- data.frame(Type='feature', Symbol=symbs[df5_up_X], X=df5_up_X, Cluster='C5_up')
sf5dn <- data.frame(Type='feature', Symbol=symbs[df5_dn_X], X=df5_dn_X, Cluster='C5_dn')

sf6up <- data.frame(Type='feature', Symbol=symbs[df6_up_X], X=df6_up_X, Cluster='C6_up')
sf6dn <- data.frame(Type='feature', Symbol=symbs[df6_dn_X], X=df6_dn_X, Cluster='C6_dn')

sf <- rbind(sf1up,sf1dn,sf2up,sf2dn,sf3up,sf3dn,sf4up,sf4dn,sf5up,sf5dn,sf6up,sf6dn)



load('F:/Work/cluster_prediction_data/geneSetSymbols.rda')

sigdf1 <- data.frame(Type='signature', Symbol=genesetsymbols$CHANG_CORE_SERUM_RESPONSE_UP, X=0, Cluster='chang_core_serum_response_up')
sigdf2 <- data.frame(Type='signature', Symbol=genesetsymbols$CSF1_response, X=0, Cluster='csf1_reponse')
sigdf3 <- data.frame(Type='signature', Symbol=genesetsymbols$LIexpression_score, X=0, Cluster='li_expression')
sigdf4 <- data.frame(Type='signature', Symbol=genesetsymbols$Module3_IFN_score, X=0, Cluster='module3_ifn_score')
sigdf5 <- data.frame(Type='signature', Symbol=genesetsymbols$TGFB_score_21050467, X=0, Cluster='tgfb_score_2105467')

sf <- rbind(sf,sigdf1,sigdf2,sigdf3,sigdf4,sigdf5)\

write.csv(sf,file='selected_features_11_10_22.csv',row.names = F)
