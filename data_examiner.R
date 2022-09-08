
library(readr)
library(stringr)

# The immune subtypes
### Use: ClusterModel1 as the subtype
subtypes <- readr::read_delim('data/five_signature_mclust_ensemble_results.tsv.gz', delim = '\t')
subtypes$AliquotBarcode <- str_replace_all(string = subtypes$AliquotBarcode, pattern = '\\.', replacement = '-')
subtypes$SampleBarcode2 <- str_sub(subtypes$SampleBarcode, start = 1, end = 15)
subtypes$SampleBarcode2 <- str_sub(subtypes$SampleBarcode, start = 1, end = 15)


# The Pan-Cancer EB++ (named EBpp) batch corrected data set 
ebpp <- readr::read_delim('data/EBpp/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv',delim = '\t')
ebpp_genes <- readr::read_csv('data/EBpp/ebpp_genes.txt')


# Rsubread data in both FPKM and TPM normalizations
Sys.setenv(VROOM_CONNECTION_SIZE=500072)
rsubread_fpkm <- readr::read_delim('data/Rsubread/GSM1536837_01_27_15_TCGA_20.Illumina.tumor_Rsubread_FPKM.tsv.gz',delim = '\t')
rsubread_tpm <- readr::read_delim('data/Rsubread/GSM1536837_01_27_15_TCGA_20.Illumina.tumor_Rsubread_TPM.tsv.gz', delim = '\t')


# Xena Kallisto TPM
xena_kallisto_tpm <- readr::read_delim('data/xena_kallisto_tpm/tcga_Kallisto_tpm.tsv.gz', delim='\t')
xena_gene_id_map <- readr::read_csv('data/xena_kallisto_tpm/gene_id_map.csv')


# xena rsem TPM
xena_rsem_tpm <- readr::read_delim('data/xena_rsem_tpm/xena_tcga_RSEM_gene_tpm.tsv.gz')
xena_gene_id_map <- readr::read_csv('data/xena_rsem_tpm/gene_id_map.csv')


# xena rsem FPKM
xena_rsem_fpkm <- readr::read_delim('data/xena_rsem_fpkm/tcga_RSEM_gene_fpkm.tsv.gz')
xena_gene_id_map <- readr::read_csv('data/xena_rsem_fpkm/gene_id_map.csv')


# The xena data have sample barcodes of length 15
length(intersect(subtypes$SampleBarcode2, colnames(xena_rsem_fpkm)))  # 8500


# The Rsubread data have columns as TCGA aliquot barcodes
length(intersect(colnames(rsubread_tpm), subtypes$AliquotBarcode))  # 6680

