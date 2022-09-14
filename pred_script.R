#First to install the package.
library(readr)
library(stringr)
library(devtools)
require(robencla)


run_pred <- function(train_data, test_data, data_dir, out_dir, file_prefix, sigs) {

    ### MAKE SURE OUTPUT PATH IS OK TO WRITE ###
    if (!dir.exists(out_dir)){
      dir.create(out_dir)
    }

    #Now we're ready to do some training
    # our classifier object named Susan
    susan <- Robencla$new("susan")

    # xgboost parameters to pass to each sub-classifier in the ensembles
    params <- list(max_depth=12,
                   eta=0.2,
                   nrounds=24,
                   nthreads=4,
                   gamma=1,
                   lambda=1.25,
                   alpha=0.25,
                   verbose=0)

    # split the data, train and test
    susan$autotrain(data_file=paste0(data_dir,train_data),
                   label_name='Label',
                   sample_id = 'Barcode',
                   data_mode=c('sigpairs','pairs','quartiles'), # pairs,sigpairs,quartiles,tertiles,binarize,ranks,original #
                   signatures=sigs,
                   size=3,
                   params=params,
                   train_perc=0.8,
                   combine_function='median')

    susan$autotest(paste0(data_dir,test_data),
                  label_name='Label',
                  sample_id = 'Barcode')

    # Write out the predictions for each sample
    pred_results <- susan$results(include_label = T)
    pred_results$RunName <- file_prefix
    pred_results <- apply(pred_results,2,as.character)
    write.csv(pred_results, file=paste0(out_dir,file_prefix,'_pred_results.csv'), row.names = F, quote=F)

    # metrics on the test set predictions
    pred_metrics <- susan$classification_metrics(use_cv_results = F)  #
    pred_metrics$RunName <- file_prefix
    write_csv(pred_metrics,file=paste0(out_dir,file_prefix,'_pred_metrics.csv'))

    # and get the importance of features in each ensemble member
    imp_list <- susan$importance() # uses the last fold trained.
    imp_table <- data.frame()
    for (li in names(imp_list)){
      imp_list[[li]]$Label <- li
      imp_table <- rbind(imp_table, imp_list[[li]])
    }
    imp_table$RunName <- file_prefix
    write_csv(imp_table,file=paste0(out_dir,file_prefix,'_importance.csv'))


    ### Write out an info table
    info_df <- data.frame(Info=c('Run_Date', 'Training_Data', 'Test_Data', 'Data_Dir', 'Output_Dir', 'Prefix'),
                          Value=c(Sys.time(), train_data, test_data, data_dir, out_dir, prefix))

    return()
}


## EDIT THESE LOCATIONS##

###### SHOULD READ FROM CMD LINE ######

train_data <- 'EBPP_subset.csv'

test_data <- 'xena_rsem_fpkm_subset.csv'

### WHERE THE FORMATTED DATA DIR LIVES ###
data_dir <- 'D:/Work/cluster_prediction_data/formatted/'

### WHERE TO WRITE THE NEW TABLES ###
out_dir <- 'D:/Work/cluster_prediction_data/results4/'

### WHAT TO NAME THE FILES, WILL BE A COLUMN
prefix <- "Random_Genes_Sept13_run5"

### GENES TO INCLUDE AS A FILE ###
gene_table <- read.csv(paste0(data_dir,'random_genes_table.csv'), header=T, stringsAsFactors = F)
all_genes <- unique(gene_table$Gene)

# then create some gene signatures:
random_sigs <- list(
sig1=gene_table$Gene[gene_table$Type == 'Set1'],
sig2=gene_table$Gene[gene_table$Type == 'Set2'],
sig3=gene_table$Gene[gene_table$Type == 'Set3']
)

run_pred(train_data, test_data, data_dir, out_dir, prefix, random_sigs)