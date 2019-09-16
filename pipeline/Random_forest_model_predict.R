#!/usr/bin/env Rscript
###############################################################################################
##   Author: Peng Li
##   Email: peng-li@outlook.com
##   Created: 2019/05/27
##   Version: 0.0.1
##   Description: This script is used to predict tumor type by RandomForest model.
##   Requirments: dplyr, optparse, randomForest, caret, e1071
##   usage example:
##          
##   Output format:
###############################################################################################
library(optparse)
library(dplyr)

library(randomForest)
library(caret)
library(e1071)

set.seed(25)

option_list = list(
  make_option(c("-s", "--sample"), type="character", default=NULL,
              help="The sample's name ( required ).", metavar="character"),
  make_option(c("-d", "--data"), type="character", default=NULL,
              help="The sample's FPKM  table ( abs path, required ).", metavar="character"),
  make_option(c("-m", "--model"), type="character", default=NULL,
              help="The random forest model ( abs path, required ).", metavar="character"),
  make_option(c("-o", "--out_prefix"), type="character", default=NULL,
              help="The output prefix.( abs path, required)", metavar="character")
  
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
# check input sample name
if (is.null(opt$sample)){
  print_help(opt_parser)
  stop("The sample's name  must be supplied.", call.=FALSE)
}

# check input FPKM table 
if (is.null(opt$data)){
  print_help(opt_parser)
  stop("The sample's FPKM  table must be supplied.", call.=FALSE)
}

# check input model
if (is.null(opt$model)){
  print_help(opt_parser)
  stop("The random forest model must be supplied.", call.=FALSE)
}

# check output prefix
if (is.null(opt$out_prefix)){
  print_help(opt_parser)
  stop("The output prefix( abs path, required)", call.=FALSE)
}


# get args
sample_name <- opt$sample
FPKM_file <- opt$data
model_file <- opt$model
output_prefix <- opt$out_prefix

# output file name
output_predict <- paste0(output_prefix,"_RF_pred_tumor_type.csv")
output_predict_prob <- paste0(output_prefix,"_RF_pred_probs.csv")

model <- readRDS(model_file)
selected_features <- names(model$TOP200$RF$cross_validation$preprocessParams$mean)
preprocessParams <- model$TOP200$RF$cross_validation$preprocessParams
train.model <- model$TOP200$RF$cross_validation$train.model$finalModel

# read FPKM table
fpkm_df <- read.table(FPKM_file,header = T,sep = "\t",stringsAsFactors = F)
test_df <- data.frame(fpkm_df$fpkm)
rownames(test_df) <- fpkm_df$gene_name
colnames(test_df) <- sample_name
test_df <- as.data.frame(t(test_df),stringsAsFactors = F)

test_df_select_feature <- test_df[,selected_features]
test_df_select_feature_log <- log2(test_df_select_feature+1)
test_df_select_feature_log_transformed <- predict(preprocessParams,test_df_select_feature_log)
RF.pred <- predict(train.model, newdata = test_df_select_feature_log_transformed)
RF.pred_df <- as.data.frame(RF.pred)
RF.pred_df <- data.frame(sample=rownames(RF.pred_df),predict_class=RF.pred_df$RF.pred,stringsAsFactors = F)
write.csv(RF.pred_df,file = output_predict,row.names = T)

RF.pred.prob <- predict(train.model, newdata = test_df_select_feature_log_transformed, type = "prob")
write.csv(RF.pred.prob,file = output_predict_prob,row.names = T)
