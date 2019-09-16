#rm(list = ls())
library(data.table)
library(dplyr)


result_dir <- "./output/04_make_R_dataset/"

dir.create(result_dir)

# load FPKM data all
all33_fpkm_df_tumor <- fread('./output/02_tumor_fpkm/all33_fpkm_df_tumor.csv', header = T, sep = ',',stringsAsFactors=F)
all33_fpkm_df_tumor<- as.data.frame(all33_fpkm_df_tumor)
rownames(all33_fpkm_df_tumor) <- all33_fpkm_df_tumor$V1
all33_fpkm_df_tumor$V1 <- NULL


# load FPKM data train
all33_fpkm_df_tumor_train <- fread('./output/02_tumor_fpkm/all33_fpkm_df_tumor_train.csv', header = T, sep = ',',stringsAsFactors=F)
all33_fpkm_df_tumor_train <- as.data.frame(all33_fpkm_df_tumor_train)
rownames(all33_fpkm_df_tumor_train) <- all33_fpkm_df_tumor_train$V1
all33_fpkm_df_tumor_train$V1 <- NULL

# load FPKM data test
all33_fpkm_df_tumor_test <- fread('./output/02_tumor_fpkm/all33_fpkm_df_tumor_test.csv', header = T, sep = ',',stringsAsFactors=F)
all33_fpkm_df_tumor_test <- as.data.frame(all33_fpkm_df_tumor_test)
rownames(all33_fpkm_df_tumor_test) <- all33_fpkm_df_tumor_test$V1
all33_fpkm_df_tumor_test$V1 <- NULL

all33_sample_info_tumor <- read.table("./output/02_tumor_fpkm/all33_sample_info_tumor.tsv",header=T,sep="\t",stringsAsFactors=F)
all33_sample_info_tumor_train <- read.table("./output/02_tumor_fpkm/all33_sample_info_tumor_train.tsv",header=T,sep="\t",stringsAsFactors=F)
all33_sample_info_tumor_test <- read.table("./output/02_tumor_fpkm/all33_sample_info_tumor_test.tsv",header=T,sep="\t",stringsAsFactors=F)

all33_data <- list()

all33_data$sample_info <- all33_sample_info_tumor
all33_data$sample_info_train <- all33_sample_info_tumor_train
all33_data$sample_info_test <- all33_sample_info_tumor_test



all33_data$all_data <- all33_fpkm_df_tumor
all33_data$all_data_train <- all33_fpkm_df_tumor_train
all33_data$all_data_test <- all33_fpkm_df_tumor_test



all_selected_feature_names <- c("TOP100","TOP150","TOP200","log2FC_2","log2FC_2.5","log2FC_3")

for(f in all_selected_feature_names){
  gene_table <- read.csv(file = paste0("./output/03_featureSelection/","feature_selection_unique_",f,".csv"),header = T,stringsAsFactors = F)
  selected_features <- gene_table$x
  # random sample selected features
  selected_features <- selected_features[sample(length(selected_features))]
  selected_df_train <- all33_data$all_data_train[selected_features,]
  selected_df_test <- all33_data$all_data_test[selected_features,]
  all33_data[[f]] <- list()
  all33_data[[f]][["train"]] <- as.data.frame(t(selected_df_train))
  all33_data[[f]][["test"]] <- as.data.frame(t(selected_df_test))
  
}

saveRDS(all33_data, file = paste0(result_dir,"all_data.rds"))
