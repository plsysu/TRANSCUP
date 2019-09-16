rm(list = ls())
library(data.table)
library(dplyr)
set.seed(1234)
library(splitstackshape)
library(caret)

result_dir <- "./output/02_tumor_fpkm/"

dir.create(result_dir)

all33_fpkm_df <- fread('./output/01_raw_fpkm/all33_fpkm.csv', header = T, sep = ',',stringsAsFactors=F)
all33_fpkm_df <- as.data.frame(all33_fpkm_df)
rownames(all33_fpkm_df) <- all33_fpkm_df$V1
all33_fpkm_df$V1 <- NULL
all33_sample_info <- read.table("./output/01_raw_fpkm/all33_sample_info.tsv",header=T,sep="\t",stringsAsFactors=F)


# filter samples only keep: tumor
sample_names_keep <- NULL
sample_names_keep <- all33_sample_info$cases[all33_sample_info$tissue.definition=="Primary solid Tumor"]
sample_names_keep <- c(sample_names_keep,all33_sample_info$cases[all33_sample_info$tissue.definition=="Primary Blood Derived Cancer - Peripheral Blood"])
sample_names_keep <- c(sample_names_keep,all33_sample_info$cases[all33_sample_info$tissue.definition=="Recurrent Solid Tumor"])
sample_names_keep <- c(sample_names_keep,all33_sample_info$cases[all33_sample_info$tissue.definition=="Metastatic"])
sample_names_keep <- c(sample_names_keep,all33_sample_info$cases[all33_sample_info$tissue.definition=="Additional Metastatic"])
sample_names_keep <- c(sample_names_keep,all33_sample_info$cases[all33_sample_info$tissue.definition=="Additional - New Primary"])

all33_fpkm_df_tumor <- all33_fpkm_df %>% select(sample_names_keep)
all33_sample_info_tumor <- all33_sample_info[which(all33_sample_info$cases %in% sample_names_keep ),]


all33_sample_info_tumor$tumor_type2 <- all33_sample_info_tumor$tumor_type
all33_sample_info_tumor$tumor_type2[which(all33_sample_info_tumor$tumor_type2=="COAD")]="CRC"
all33_sample_info_tumor$tumor_type2[which(all33_sample_info_tumor$tumor_type2=="READ")]="CRC"



train.index <- createDataPartition(all33_sample_info_tumor$tumor_type, p = .8, list = FALSE)
all33_sample_info_tumor_train <- all33_sample_info_tumor[train.index,]
all33_sample_info_tumor_test <- all33_sample_info_tumor[-train.index,]
# order FPKM by sample info
all33_fpkm_df_tumor <- all33_fpkm_df_tumor[,all33_sample_info_tumor$cases]
all33_fpkm_df_tumor_train <- all33_fpkm_df_tumor[,all33_sample_info_tumor_train$cases]
all33_fpkm_df_tumor_test <- all33_fpkm_df_tumor[,all33_sample_info_tumor_test$cases]




write.table(all33_sample_info_tumor,file = paste0(result_dir,"all33_sample_info_tumor.tsv"),sep = "\t",row.names = F)
write.csv(all33_fpkm_df_tumor,file = paste0(result_dir,"all33_fpkm_df_tumor.csv"),row.names = T)


write.table(all33_sample_info_tumor_train,file = paste0(result_dir,"all33_sample_info_tumor_train.tsv"),sep = "\t",row.names = F)
write.csv(all33_fpkm_df_tumor_train,file = paste0(result_dir,"all33_fpkm_df_tumor_train.csv"),row.names = T)


write.table(all33_sample_info_tumor_test,file = paste0(result_dir,"all33_sample_info_tumor_test.tsv"),sep = "\t",row.names = F)
write.csv(all33_fpkm_df_tumor_test,file = paste0(result_dir,"all33_fpkm_df_tumor_test.csv"),row.names = T)

