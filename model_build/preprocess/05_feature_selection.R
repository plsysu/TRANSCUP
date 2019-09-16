#rm(list = ls())
library(data.table)
library(dplyr)

result_dir <- "./output/03_featureSelection/"

dir.create(result_dir)
# for train dataset:
#     log2 transform gene exp
#     filter low expression genes
#     filter low variance genes
#     feature selectation

# reorder exp matrix column match sample info's order
all33_fpkm_df_tumor_train <- fread('./output/02_tumor_fpkm/all33_fpkm_df_tumor.csv', header = T, sep = ',',stringsAsFactors=F)
all33_fpkm_df_tumor_train <- as.data.frame(all33_fpkm_df_tumor_train)
rownames(all33_fpkm_df_tumor_train) <- all33_fpkm_df_tumor_train$V1
all33_fpkm_df_tumor_train$V1 <- NULL

all33_sample_info_tumor_train <- read.table("./output/02_tumor_fpkm/all33_sample_info_tumor.tsv",header=T,sep="\t",stringsAsFactors=F)
all33_fpkm_df_tumor_train.log <- log2(all33_fpkm_df_tumor_train+1)

# filter genes
gene.max <- apply(all33_fpkm_df_tumor_train.log,1,max)
gene.median <- apply(all33_fpkm_df_tumor_train.log,1,median)
gene.mean <- apply(all33_fpkm_df_tumor_train.log,1,mean)
gene.sd  <- apply(all33_fpkm_df_tumor_train.log,1,sd)

# select  1st Qu.
all33_fpkm_df_tumor_train.log.filter <- all33_fpkm_df_tumor_train.log[gene.max > 2.5 & gene.sd > 0.1,]
# feature selection

find.tumor.specific.features <- function(exp,sample.type, type="ACC")
{
  idx1 <- which(sample.type == type)
  idx2 <- which(sample.type != type)	
  
  t.test.fun <- function(x, idx1,idx2)
  {
    return(t.test(x[idx1],x[idx2])$p.value)
  }
  
  
  logFC.fun <- function(x, idx1,idx2)
  {
    return(mean(x[idx1]) -  mean(x[idx2]))
  }
  logFC <- apply(exp,1,logFC.fun,idx1,idx2)
  pvalue <- apply(exp,1,t.test.fun,idx1,idx2)
  res <- data.frame(pvalue,logFC)
  rownames(res) <- rownames(exp)
  #res.sel <- res[res$pvalue < 0.001,]
  res.sel <- res
  res.ord <- res.sel[order(res.sel$logFC,decreasing=T),]
  res.ord$tumpr_type <- type
  return(res.ord)
}



types <- unique(all33_sample_info_tumor_train$tumor_type)
class.labels <- all33_sample_info_tumor_train$tumor_type
feature_res <- NULL
for(i in 1:length(types))
{
  
  res <- find.tumor.specific.features(all33_fpkm_df_tumor_train.log.filter,class.labels, type=types[i])
  # select top 100 overexpression genes
  #res <- res[1:num,]
  res$gene_name <- rownames(res)
  feature_res <- rbind(feature_res,res)
  #sel.features <- union(sel.features,res$gene_name)
}
write.csv(feature_res,file = paste0(result_dir,"feature_selection_res_tab_","all.RAW",".csv"),row.names = F) 

# select features TOP100,150,200
for(n in c(100,150,200)){
  selected_all <- NULL
  for(i in 1:length(types)){
    feature_res_sub <- feature_res[which(feature_res$tumpr_type==types[i]),]
    feature_res_sub <- feature_res_sub[which(feature_res_sub$pvalue < 0.001),]
    res <- feature_res_sub[1:n,]
    selected_all <- rbind(selected_all,res)
    
    
  }
  unique_features <- unique(selected_all$gene_name)
  write.csv(selected_all,file = paste0(result_dir,"feature_selection_res_tab_TOP",n,".csv"),row.names = F) 
  write.csv(unique_features,file = paste0(result_dir,"feature_selection_unique_TOP",n,".csv"),row.names = F) 
}

# select features log2FC > 2,2.5,3
for(n in c(2,2.5,3)){
  selected_all <- NULL
  for(i in 1:length(types)){
    feature_res_sub <- feature_res[which(feature_res$tumpr_type==types[i]),]
    feature_res_sub <- feature_res_sub[which(feature_res_sub$pvalue < 0.001),]
    res <- feature_res_sub[which(feature_res_sub$logFC > n),]
    selected_all <- rbind(selected_all,res)
    
    
  }
  unique_features <- unique(selected_all$gene_name)
  write.csv(selected_all,file = paste0(result_dir,"feature_selection_res_tab_log2FC_",n,".csv"),row.names = F) 
  write.csv(unique_features,file = paste0(result_dir,"feature_selection_unique_log2FC_",n,".csv"),row.names = F)
}

