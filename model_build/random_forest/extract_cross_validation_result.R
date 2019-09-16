out_prefix <- "D:/Project/产品/TRANSCUP/model_build/random_forest/output/extract_CV/"

extract_confusion_matrix <- function(cm,ID,out_prefix_cm){
  overal_all_df <- as.data.frame(cm$overall)
  overal_all_df <- as.data.frame(t(overal_all_df))
  rownames(overal_all_df) <- f
  write.csv(overal_all_df,file = paste0(out_prefix_cm,"_overall_stat.csv"))
  
  # confusion matrix table
  cm_table_df <- as.matrix(cm$table)
  write.csv(cm_table_df,file = paste0(out_prefix_cm,"_table.csv"))
  
  # by class stat
  byClass_mat <- cm$byClass
  rownames(byClass_mat) <- gsub("Class: ","",rownames(byClass_mat))
  write.csv(byClass_mat,file = paste0(out_prefix_cm,"_byClass_stat.csv"))
}

all_selected_feature_names <- c("TOP200","TOP150","TOP100","log2FC_2.5","log2FC_3")

# extract default
for(f in all_selected_feature_names){
  model_CV <- readRDS(paste0("D:/Project/产品/TRANSCUP/model_build/random_forest/output/Model_RF_corss_validation_",f,".rds"))
  out_prefix_cm <- paste0(out_prefix,f,"_CV_confusionMatrix")
  cm <- model_CV[[f]][["RF"]][["cross_validation"]][["train.model"]][["finalModel"]][["confusion"]]
  extract_confusion_matrix(cm,f,out_prefix_cm)
  
  
}
