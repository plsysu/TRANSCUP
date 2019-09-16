library(dplyr)

library(randomForest)
library(caret)
library(e1071)



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




f <- "TOP200"
model_cv_T200 <- readRDS(paste0("D:/Project/产品/TRANSCUP/model_build/random_forest/output/Model_RF_corss_validation_",f,".rds"))

true_class <- model_cv_T200$TOP200$RF$cross_validation$train.model$finalModel$y
predicted_class <- model_cv_T200$TOP200$RF$cross_validation$train.model$finalModel$predicted
confusion_matrix_T200 <- confusionMatrix(predicted_class,true_class)

out_prefix <- "D:/Project/产品/TRANSCUP/model_build/random_forest/output/extract_CV/"

out_prefix_cm <- paste0(out_prefix,f,"_CV_confusionMatrix")
extract_confusion_matrix(confusion_matrix_T200,f,out_prefix_cm)
