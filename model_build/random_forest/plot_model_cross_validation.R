
out_prefix <- "D:/Project/产品/TRANSCUP/model_build/random_forest/output/"
library(randomForest)
library(caret)
library(mlbench)
library(e1071)

set.seed(1234)


all_selected_feature_names <- c("TOP100","TOP150","TOP200","log2FC_2.5","log2FC_3")
for(f in all_selected_feature_names){
  model_data <-readRDS(paste0("D:/Project/产品/TRANSCUP/model_build/random_forest/output/Model_RF_corss_validation_",f,".rds"))
  pdf(paste0(out_prefix,f,"_accuracy_CV.pdf"),width=10,height = 7)
  print(plot(model_data[[f]][['RF']][['cross_validation']][['train.model']]))
  dev.off()
  
  
  bestTune <- model_data[[f]][['RF']][['cross_validation']][['train.model']][['bestTune']]
  write.csv(bestTune,file = paste0(out_prefix,f,"_best_tune.csv"))
  
  tune_results <- model_data[[f]][['RF']][['cross_validation']][['train.model']][['results']]
  write.csv(tune_results,file = paste0(out_prefix,f,"_tune_results.csv"))
  

  
  
}
