# train a model 

result_dir <- "./output/"

dir.create(result_dir)

library(randomForest)
library(caret)
set.seed(1234)
all33_data <- readRDS("../preprocess/output/04_make_R_dataset/all_data.rds")
#registerDoParallel(cores = 48)

#all_selected_feature_names <- c("TOP200","TOP150","TOP100","log2FC_1.5","log2FC_2","log2FC_2.5","log2FC_3")
all_selected_feature_names <- c("TOP150")

for(f in all_selected_feature_names){
  all33_data_result <- list()
  all_data <- all33_data$all_data
  all_data <- all_data[colnames(all33_data[[f]][["train"]]),]
  all_data <- t(all_data)
  all_data <- log2(all_data + 1)
  preprocessParams <- preProcess(all_data, method=c("center", "scale"))
  all33_data_result[[f]] <- list()
  all33_data_result[[f]][["RF"]] <- list()
  all33_data_result[[f]][["RF"]][["final"]] <- list()
  all33_data_result[[f]][["RF"]][["final"]][["preprocessParams"]] <- preprocessParams
  train_transformed <- predict(preprocessParams, all_data)
  train_transformed <- as.data.frame(train_transformed)
  train_transformed$tumor_type <- as.factor(all33_data$sample_info$tumor_type2)
  ntrees <- 1000
  ntrys <- 54
  train.model <- randomForest(tumor_type ~ .,data = train_transformed, ntree = ntrees, mtry = ntrys)
  all33_data_result[[f]][["RF"]][["final"]][["ntree"]] <- ntrees
  all33_data_result[[f]][["RF"]][["final"]][["ntry"]] <- ntrys
  all33_data_result[[f]][["RF"]][["final"]][["train.model"]] <- train.model
  saveRDS(all33_data_result, file = paste0(result_dir,"Model_RF_final_",f,".rds"))
}
