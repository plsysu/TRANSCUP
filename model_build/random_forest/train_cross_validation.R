# train a model 

result_dir <- "./output/"

dir.create(result_dir)

library(randomForest)
library(caret)
library(mlbench)
library(e1071)
library(doParallel)

set.seed(1234)
registerDoParallel(cores = 48)

all33_data <- readRDS("../preprocess/output/04_make_R_dataset/all_data.rds")


# Extend Caret
customRF <- list(type = "Classification", library = "randomForest", loop = NULL)
customRF$parameters <- data.frame(parameter = c("mtry", "ntree"), class = rep("numeric", 2), label = c("mtry", "ntree"))
customRF$grid <- function(x, y, len = NULL, search = "grid") {}
customRF$fit <- function(x, y, wts, param, lev, last, weights, classProbs, ...) {
  randomForest(x, y, mtry = param$mtry, ntree=param$ntree, ...)
}
customRF$predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata)
customRF$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata, type = "prob")
customRF$sort <- function(x) x[order(x[,1]),]
customRF$levels <- function(x) x$classes

# train model



all_selected_feature_names <- c("TOP100","TOP150","TOP200","log2FC_2","log2FC_2.5","log2FC_3")
#all_selected_feature_names <- c("TOP150")
for(f in all_selected_feature_names){
  all_data <- all33_data$all_data
  all_data <- all_data[colnames(all33_data[[f]][["train"]]),]
  all_data <- t(all_data)
  all_data <- log2(all_data + 1)
  
  all33_data_result <- list()
  all33_data_result[[f]] <- list()
  preprocessParams <- preProcess(all_data, method=c("center", "scale"))
  all33_data_result[[f]][["RF"]] <- list()
  all33_data_result[[f]][["RF"]][["cross_validation"]] <- list()
  all33_data_result[[f]][["RF"]][["cross_validation"]][["preprocessParams"]] <- preprocessParams
  train_transformed <- predict(preprocessParams, all_data)
  train_transformed <- as.data.frame(train_transformed)
  train_transformed$tumor_type <- as.factor(all33_data$sample_info$tumor_type2)
  ntree <- as.integer(round(((ncol(train_transformed)-1)/2)/1000)*1000)
  ntry <- as.integer(floor(sqrt(ncol(train_transformed)-1)))
  
  ntrees <- c(ntree-500,ntree,ntree+500,ntree+1000)
  ntrys <- c(ntry-4,ntry-2,ntry,ntry+2,ntry+4)
  
  tunegrid <- expand.grid(.mtry=ntrys, .ntree=ntrees)
  control <- trainControl(method="repeatedcv", number=10, repeats=3)
  
  all33_data_result[[f]][["RF"]][["cross_validation"]][["ntree"]] <- ntrees
  all33_data_result[[f]][["RF"]][["cross_validation"]][["ntry"]] <- ntrys
  all33_data_result[[f]][["RF"]][["cross_validation"]][["tunegrid"]] <- tunegrid
  all33_data_result[[f]][["RF"]][["cross_validation"]][["control"]] <- control
  
  metric <- "Accuracy"
  
  train.model <- train(tumor_type ~ ., data=train_transformed, method=customRF, metric=metric, tuneGrid=tunegrid, trControl=control)
  #train.model <- randomForest(tumor_type ~ .,data = train_transformed, ntree = ntree, mtry = ntry)
  
  all33_data_result[[f]][["RF"]][["cross_validation"]][["train.model"]] <- train.model
  saveRDS(all33_data_result, file = paste0(result_dir,"Model_RF_corss_validation_",f,".rds"))
}
