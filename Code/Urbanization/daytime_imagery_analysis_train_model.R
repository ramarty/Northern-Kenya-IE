# Northern Kenya IE

# Setup ------------------------------------------------------------------------
if(Sys.info()["user"] == "WB521633") source("C:/Users/wb521633/Dropbox/World Bank/IEs/Northern Kenya IE/Code/northern_kenya_ie_master.R")
if(Sys.info()["user"] == "robmarty") source("~/Dropbox/World Bank/IEs/Northern Kenya IE/Code/northern_kenya_ie_master.R")
if(Sys.info()["user"] == "r521633") source("/home/wb521633/IEs/Northern-Kenya-IE/Code/northern_kenya_ie_master.R")

library(randomForest)
library(dplyr)

set.seed(42)

# Load and Prep Data -----------------------------------------------------------
#### Load Data
load(file.path(intermediate_data_file_path, "Africover Random Points with Landsat", "random_points_landsat.Rda"))

#### Create Variables 
random_points$built_up <- as.factor(random_points$africover_kenya == 8) # %>% factor(labels = c("non-built-up", "built-up"))
random_points$ndvi <- (random_points$band_5 - random_points$band_4) / (random_points$band_5 + random_points$band_4)
random_points$ndwi <- (random_points$band_3 - random_points$band_5) / (random_points$band_3 + random_points$band_5)
random_points$ndbi <- (random_points$band_6 - random_points$band_5) / (random_points$band_6 + random_points$band_5)
random_points$evi <- 2.5 * ((random_points$band_5/random_points$band_4)/(random_points$band_5 +6*random_points$band_4 - 7.5*random_points$band_2 +1))
random_points$ui <- (random_points$band_7 - random_points$band_5) / (random_points$band_7 + random_points$band_5)

#### Nighttime Lights
for(i in 1:12) random_points[[paste0("viirs_rad_",i)]][random_points[[paste0("viirs_cf_cvg_",i)]] %in% 0] <- NA # if no cloud free days, set radiance as 0.
random_points$viirs_rad_mean <- random_points[,names(random_points)[grepl(("viirs_rad_"), names(random_points))]] %>% rowMeans(na.rm=T) # apply(1,median)

#### Limit variables to only those used in model
training_variables <- c(paste0("band_",1:11),"ndvi","ndwi","ndbi","evi","ui","viirs_rad_mean")
random_points$na_in_train_var <- random_points[,training_variables] %>% rowMeans %>% is.na
random_points <- random_points[!random_points$na_in_train_var,]

#### Standardize Variables
#for(var in training_variables) random_points[[var]] <- scale(random_points[[var]])

# Train and Test Sets ----------------------------------------------------------
percent_train <- 0.01
urban_percent_total_truth <- 0.0013308

#### Create validate set
# Make validation set look like composition of truth data (in terms of distribution of classes)
random_points$validate_nonbuiltup[random_points$built_up == F] <- sample(x=c(FALSE,TRUE), 
                                                                  prob=c(percent_train,1-percent_train), 
                                                                  size=sum(random_points$built_up == F), 
                                                                  replace=TRUE)
number_urban_cells_validate <- ((sum(random_points$validate_nonbuiltup, na.rm=T) * urban_percent_total_truth) / (1-urban_percent_total_truth)) %>% ceiling
random_points$validate_builtup[random_points$built_up == T] <- sample(x=c(FALSE,TRUE), 
                                                                         prob=c(1-number_urban_cells_validate/sum(random_points$built_up == T),
                                                                                number_urban_cells_validate/sum(random_points$built_up == T)), 
                                                                         size=sum(random_points$built_up == T), 
                                                                         replace=TRUE)
random_points$validate <- (random_points$validate_nonbuiltup %in% T) | (random_points$validate_builtup %in% T)

random_points$train <- random_points$validate %in% F

#### Create train set that looks like composition of truth data
number_urban_cells_train <- ((sum(random_points$train[random_points$built_up == F]) * urban_percent_total_truth) / (1-urban_percent_total_truth)) %>% ceiling
random_points$train_builtup[(random_points$built_up == T) & (random_points$train == T)] <- sample(x=c(FALSE,TRUE),
                                                                                                  prob=c(1-sum(number_urban_cells_train)/sum((random_points$built_up == T) & (random_points$train == T)),
                                                                                                         sum(number_urban_cells_train)/sum((random_points$built_up == T) & (random_points$train == T))),
                                                                                                  size=sum((random_points$built_up == T) & (random_points$train == T)),
                                                                                                  replace=T)

random_points$train_classprop <- FALSE
random_points$train_classprop[random_points$train == T & random_points$built_up == F] <- TRUE
random_points$train_classprop[random_points$train_builtup == T] <- TRUE

# Confirm that train and validation are separate
table((random_points$validate==T & random_points$train == F) | (random_points$validate==F & random_points$train == T))
sum(!is.na(random_points$train_builtup) & random_points$train == T)
sum(!is.na(random_points$train_builtup))

#### Separate into train and validation sets
TrainSet <- random_points[random_points$train %in% T,]
TrainSet_classprop <- random_points[random_points$train_classprop  %in% T,]
ValidSet <- random_points[random_points$train %in% F,]

TrainSet <- TrainSet[order(runif(nrow(TrainSet))),] 
TrainSet_50000 <- rbind(
                    TrainSet[TrainSet$built_up == T,][1:25000,],
                    TrainSet[TrainSet$built_up == F,][1:25000,]
                  )
TrainSet_100000 <- rbind(
                    TrainSet[TrainSet$built_up == T,][1:50000,],
                    TrainSet[TrainSet$built_up == F,][1:50000,]
                    )

# Validation Accuracy Dataframe ------------------------------------------------
calc_validation_accuracy <- function(rf_model, training_vars){
  predValid <- predict(rf_model, ValidSet, type = "class")
  
  total_builtup <- sum(ValidSet$built_up == TRUE) # Number of built up pixels (truth)
  predict_builtup_accurate <- sum(predValid == TRUE & ValidSet$built_up == TRUE) # Accurately predicted urban
  predict_builtup_inaccurate <- sum(predValid == TRUE & ValidSet$built_up == FALSE) # False positive
  
  predict_builtup_accurate_perBU <- sum(predValid == TRUE & ValidSet$built_up == TRUE) / sum(ValidSet$built_up == TRUE)
  predict_builtup_inaccurate_perBU <- sum(predValid == TRUE & ValidSet$built_up == FALSE) / sum(ValidSet$built_up == TRUE)
  
  out <- c(
      total_builtup %>% as.character,
      predict_builtup_accurate %>% as.character,
      predict_builtup_inaccurate %>% as.character,
      predict_builtup_accurate_perBU %>% round(4) %>% as.character,
      predict_builtup_inaccurate_perBU %>% round(4) %>% as.character,
      training_vars %>% as.character,
      rf_model$ntree %>% as.character,
      rf_model$mtry %>% as.character
      )
  return(out)
}

results_df <- c("Total Built Up",
                "Predict Built Up: Accurate",
                "Predict Built Up: Inaccurate",
                "Predict Build Up Accurate (% Built Up)",
                "Predict Build Up Inaccurate (% Built Up)",
                "Training Vars",
                "ntree",
                "mtry") %>% as.data.frame %>% rename(categroy=".")



# Implement RF Model -----------------------------------------------------------
ntree <- rep(c(10,50,100),2)
mtry <- c(rep(5,3), rep(10,3))

for(i in 1:6){
  rf_model <- randomForest(as.formula(paste0("built_up ~ ", paste(training_variables, collapse=" + "))), 
                           data = TrainSet_classprop, 
                           ntree = ntree[i], mtry = mtry[i], importance = FALSE)
  results_df[[paste0("classprop_",i)]] <- calc_validation_accuracy(rf_model, "all")
  
  rf_model <- randomForest(as.formula(paste0("built_up ~ ", paste(training_variables, collapse=" + "))), 
                           data = TrainSet_50000, 
                           ntree = ntree[i], mtry = mtry[i], importance = FALSE)
  results_df[[paste0("train50k_",i)]] <- calc_validation_accuracy(rf_model, "all")
  
  print(i)
}

rf_model <- randomForest(as.formula(paste0("built_up ~ ", paste(training_variables, collapse=" + "))), 
                         data = TrainSet_classprop, 
                         ntree = 20, mtry = 10, importance = FALSE)
varImpPlot(rf_model)
results_df$model1 <- calc_validation_accuracy(rf_model)$value

rf_model <- randomForest(as.formula(paste0("built_up ~ ", paste(training_variables, collapse=" + "))), 
                         data = TrainSet_50000, 
                         ntree = 20, mtry = 10, importance = FALSE)
varImpPlot(rf_model)
results_df$model2 <- calc_validation_accuracy(rf_model)$value

results_df

# Analysis ---------------------------------------------------------------------
#### Train Set
predTrain <- predict(rf_model1, TrainSet, type = "class")
table(predTrain, TrainSet$built_up) 

#### Validation Set
predValid <- predict(rf_model1, ValidSet, type = "class")
table(predValid, ValidSet$built_up)

sum(ValidSet$built_up == TRUE) # Number of built up pixels (truth)
sum(predValid == TRUE & ValidSet$built_up == TRUE) # Accurately predicted urban
sum(predValid == TRUE & ValidSet$built_up == FALSE) # False positive

sum(predValid == TRUE & ValidSet$built_up == TRUE) / sum(ValidSet$built_up == TRUE)
sum(predValid == TRUE & ValidSet$built_up == FALSE) / sum(ValidSet$built_up == TRUE)

#### False Positives by Class
false_positives_by_class <- function(class){
  sum(predValid[ValidSet$africover_kenya == class] == TRUE & ValidSet$built_up[ValidSet$africover_kenya == class] == FALSE) %>% return()
}

false_positives_class_df <- lapply(0:10, false_positives_by_class) %>% 
  unlist %>% 
  as.data.frame %>% 
  rename(N = ".")
false_positives_class_df$percent_total <- false_positives_class_df$N / sum(predValid == TRUE & ValidSet$built_up == FALSE)
false_positives_class_df$class <- 0:10

class_proportions <- read.csv(file.path(intermediate_data_file_path, "Africover Random Points with Landsat", "africover_class_proportions.csv"))
class_proportions$relative_percent_urban <- class_proportions$Proportion / class_proportions$Proportion[class_proportions$class == 8]
class_proportions <- subset(class_proportions, select=c(class,relative_percent_urban))

false_positives_class_df <- merge(false_positives_class_df, class_proportions, by="class")
false_positives_class_df$N_weighted <- false_positives_class_df$N * false_positives_class_df$relative_percent_urban

#### Variable Importance
importance(rf_model1)        
varImpPlot(rf_model1)  










