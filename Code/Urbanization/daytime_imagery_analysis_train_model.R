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
# Indices
random_points$built_up <- as.factor(random_points$africover_kenya == 8) # %>% factor(labels = c("non-built-up", "built-up"))
random_points$ndvi <- (random_points$band_5 - random_points$band_4) / (random_points$band_5 + random_points$band_4)
random_points$ndwi <- (random_points$band_3 - random_points$band_5) / (random_points$band_3 + random_points$band_5)
random_points$ndbi <- (random_points$band_6 - random_points$band_5) / (random_points$band_6 + random_points$band_5)
random_points$evi <- 2.5 * ((random_points$band_5/random_points$band_4)/(random_points$band_5 +6*random_points$band_4 - 7.5*random_points$band_2 +1))
random_points$ui <- (random_points$band_7 - random_points$band_5) / (random_points$band_7 + random_points$band_5)

# NDSV Variables
create_ndsv_band_combinations <- function(i){
  band1 <- rep(i,(11-i))
  band2 <- rep((i+1):11) 
  
  df_out <- cbind(band1, band2) %>% as.data.frame
  
  return(df_out)
}

ndsv_combs_df <- lapply(1:10, create_ndsv_band_combinations) %>% bind_rows

for(i in 1:nrow(ndsv_combs_df)){
  random_points[[paste0("ndsv_",i)]] <- (random_points[[paste0("band_",ndsv_combs_df$band1[i])]] - random_points[[paste0("band_",ndsv_combs_df$band2[i])]]) / (random_points[[paste0("band_",ndsv_combs_df$band1[i])]] + random_points[[paste0("band_",ndsv_combs_df$band2[i])]])
}


# Nighttime Lights
for(i in 1:12) random_points[[paste0("viirs_rad_",i)]][random_points[[paste0("viirs_cf_cvg_",i)]] %in% 0] <- NA # if no cloud free days, set radiance as 0.
random_points$viirs_rad_mean <- random_points[,names(random_points)[grepl(("viirs_rad_"), names(random_points))]] %>% rowMeans(na.rm=T) # apply(1,median)

#### Limit variables to only those used in model
training_variables <- c(paste0("band_",1:11),"ndvi","ndwi","ndbi","evi","ui")
random_points$na_in_train_var <- random_points[,c(paste0("band_",1:11),"ndvi","ndwi","ndbi","evi","ui","viirs_rad_mean")] %>% rowMeans %>% is.na
random_points <- random_points[!random_points$na_in_train_var,]

#### Standardize Variables
#for(var in training_variables) random_points[[var]] <- scale(random_points[[var]])

# Train and Test Sets ----------------------------------------------------------
create_train_valid_sets <- function(train_set_bu_size, train_set_nbu_size, grid_id){
  
  # Select train_set_bu_size or 90% of train (whichever smallest)
  random_points_i <- random_points[random_points$grid_id %in% grid_id,]
  random_points_i$train <- F
  random_points_i$train[sample(which(random_points_i$built_up == F), train_set_bu_size, replace=F)] <- T
  random_points_i$train[sample(which(random_points_i$built_up == T), 
                               min(train_set_nbu_size, sum(random_points_i$built_up == T)*.9), 
                               replace=F)] <- T
  
  # Make sure non-built-up areas are definitely not urban
  random_points_i$train[(random_points_i$built_up == F) & (random_points_i$viirs_rad_mean > .1)] <- F
  
  trainset <- random_points_i[random_points_i$train == T,]
  validset <- random_points_i[random_points_i$train == F,]
  
  return(list(trainset=trainset,
              validset=validset))
}

grid1_50000_traintest <- create_train_valid_sets(20000,25000, 1)
grid2_50000_traintest <- create_train_valid_sets(20000,25000, 2)
grid3_50000_traintest <- create_train_valid_sets(20000,25000, 3)
grid4_50000_traintest <- create_train_valid_sets(20000,25000, 4)
grid5_50000_traintest <- create_train_valid_sets(20000,25000, 5)
grid6_50000_traintest <- create_train_valid_sets(20000,25000, 6)

grid1_50000_traintest <- grid5_50000_traintest

# Implement RF Model -----------------------------------------------------------
training_variables <- c(paste0("band_",1:11),"ndvi","ndwi","ndbi","evi","ui")
training_variables <- c(paste0("band_",1:11),paste0("ndsv_",1:55))

rf_model <- randomForest(as.formula(paste0("built_up ~ ", paste(training_variables, collapse=" + "))), 
                         data = grid1_50000_traintest$trainset, 
                         ntree = 600, importance = FALSE)

rf_model
varImpPlot(rf_model)  

grid1_50000_traintest$validset$built_up %>% table
predValid <- predict(rf_model, grid1_50000_traintest$validset, type="prob")

thrsh <- .985
table(predValid[,2]>thrsh, grid1_50000_traintest$validset$built_up)


plot(grid1_50000_traintest$validset$coord_x[grid1_50000_traintest$validset$built_up == F & predValid[,2]>thrsh],
     grid1_50000_traintest$validset$coord_y[grid1_50000_traintest$validset$built_up == F & predValid[,2]>thrsh])

points(grid1_50000_traintest$validset$coord_x[grid1_50000_traintest$validset$built_up == T],
       grid1_50000_traintest$validset$coord_y[grid1_50000_traintest$validset$built_up == T],
       col="red")



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

table(predValid)


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

random_points_BU <- random_points[random_points$built_up == T,]
random_points_NBU <- random_points[random_points$built_up == F,]


plot(random_points$coord_x[random_points$built_up == T], random_points$coord_y[random_points$built_up == T], pch=16,cex=.2)
plot(random_points$coord_x[random_points$train == F][predValid == T], random_points$coord_y[random_points$train == F][predValid == T], pch=16,cex=.2,col="red")

random_points$coord_x[random_points$train == F][predValid == T]


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






# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# ---------------------------- Purgatory ---------------------------------------
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

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

#TrainSet <- TrainSet[(TrainSet$viirs_rad_mean > .75 & TrainSet$built_up == F),]