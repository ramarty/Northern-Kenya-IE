# Northern Kenya IE

# I MADE A CHANGE. CHANGE CHANGE CHANGE.

# Setup ------------------------------------------------------------------------
if(Sys.info()["user"] == "WB521633") source("C:/Users/wb521633/Dropbox/World Bank/IEs/Northern Kenya IE/Code/northern_kenya_ie_master.R")
if(Sys.info()["user"] == "robmarty") source("~/Dropbox/World Bank/IEs/Northern Kenya IE/Code/northern_kenya_ie_master.R")

library(randomForest)

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
training_variables <- c(paste0("band_",1:11),"ndvi","ndwi","ndbi","evi","ui")
random_points$na_in_train_var <- random_points[,training_variables] %>% rowMeans %>% is.na
random_points <- random_points[!random_points$na_in_train_var,]

#### Standardize Variables
#for(var in training_variables) random_points[[var]] <- scale(random_points[[var]])

# Train and Test Sets ----------------------------------------------------------
#### Generate train variable
random_points$train <- sample(x=c(FALSE,TRUE), 
                              prob=c(.95,.05), 
                              size=nrow(random_points), 
                              replace=TRUE) 
if(F){
  random_points$train[random_points$built_up==T] <- sample(x=c(FALSE,TRUE), 
                                                           prob=c(.8,.2), 
                                                           size=length(random_points$train[random_points$built_up==T]), 
                                                           replace=TRUE) 
}

#### Separate into train and validation sets
TrainSet <- random_points[random_points$train==T,]
ValidSet <- random_points[random_points$train==F,]

TrainSet$built_up %>% table
ValidSet$built_up %>% table

# Implement RF Model -----------------------------------------------------------
rf_model1 <- randomForest(as.formula(paste0("built_up ~ ", paste(training_variables, collapse=" + "))), 
                          data = TrainSet, ntree = 20, mtry = 10, importance = FALSE)

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










