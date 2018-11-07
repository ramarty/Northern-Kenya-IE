# Northern Kenya IE

# Nighttime light analysis
# 1. Use threshold approach to examine how NTL have changed over time.
#  1.1 Increase in number of cells above threshold
#  1.2 Increase in NTL for cells above threshold
# 2. Do NTL thresholds exclude urban areas?

# Setup ------------------------------------------------------------------------
if(Sys.info()["user"] == "WB521633") source("C:/Users/wb521633/Dropbox/World Bank/IEs/Northern Kenya IE/Code/northern_kenya_ie_master.R")

library(raster)
library(rgdal)
library(velox)
library(doBy)
library(ggplot2)
library(gridExtra)
library(ggpubr)

buffer_size <- 20 #km
create_viirs_shapefile <- FALSE

# Create/Load VIIRS-Level Polygon with Africover Data --------------------------
if(create_viirs_shapefile){
  convert_africover_builtup <- function(r){
    r[] <- as.numeric(r[] %in% 8)
    return(r)
  }
  
  convert_globcover_builtup <- function(r){
    r[] <- as.numeric(r[] %in% 190)
    return(r)
  }
  
  #### Study Area
  setwd(file.path(intermediate_data_file_path, "Treatment Roads Buffer"))
  roads_study_area <- readOGR(dsn=".", layer=paste0("kenya_treat_roads_",buffer_size,"km_buff"))
  roads_study_area$study_area <- 1
  
  #### Kenya ADM
  #setwd(file.path(raw_data_file_path,"GADM"))
  #ken_adm0 <- getData('GADM', country='KEN', level=0)
  
  #### Africover
  africover <- raster(file.path(intermediate_data_file_path, "Africover", "africover_kenya.tif"))
  
  #### Globcover
  globcover_2015 <- raster(file.path(raw_data_file_path, "esa_globcover", "scratch", "ESACCI-LC-L4-LCCS-Map-300m-P1Y-1992_2015-v2.0.7.tif"), band=24) %>% crop(roads_study_area)
  
  #### Create VIIRS-Level Polygon with VIIRS Data
  
  # Load VIIRS Data
  ken_viirs_2012_median <- raster(file.path(raw_data_file_path,"NTL Rasters Annual","ken_viirs_2012_median.tif")) %>% crop(roads_study_area) %>% mask(roads_study_area)
  ken_viirs_2013_median <- raster(file.path(raw_data_file_path,"NTL Rasters Annual","ken_viirs_2013_median.tif")) %>% crop(roads_study_area) %>% mask(roads_study_area)
  ken_viirs_2014_median <- raster(file.path(raw_data_file_path,"NTL Rasters Annual","ken_viirs_2014_median.tif")) %>% crop(roads_study_area) %>% mask(roads_study_area)
  ken_viirs_2015_median <- raster(file.path(raw_data_file_path,"NTL Rasters Annual","ken_viirs_2015_median.tif")) %>% crop(roads_study_area) %>% mask(roads_study_area)
  ken_viirs_2016_median <- raster(file.path(raw_data_file_path,"NTL Rasters Annual","ken_viirs_2016_median.tif")) %>% crop(roads_study_area) %>% mask(roads_study_area)
  ken_viirs_2017_median <- raster(file.path(raw_data_file_path,"NTL Rasters Annual","ken_viirs_2017_median.tif")) %>% crop(roads_study_area) %>% mask(roads_study_area)
  ken_viirs_2018_median <- raster(file.path(raw_data_file_path,"NTL Rasters Annual","ken_viirs_2018_median.tif")) %>% crop(roads_study_area) %>% mask(roads_study_area)
  
  # Create VIIRS Polygon
  viir_sdf <- rasterToPolygons(ken_viirs_2012_median)
  viir_sdf$ken_viirs_2013_median <- ken_viirs_2013_median[][!is.na(ken_viirs_2013_median[])]
  viir_sdf$ken_viirs_2014_median <- ken_viirs_2014_median[][!is.na(ken_viirs_2014_median[])]
  viir_sdf$ken_viirs_2015_median <- ken_viirs_2015_median[][!is.na(ken_viirs_2015_median[])]
  viir_sdf$ken_viirs_2016_median <- ken_viirs_2016_median[][!is.na(ken_viirs_2016_median[])]
  viir_sdf$ken_viirs_2017_median <- ken_viirs_2017_median[][!is.na(ken_viirs_2017_median[])]
  viir_sdf$ken_viirs_2018_median <- ken_viirs_2018_median[][!is.na(ken_viirs_2018_median[])]
  
  #### Extract Africover and Globcover to Polygon
  # Africover/Globcover Urban Raster
  africover_urban <- calc(africover,  fun=convert_africover_builtup)
  globcover_urban <- calc(globcover_2015,  fun=convert_globcover_builtup)
  
  # Extract
  viir_sdf$africover_urban <- NA
  viir_sdf$globcover_urban <- NA
  chunk_size <- 2000
  start_ids <- seq(1, nrow(viir_sdf), by=chunk_size)
  for(start_i in start_ids){
    end_i <- min(start_i + chunk_size - 1, nrow(viir_sdf))
    
    num_range <- start_i:end_i
    viir_sdf_i <- viir_sdf[num_range,]
    africover_urban_i <- africover_urban %>% crop(viir_sdf_i)
    globcover_urban_i <- globcover_urban %>% crop(viir_sdf_i)
    
    viir_sdf$africover_urban[num_range] <- as.numeric(velox(africover_urban_i)$extract(sp=viir_sdf_i, fun=function(x){mean(x, na.rm=TRUE)}))
    viir_sdf$globcover_urban[num_range] <- as.numeric(velox(globcover_urban_i)$extract(sp=viir_sdf_i, fun=function(x){mean(x, na.rm=TRUE)}))
    
    print(start_i)
  }
  
  #### Study Area
  viir_sdf_OVER_roads_study_area <- over(viir_sdf, roads_study_area)
  viir_sdf$road <- viir_sdf_OVER_roads_study_area$road
  
  save(viir_sdf, file=file.path(intermediate_data_file_path, "viirs_shapefile", "viirslevel_studyarea_sdf.Rda"))
} else{
  load(file.path(intermediate_data_file_path, "viirs_shapefile", "viirslevel_studyarea_sdf.Rda"))
}

# Construct Threshold Dummies --------------------------------------------------
for(threshold in c(0.5,1,2,5)){
  for(year in 2012:2018){
    viir_sdf[[paste0("ntl",year,"_thresh",threshold)]] <- as.numeric(viir_sdf[[paste0("ken_viirs_",year,"_median")]] >= threshold) * 750 / 1000
  }
}

# NTL Thresholds Over Time -----------------------------------------------------
# Collapse
viir_sdf_collapsed <- summaryBy(. ~ road, data=viir_sdf@data, FUN=sum, keep.names=T)

var_names <- names(viir_sdf_collapsed)[grepl("2012", names(viir_sdf_collapsed))]
var_names <- gsub("2012","",var_names)
var_names <- gsub("__","_",var_names)

subset_data_for_stacking <- function(year){
  df_out <- viir_sdf_collapsed[,names(viir_sdf_collapsed)[grepl(paste0(year), names(viir_sdf_collapsed))]]
  names(df_out) <- var_names
  df_out$year <- year
  df_out$road <- viir_sdf_collapsed$road 
  return(df_out)
}

viirs_stacked <- lapply(2012:2017, subset_data_for_stacking) %>% bind_rows

viirs_stacked$road <- as.character(viirs_stacked$road)
viirs_stacked$road[viirs_stacked$road %in% "a1"] <- "A1"
viirs_stacked$road[viirs_stacked$road %in% "a2"] <- "A2"
viirs_stacked$road[viirs_stacked$road %in% "b9"] <- "B9"

ntl_urban_thrs0.5 <- ggplot(viirs_stacked) + 
  geom_line(aes(x=year, y=ntl_thresh0.5, color=road),size=1.5) + 
  labs(x="",
       y="Urban Area (Kilometers Squared)",
       color="Road",
       title="Radiance > 0.5 Considered Urban") +
  scale_y_continuous(limits = c(0, max(viirs_stacked$ntl_thresh0.5))) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) 

ntl_urban_thrs1 <- ggplot(viirs_stacked) + 
  geom_line(aes(x=year, y=ntl_thresh1, color=road),size=1.5) + 
  labs(x="",
       y="Urban Area (Kilometers Squared)",
       color="Road",
       title="Radiance > 1 Considered Urban") +
  scale_y_continuous(limits = c(0, max(viirs_stacked$ntl_thresh1))) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) 

ntl_urban_thrs2 <- ggplot(viirs_stacked) + 
  geom_line(aes(x=year, y=ntl_thresh2, color=road),size=1.5) + 
  labs(x="",
       y="Urban Area (Kilometers Squared)",
       color="Road",
       title="Radiance > 2 Considered Urban") +
  scale_y_continuous(limits = c(0, max(viirs_stacked$ntl_thresh2))) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) 

ntl_urban_thrs5 <- ggplot(viirs_stacked) + 
  geom_line(aes(x=year, y=ntl_thresh5, color=road),size=1.5) + 
  labs(x="",
       y="Urban Area (Kilometers Squared)",
       color="Road",
       title="Radiance > 5 Considered Urban") +
  scale_y_continuous(limits = c(0, max(viirs_stacked$ntl_thresh5))) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) 

ntl_urban_thrs <- ggarrange(ntl_urban_thrs0.5, ntl_urban_thrs1, ntl_urban_thrs2, ntl_urban_thrs5, ncol=2, nrow=2, common.legend = TRUE, legend="right")
ggsave(ntl_urban_thrs, filename=file.path(figures_file_path, "viirs_threshold_urban_trends.png"),height=8,width=9)

# Does VIIRS Exclude Urban Areas: Africover ------------------------------------
determine_ntl_threshold_accuracy <- function(ntl_threshold, true_urban_threshold){
  pixels_with_true_urban_area <- sum(viir_sdf$africover_urban > true_urban_threshold)
  
  correctly_identified_urban <- sum((viir_sdf$africover_urban > true_urban_threshold) & (viir_sdf$ken_viirs_2016_median >= ntl_threshold))
  incorrectly_identified_urban <- sum((viir_sdf$africover_urban <= true_urban_threshold) & (viir_sdf$ken_viirs_2016_median >= ntl_threshold))
  
  df_out <- c(correctly_identified_urban, incorrectly_identified_urban) %>% t %>% as.data.frame
  names(df_out) <- c("correctly_identified_urban", "incorrectly_identified_urban")
  
  df_out$ntl_threshold <- ntl_threshold 
  df_out$pixels_with_true_urban_area <- pixels_with_true_urban_area 
  return(df_out)
}

accuracy_df <- lapply(seq(from=.15,to=5,by=.05), determine_ntl_threshold_accuracy, 0) %>% bind_rows
accuracy_africover_0 <- ggplot(accuracy_df) + 
  geom_line(aes(x=ntl_threshold, y=correctly_identified_urban, color="Correctly Identified as Urban"),size=1) +
  geom_line(aes(x=ntl_threshold, y=incorrectly_identified_urban, color="Incorrectly Identified as Urban"),size=1) +
  geom_hline(yintercept=accuracy_df$pixels_with_true_urban_area[1], color="black", size=1) +
  geom_text(label="Urban Pixels, According to Africover",
            x=2.8,y=accuracy_df$pixels_with_true_urban_area[1]-60) +
  scale_color_manual(values=c("blue","red")) +
  theme_minimal() +
  labs(x="Nighttime Lights Threshold for Urban", 
       y="Number of Pixels",
       color="",
       title="A. Defining Pixels with > 0% Urban Area (According to Africover) as Urban")

accuracy_df <- lapply(seq(from=.15,to=5,by=.05), determine_ntl_threshold_accuracy, .1) %>% bind_rows
accuracy_africover_0.1 <- ggplot(accuracy_df) + 
  geom_line(aes(x=ntl_threshold, y=correctly_identified_urban, color="Correctly Identified as Urban"),size=1) +
  geom_line(aes(x=ntl_threshold, y=incorrectly_identified_urban, color="Incorrectly Identified as Urban"),size=1) +
  geom_hline(yintercept=accuracy_df$pixels_with_true_urban_area[1], color="black", size=1) +
  geom_text(label="Urban Pixels, According to Africover",
            x=2.8,y=accuracy_df$pixels_with_true_urban_area[1]+70) +
  scale_color_manual(values=c("blue","red")) +
  theme_minimal() +
  labs(x="Nighttime Lights Threshold for Urban", 
       y="Number of Pixels",
       color="",
       title="B. Defining Pixels with > 10% Urban Area (According to Africover) as Urban")

accuracy_africover <- grid.arrange(accuracy_africover_0, accuracy_africover_0.1, nrow = 2, ncol=1)
ggsave(accuracy_africover, filename=file.path(figures_file_path, "accuracy_ntlthresh_africover.png"),height=10,width=9)

# Does VIIRS Exclude Urban Areas: Globcover ------------------------------------
determine_ntl_threshold_accuracy <- function(ntl_threshold, true_urban_threshold){
  pixels_with_true_urban_area <- sum(viir_sdf$globcover_urban > true_urban_threshold)
  
  correctly_identified_urban <- sum((viir_sdf$globcover_urban > true_urban_threshold) & (viir_sdf$ken_viirs_2015_median >= ntl_threshold))
  incorrectly_identified_urban <- sum((viir_sdf$globcover_urban <= true_urban_threshold) & (viir_sdf$ken_viirs_2015_median >= ntl_threshold))
  
  df_out <- c(correctly_identified_urban, incorrectly_identified_urban) %>% t %>% as.data.frame
  names(df_out) <- c("correctly_identified_urban", "incorrectly_identified_urban")
  
  df_out$ntl_threshold <- ntl_threshold 
  df_out$pixels_with_true_urban_area <- pixels_with_true_urban_area 
  return(df_out)
}

accuracy_df <- lapply(seq(from=.25,to=5,by=.05), determine_ntl_threshold_accuracy, 0) %>% bind_rows
accuracy_globcover_0 <- ggplot(accuracy_df) + 
  geom_line(aes(x=ntl_threshold, y=correctly_identified_urban, color="Correctly Identified as Urban"),size=1) +
  geom_line(aes(x=ntl_threshold, y=incorrectly_identified_urban, color="Incorrectly Identified as Urban"),size=1) +
  geom_hline(yintercept=accuracy_df$pixels_with_true_urban_area[1], color="black", size=1) +
  geom_text(label="Urban Pixels, According to Globcover",
            x=1.5,y=accuracy_df$pixels_with_true_urban_area[1]+30) +
  scale_color_manual(values=c("blue","red")) +
  theme_minimal() +
  labs(x="Nighttime Lights Threshold for Urban", y="Number of Pixels",color="")
ggsave(accuracy_globcover_0, filename=file.path(figures_file_path, "accuracy_ntlthresh_globcover.png"),height=5,width=9)
