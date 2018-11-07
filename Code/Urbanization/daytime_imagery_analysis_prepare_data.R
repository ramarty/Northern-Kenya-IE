# Northern Kenya IE

# Setup ------------------------------------------------------------------------
if(Sys.info()["user"] == "WB521633") source("C:/Users/wb521633/Dropbox/World Bank/IEs/Northern Kenya IE/Code/northern_kenya_ie_master.R")
if(Sys.info()["user"] == "robmarty") source("~/Dropbox/World Bank/IEs/Northern Kenya IE/Code/northern_kenya_ie_master.R")

library(raster)
library(rgdal)
library(velox)
library(doBy)
library(ggplot2)
library(gridExtra)

buffer_size <- 20

set.seed(42)

# Load Data --------------------------------------------------------------------
  
# Study Area
setwd(file.path(intermediate_data_file_path, "Treatment Roads Buffer"))
roads_study_area <- readOGR(dsn=".", layer=paste0("kenya_treat_roads_",buffer_size,"km_buff"))
roads_study_area$road <- roads_study_area$road %>% as.character

# Africover  
africover <- raster(file.path(intermediate_data_file_path, "Africover", "africover_kenya.tif")) 
africover_a1 <- raster(file.path(intermediate_data_file_path, "Africover", "africover_kenya.tif")) %>% raster::crop(roads_study_area[roads_study_area$road %in% "a1",])
africover_a2 <- raster(file.path(intermediate_data_file_path, "Africover", "africover_kenya.tif")) %>% crop(roads_study_area[roads_study_area$road %in% "a2",])
africover_b9 <- raster(file.path(intermediate_data_file_path, "Africover", "africover_kenya.tif")) %>% crop(roads_study_area[roads_study_area$road %in% "b9",])

# Create Grids to Random Sample Points In --------------------------------------

# Grid in Study Area
grid_within_shp <- function(roads_study_area){
  grdpts <- makegrid(roads_study_area, cellsize = .1)
  spgrd <- SpatialPoints(grdpts, proj4string = CRS(proj4string(roads_study_area)))
  spgrdWithin <- SpatialPixels(spgrd[roads_study_area,])
  spgrdWithin <- as(spgrdWithin, "SpatialPolygons")
  return(spgrdWithin)
}

africover_a1_grid <- grid_within_shp(roads_study_area[roads_study_area$road %in% "a1",])
africover_a2_grid <- grid_within_shp(roads_study_area[roads_study_area$road %in% "a2",])
africover_b9_grid <- grid_within_shp(roads_study_area[roads_study_area$road %in% "b9",])

# Create Random Points ---------------------------------------------------------
random_sample_grid <- function(i, africover, africover_grid, road, stratefied){
  
  is_error <- tryCatch({
    
    africover_points_i <- africover %>% 
      crop(africover_grid[i,]) %>%
      #mask(roads_study_area[roads_study_area$road %in% road,]) %>%
      rasterToPoints %>%
      as.data.frame
    
  },
  error = function(e) return("error"))
  
  if(!is.character(is_error)){
    
    # Maybe a way to use velox // over to further subset
    
    sample_dataframe <- function(df, var, i, sample){
      df_i <- df[df[[var]] %in% i,]
      
      if(nrow(df_i) <= sample){
        return(df_i)
      } else{
        df_i <- df_i %>% sample_n(sample)
        return(df_i)
      }
    }
    
    sample_n_mult <- 15
    
    if(stratefied){
    df_out <- rbind(
      sample_dataframe(africover_points_i, "africover_kenya", 0, 5*sample_n_mult),
      sample_dataframe(africover_points_i, "africover_kenya", 1, 70*sample_n_mult),
      sample_dataframe(africover_points_i, "africover_kenya", 2, 300*sample_n_mult),
      sample_dataframe(africover_points_i, "africover_kenya", 3, 500*sample_n_mult),
      sample_dataframe(africover_points_i, "africover_kenya", 4, 120*sample_n_mult),
      sample_dataframe(africover_points_i, "africover_kenya", 5, 5*sample_n_mult),
      sample_dataframe(africover_points_i, "africover_kenya", 6, 25*sample_n_mult),
      sample_dataframe(africover_points_i, "africover_kenya", 7, 50*sample_n_mult),
      africover_points_i[africover_points_i$africover_kenya %in% 8,],
      sample_dataframe(africover_points_i, "africover_kenya", 9, 5*sample_n_mult),
      sample_dataframe(africover_points_i, "africover_kenya", 10, 5*sample_n_mult)
    )
    } else{
      df_out <- rbind(
        sample_dataframe(africover_points_i, "africover_kenya", i=c(0:7,9:10), 15000),
        africover_points_i[africover_points_i$africover_kenya %in% 8,]
      )
      
    }
    
    df_out$road <- road
    
    
  } else{
    df_out <- NULL %>% as.data.frame
  }
  
  print(i)
  return(df_out)
}

random_points_a1 <- lapply(1:length(africover_a1_grid), random_sample_grid, africover_a1, africover_a1_grid, "a1", F) %>% bind_rows
random_points_a2 <- lapply(1:length(africover_a2_grid), random_sample_grid, africover_a2, africover_a2_grid, "a2", F) %>% bind_rows
random_points_b9 <- lapply(1:length(africover_b9_grid), random_sample_grid, africover_b9, africover_b9_grid, "b9", F) %>% bind_rows

random_points <- rbind(random_points_a1, random_points_a2, random_points_b9)
coordinates(random_points) <- ~x+y
crs(random_points) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

random_points <- random_points[!is.na(over(random_points, roads_study_area)$id),]

# Extract Landsat Values to Random Points --------------------------------------
for(band in 1:12){
  
  i=1
  for(raster_file in c("kenya_road_1_20km_landsat_2016_median-0000000000-0000000000.tif",
                       "kenya_road_1_20km_landsat_2016_median-0000009472-0000000000.tif",
                       "kenya_road_2_20km_landsat_2016_median-0000009472-0000000000.tif",
                       "kenya_road_2_20km_landsat_2016_median-0000000000-0000000000.tif",
                       "kenya_road_3_20km_landsat_2016_median-0000000000-0000000000.tif",
                       "kenya_road_3_20km_landsat_2016_median-0000000000-0000009472.tif",
                       "kenya_road_3_20km_landsat_2016_median-0000009472-0000000000.tif",
                       "kenya_road_3_20km_landsat_2016_median-0000009472-0000009472.tif")){
    
    landsat_i <- raster(file.path(raw_data_file_path, "Landsat", raster_file), band=band)
    random_points[[paste0("b",band,"_",i)]] <- extract(landsat_i, random_points)
    i <- i + 1
    print("next")
  }
  
  band_vars_df <- random_points@data[,names(random_points)[grepl(paste0("b",band,"_"),names(random_points))]]
  random_points[[paste0("band_",band)]] <- rowSums(band_vars_df, na.rm=TRUE) * ifelse(rowSums(is.na(band_vars_df)) == ncol(band_vars_df), NA, 1)
  random_points@data <- subset(random_points@data, select=c("africover_kenya", 
                                                            "road",
                                                            names(random_points)[grepl("band_",names(random_points))]))
  
  print(band)
  print(names(random_points))
}

# Extract Nighttime Lights -----------------------------------------------------
years <- c(rep(2012, length(4:12)), rep(2013,12), rep(2014,12), rep(2015,12), rep(2016,12), rep(2017,12), rep(2018,12))
months <- c(4:12                  , 1:12        , 1:12        , 1:12        , 1:12        , 1:12        , 1:12)
viirs_key_df <- cbind(years, months) %>% as.data.frame
viirs_key_df$id <- 1:nrow(viirs_key_df)

viirs_avg_rad_1 <- raster(file.path(raw_data_file_path, "viirs_monthly", "kenya_viirs_2012_to_2018_avg_rad.tif"), band=46) %>% crop(roads_study_area)
viirs_avg_rad_2 <- raster(file.path(raw_data_file_path, "viirs_monthly", "kenya_viirs_2012_to_2018_avg_rad.tif"), band=47) %>% crop(roads_study_area)
viirs_avg_rad_3 <- raster(file.path(raw_data_file_path, "viirs_monthly", "kenya_viirs_2012_to_2018_avg_rad.tif"), band=48) %>% crop(roads_study_area)
viirs_avg_rad_4 <- raster(file.path(raw_data_file_path, "viirs_monthly", "kenya_viirs_2012_to_2018_avg_rad.tif"), band=49) %>% crop(roads_study_area)
viirs_avg_rad_5 <- raster(file.path(raw_data_file_path, "viirs_monthly", "kenya_viirs_2012_to_2018_avg_rad.tif"), band=50) %>% crop(roads_study_area)
viirs_avg_rad_6 <- raster(file.path(raw_data_file_path, "viirs_monthly", "kenya_viirs_2012_to_2018_avg_rad.tif"), band=51) %>% crop(roads_study_area)
viirs_avg_rad_7 <- raster(file.path(raw_data_file_path, "viirs_monthly", "kenya_viirs_2012_to_2018_avg_rad.tif"), band=52) %>% crop(roads_study_area)
viirs_avg_rad_8 <- raster(file.path(raw_data_file_path, "viirs_monthly", "kenya_viirs_2012_to_2018_avg_rad.tif"), band=53) %>% crop(roads_study_area)
viirs_avg_rad_9 <- raster(file.path(raw_data_file_path, "viirs_monthly", "kenya_viirs_2012_to_2018_avg_rad.tif"), band=54) %>% crop(roads_study_area)
viirs_avg_rad_10 <- raster(file.path(raw_data_file_path, "viirs_monthly", "kenya_viirs_2012_to_2018_avg_rad.tif"), band=55) %>% crop(roads_study_area)
viirs_avg_rad_11 <- raster(file.path(raw_data_file_path, "viirs_monthly", "kenya_viirs_2012_to_2018_avg_rad.tif"), band=56) %>% crop(roads_study_area)
viirs_avg_rad_12 <- raster(file.path(raw_data_file_path, "viirs_monthly", "kenya_viirs_2012_to_2018_avg_rad.tif"), band=57) %>% crop(roads_study_area)

viirs_cf_cvg_1 <- raster(file.path(raw_data_file_path, "viirs_monthly", "kenya_viirs_2012_to_2018_cf_cvg.tif"), band=46) %>% crop(roads_study_area)
viirs_cf_cvg_2 <- raster(file.path(raw_data_file_path, "viirs_monthly", "kenya_viirs_2012_to_2018_cf_cvg.tif"), band=47) %>% crop(roads_study_area)
viirs_cf_cvg_3 <- raster(file.path(raw_data_file_path, "viirs_monthly", "kenya_viirs_2012_to_2018_cf_cvg.tif"), band=48) %>% crop(roads_study_area)
viirs_cf_cvg_4 <- raster(file.path(raw_data_file_path, "viirs_monthly", "kenya_viirs_2012_to_2018_cf_cvg.tif"), band=49) %>% crop(roads_study_area)
viirs_cf_cvg_5 <- raster(file.path(raw_data_file_path, "viirs_monthly", "kenya_viirs_2012_to_2018_cf_cvg.tif"), band=50) %>% crop(roads_study_area)
viirs_cf_cvg_6 <- raster(file.path(raw_data_file_path, "viirs_monthly", "kenya_viirs_2012_to_2018_cf_cvg.tif"), band=51) %>% crop(roads_study_area)
viirs_cf_cvg_7 <- raster(file.path(raw_data_file_path, "viirs_monthly", "kenya_viirs_2012_to_2018_cf_cvg.tif"), band=52) %>% crop(roads_study_area)
viirs_cf_cvg_8 <- raster(file.path(raw_data_file_path, "viirs_monthly", "kenya_viirs_2012_to_2018_cf_cvg.tif"), band=53) %>% crop(roads_study_area)
viirs_cf_cvg_9 <- raster(file.path(raw_data_file_path, "viirs_monthly", "kenya_viirs_2012_to_2018_cf_cvg.tif"), band=54) %>% crop(roads_study_area)
viirs_cf_cvg_10 <- raster(file.path(raw_data_file_path, "viirs_monthly", "kenya_viirs_2012_to_2018_cf_cvg.tif"), band=55) %>% crop(roads_study_area)
viirs_cf_cvg_11 <- raster(file.path(raw_data_file_path, "viirs_monthly", "kenya_viirs_2012_to_2018_cf_cvg.tif"), band=56) %>% crop(roads_study_area)
viirs_cf_cvg_12 <- raster(file.path(raw_data_file_path, "viirs_monthly", "kenya_viirs_2012_to_2018_cf_cvg.tif"), band=57) %>% crop(roads_study_area)

random_points$viirs_rad_1 <- extract(viirs_avg_rad_1, random_points)
random_points$viirs_rad_2 <- extract(viirs_avg_rad_2, random_points)
random_points$viirs_rad_3 <- extract(viirs_avg_rad_3, random_points)
random_points$viirs_rad_4 <- extract(viirs_avg_rad_4, random_points)
random_points$viirs_rad_5 <- extract(viirs_avg_rad_5, random_points)
random_points$viirs_rad_6 <- extract(viirs_avg_rad_6, random_points)
random_points$viirs_rad_7 <- extract(viirs_avg_rad_7, random_points)
random_points$viirs_rad_8 <- extract(viirs_avg_rad_8, random_points)
random_points$viirs_rad_9 <- extract(viirs_avg_rad_9, random_points)
random_points$viirs_rad_10 <- extract(viirs_avg_rad_10, random_points)
random_points$viirs_rad_11 <- extract(viirs_avg_rad_11, random_points)
random_points$viirs_rad_12 <- extract(viirs_avg_rad_12, random_points)

random_points$viirs_cf_cvg_1 <- extract(viirs_cf_cvg_1, random_points)
random_points$viirs_cf_cvg_2 <- extract(viirs_cf_cvg_2, random_points)
random_points$viirs_cf_cvg_3 <- extract(viirs_cf_cvg_3, random_points)
random_points$viirs_cf_cvg_4 <- extract(viirs_cf_cvg_4, random_points)
random_points$viirs_cf_cvg_5 <- extract(viirs_cf_cvg_5, random_points)
random_points$viirs_cf_cvg_6 <- extract(viirs_cf_cvg_6, random_points)
random_points$viirs_cf_cvg_7 <- extract(viirs_cf_cvg_7, random_points)
random_points$viirs_cf_cvg_8 <- extract(viirs_cf_cvg_8, random_points)
random_points$viirs_cf_cvg_9 <- extract(viirs_cf_cvg_9, random_points)
random_points$viirs_cf_cvg_10 <- extract(viirs_cf_cvg_10, random_points)
random_points$viirs_cf_cvg_11 <- extract(viirs_cf_cvg_11, random_points)
random_points$viirs_cf_cvg_12 <- extract(viirs_cf_cvg_12, random_points)

# Save Data --------------------------------------------------------------------
# Convert to Dataframe
random_points$coord_x <- coordinates(random_points)[,1]
random_points$coord_y <- coordinates(random_points)[,2]
random_points <- random_points@data

save(random_points, file=file.path(intermediate_data_file_path, "Africover Random Points with Landsat", "random_points_landsat.Rda"))

# Proportaion Classes ----------------------------------------------------------
if(F){
  random_points_a1 <- lapply(1:length(africover_a1_grid), random_sample_grid, africover_a1, africover_a1_grid, "a1", FALSE) %>% bind_rows
  random_points_a2 <- lapply(1:length(africover_a2_grid), random_sample_grid, africover_a2, africover_a2_grid, "a2", FALSE) %>% bind_rows
  random_points_b9 <- lapply(1:length(africover_b9_grid), random_sample_grid, africover_b9, africover_b9_grid, "b9", FALSE) %>% bind_rows
  
  random_points <- rbind(random_points_a1, random_points_a2, random_points_b9)
  coordinates(random_points) <- ~x+y
  crs(random_points) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  
  random_points <- random_points[!is.na(over(random_points, roads_study_area)$id),]
  
  class_proportions <- random_points$africover_kenya %>% 
    table %>% 
    as.data.frame %>% 
    rename(class=".")
  class_proportions$Proportion <- class_proportions$Freq / sum(class_proportions$Freq)
  write.csv(class_proportions, file.path(intermediate_data_file_path, "Africover Random Points with Landsat", "africover_class_proportions.csv"), row.names=F)
}