# Create Urban Rasters from VIIRS
# Northern Kenya IE

# https://pdfs.semanticscholar.org/cf1a/5c9a63f1a5cd0bf4513fa84a3dfc14b43f60.pdf
# Gradient Magnitude: https://cran.r-project.org/web/packages/imager/vignettes/gettingstarted.html
# https://www.bioconductor.org/help/course-materials/2015/BioC2015/BioC2015Oles.html#38

# Setup ------------------------------------------------------------------------
if(Sys.info()["user"] == "WB521633") source("C:/Users/wb521633/Dropbox/World Bank/IEs/Northern Kenya IE/Code/northern_kenya_ie_master.R")
if(Sys.info()["user"] == "robmarty") source("~/Dropbox/World Bank/IEs/Northern Kenya IE/Code/northern_kenya_ie_master.R")

import_raw_africover <- FALSE
create_viirs_shapefile_raw <- FALSE

# Load Data --------------------------------------------------------------------
# Nighttime Lights
viirs_band1 <- raster(file.path(project_file_path, "Data", "RawData","viirs_monthly","kenya_viirs_2012_to_2018_avg_rad.tif"), band=1)

# VIIRS Shapefile
load(file.path(intermediate_data_file_path, "viirs_shapefile", "viirs_band1_shp_africover_urban.Rda"))
viirs_level_df <- viirs_band1_shp_africover_urban

# GADM
ken_adm0 <- getData('GADM', country='KEN', level=0)

# Nairobi
nairobi <- c(-1.283333, 36.816667) %>% t %>% as.data.frame %>% dplyr::rename(lat=V1) %>% dplyr::rename(lon=V2)
viirs_band1_shp_africover_urban$distance_nairobi <- sqrt((viirs_band1_shp_africover_urban$longitude - nairobi$lon)^2 + (viirs_band1_shp_africover_urban$latitude - nairobi$lat)^2) * 111.12

# Annual VIIRS -----------------------------------------------------------------
year <- c(rep(2012, 9), rep(2013, 12), rep(2014, 12), rep(2015, 12), rep(2016, 12), rep(2017, 12), rep(2018, 6))
month <- c(4:9,          1:12,          1:12,          1:12,          1:12,          1:12,          1:6)
id <- 1:75
raster_ids <- cbind(id, year, month) %>% as.data.frame
raster_ids$id[raster_ids$year %in% 2016]

viirs_rad_brick <- stack(file.path(project_file_path, "Data", "RawData","viirs_monthly","kenya_viirs_2012_to_2018_avg_rad.tif")) %>% subset(46:57) %>% brick() %>% crop(ken_adm0) %>% mask(ken_adm0)
viirs_cf_brick <- stack(file.path(project_file_path, "Data", "RawData","viirs_monthly","kenya_viirs_2012_to_2018_cf_cvg.tif")) %>% subset(46:57) %>% brick() %>% crop(ken_adm0) %>% mask(ken_adm0)

# Needs to have more than 1 cloud free day; or else, NA radiance
viirs_rad_brick$kenya_viirs_2012_to_2018_avg_rad.46[][viirs_cf_brick$kenya_viirs_2012_to_2018_cf_cvg.46[] <= 1] <- NA
viirs_rad_brick$kenya_viirs_2012_to_2018_avg_rad.47[][viirs_cf_brick$kenya_viirs_2012_to_2018_cf_cvg.47[] <= 1] <- NA
viirs_rad_brick$kenya_viirs_2012_to_2018_avg_rad.48[][viirs_cf_brick$kenya_viirs_2012_to_2018_cf_cvg.48[] <= 1] <- NA
viirs_rad_brick$kenya_viirs_2012_to_2018_avg_rad.49[][viirs_cf_brick$kenya_viirs_2012_to_2018_cf_cvg.49[] <= 1] <- NA
viirs_rad_brick$kenya_viirs_2012_to_2018_avg_rad.50[][viirs_cf_brick$kenya_viirs_2012_to_2018_cf_cvg.50[] <= 1] <- NA
viirs_rad_brick$kenya_viirs_2012_to_2018_avg_rad.51[][viirs_cf_brick$kenya_viirs_2012_to_2018_cf_cvg.51[] <= 1] <- NA
viirs_rad_brick$kenya_viirs_2012_to_2018_avg_rad.52[][viirs_cf_brick$kenya_viirs_2012_to_2018_cf_cvg.52[] <= 1] <- NA
viirs_rad_brick$kenya_viirs_2012_to_2018_avg_rad.53[][viirs_cf_brick$kenya_viirs_2012_to_2018_cf_cvg.53[] <= 1] <- NA
viirs_rad_brick$kenya_viirs_2012_to_2018_avg_rad.54[][viirs_cf_brick$kenya_viirs_2012_to_2018_cf_cvg.54[] <= 1] <- NA
viirs_rad_brick$kenya_viirs_2012_to_2018_avg_rad.55[][viirs_cf_brick$kenya_viirs_2012_to_2018_cf_cvg.55[] <= 1] <- NA
viirs_rad_brick$kenya_viirs_2012_to_2018_avg_rad.56[][viirs_cf_brick$kenya_viirs_2012_to_2018_cf_cvg.56[] <= 1] <- NA
viirs_rad_brick$kenya_viirs_2012_to_2018_avg_rad.57[][viirs_cf_brick$kenya_viirs_2012_to_2018_cf_cvg.57[] <= 1] <- NA

# Collapse raster brick
viirs_rad <- calc(viirs_rad_brick, fun = median, na.rm = T)

# Cut off: Areas below here cannot be urban.
viirs_rad[][viirs_rad[] <= 0.2] <- 0

# Add VIIRS to viirs level dataframe
viirs_level_df$viirs_rad <- viirs_rad[][!is.na(viirs_rad[])]

# Potential Urban Clusters -----------------------------------------------------
calc_raster_gradient_mag <- function(r){
  g <- rast.grad(r)
  r_gm <- sqrt(g$rast.grad.x^2 + g$rast.grad.y^2)
  return(r_gm)
}

# Segmentation
viirs_rad_cluster <- viirs_rad %>% 
                      calc_raster_gradient_mag %>%
                      as.cimg %>% 
                      distmap %>% 
                      watershed(tolerance=0.3,ext=5)
#plot(colorLabels(getFrame(viirs_rad_segment, 1)))

# Raster of Segmentation Values
viirs_rad_cluster_r <- viirs_rad # Make raster to replace values of
viirs_rad_cluster_r[] <- viirs_rad_cluster[] %>% as.list %>% unlist

# Extract Raster Segment IDs to VIIRS DF Grid
viirs_level_sdf <- viirs_level_df@data
coordinates(viirs_level_sdf) <- ~longitude+latitude
crs(viirs_level_sdf) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

viirs_level_sdf$cluster_id <- extract(viirs_rad_cluster_r, viirs_level_sdf)

# Determine Optimal Threshold --------------------------------------------------
viirs_level_dt <- as.data.table(viirs_level_sdf@data)
viirs_level_dt$Ncells <- 1
viirs_level_dt$seg_size <- 750

create_cluster_level_dataset <- function(cluster_id, dt){
  dt_i <- dt[dt$cluster_id %in% cluster_id,]
  proportion_urban_truth <- mean(dt_i$africover_urban)
  
  potential_thresholds <- seq(0,2,by=.02)
  proportion_urban <- lapply(potential_thresholds, function(thresh) mean(dt_i$viirs_rad > thresh)) %>% unlist

  threshold_urban_df <- cbind(potential_thresholds, proportion_urban)
  optm_thresh_df <- threshold_urban_df[which.min(abs(threshold_urban_df[,2] - proportion_urban_truth)),] %>% t %>% as.data.frame
  
  names(optm_thresh_df) <- c("ntl_threshold","proportion_urban_ntlthresh")
  
  # Add other variables
  optm_thresh_df$cluster_id <- cluster_id
  optm_thresh_df$Ncells <- sum(dt_i$Ncells)
  optm_thresh_df$seg_size <- sum(dt_i$seg_size)
  optm_thresh_df$proportion_urban_truth <- proportion_urban_truth
  optm_thresh_df$viirs_rad_mean <- mean(dt_i$viirs_rad)
  optm_thresh_df$viirs_rad_min <- min(dt_i$viirs_rad)
  optm_thresh_df$viirs_rad_max <- max(dt_i$viirs_rad)
  
  return(optm_thresh_df)
}

data_mean_sum <- pbmclapply(sort(unique(viirs_level_dt$cluster_id)), create_cluster_level_dataset, viirs_level_dt, mc.cores=1) %>% bind_rows

data_mean_sum <- data_mean_sum[data_mean_sum$cluster_id != 0,]
data_mean_sum <- data_mean_sum[data_mean_sum$proportion_urban_thresh > 0,]

a<-.4
b<-.2
data_mean_sum$x_indicator <- log(data_mean_sum$seg_size^a * data_mean_sum$viirs_rad_mean^b)

plot(data_mean_sum$x_indicator, data_mean_sum$proportion_urban_thresh)


plot(data_median_sum$viirs_rad_mean, data_median_sum$proportion_urban_thresh)
plot(data_median_sum$seg_size, data_median_sum$proportion_urban_thresh)

data_mean_sum$proportion_urban_thresh



determine_optimal_threshold_per_cluster(1,data.dt)

data_mean_sum <- data.dt[,list(africover_urban_avg=mean(africover_urban), 
                               viirs_rad_mean_avg=mean(viirs_rad_mean), 
                               seg_size=sum(seg_size),
                               Ncells=sum(Ncells)),by=segment_value_mean]
data_mean_sum <- data_mean_sum[data_mean_sum$segment_value_mean != 0,]

data_median_sum <- data.dt[,list(africover_urban_avg=mean(africover_urban), 
                                 viirs_rad_median_avg=mean(viirs_rad_median), 
                                 seg_size=sum(seg_size),
                                 Ncells=sum(Ncells)),by=segment_value_median]
data_median_sum <- data_median_sum[data_median_sum$segment_value_median != 0,]

a = .2
b= .4
x <- log((data_median_sum$seg_size^a) * (data_median_sum$viirs_rad_mean_avg^b))
