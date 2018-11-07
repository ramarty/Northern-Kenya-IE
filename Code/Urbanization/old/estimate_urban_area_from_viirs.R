# Create Urban Rasters from VIIRS
# Northern Kenya IE

# Useful Resources
# https://pdfs.semanticscholar.org/cf1a/5c9a63f1a5cd0bf4513fa84a3dfc14b43f60.pdf
# Gradient Magnitude: https://cran.r-project.org/web/packages/imager/vignettes/gettingstarted.html
# https://www.bioconductor.org/help/course-materials/2015/BioC2015/BioC2015Oles.html#38

# Setup ------------------------------------------------------------------------
if(Sys.info()["user"] == "WB521633") source("C:/Users/wb521633/Dropbox/World Bank/IEs/Northern Kenya IE/Code/northern_kenya_ie_master.R")
if(Sys.info()["user"] == "robmarty") source("~/Dropbox/World Bank/IEs/Northern Kenya IE/Code/northern_kenya_ie_master.R")

import_raw_africover <- FALSE
create_viirs_shapefile_raw <- FALSE

VIIRS_YEAR_MONTH_START <- "2016-01"
VIIRS_YEAR_MONTH_END <- "2016-12"
VIIRS_DNB_PIXEL_MIN <- 0.35
VIIRS_NUMBER_CLOUDFREE_PIXELS <- 1

# Load Data --------------------------------------------------------------------
# Nighttime Lights
viirs_band1 <- raster(file.path(project_file_path, "Data", "RawData","viirs_monthly","kenya_viirs_2012_to_2018_avg_rad.tif"), band=1)

# VIIRS Shapefile
load(file.path(intermediate_data_file_path, "viirs_shapefile", "viirs_band1_shp_africover_urban.Rda"))
viirs_urban_df <- viirs_band1_shp_africover_urban

# GADM
ken_adm0 <- getData('GADM', country='KEN', level=0)
ken_adm1 <- getData('GADM', country='KEN', level=1)
nairobi_boundary <- ken_adm1[ken_adm1$NAME_1 %in% "Nairobi",]

# Nairobi
nairobi <- c(-1.283333, 36.816667) %>% t %>% as.data.frame %>% dplyr::rename(lat=V1) %>% dplyr::rename(lon=V2)
#nairobi$id <- 1
#nairobi_sdf <- nairobi
#coordinates(nairobi_sdf) <- ~lon+lat
#crs(nairobi_sdf) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
viirs_urban_df$distance_nairobi <- sqrt((viirs_urban_df$longitude - nairobi$lon)^2 + (viirs_urban_df$latitude - nairobi$lat)^2) * 111.12

# Aggregate VIIRS Tile ---------------------------------------------------------
year <- c(rep(2012, 9), rep(2013, 12), rep(2014, 12), rep(2015, 12), rep(2016, 12), rep(2017, 12), rep(2018, 6))
month <- c(4:12,          1:12,          1:12,          1:12,          1:12,          1:12,          1:6)
id <- 1:75
raster_ids <- cbind(id, year, month) %>% as.data.frame
raster_ids$year_month <- paste0(raster_ids$year,"-",raster_ids$month)
raster_ids$year_month[raster_ids$month <= 9] <- paste0(raster_ids$year[raster_ids$month <= 9],"-0",raster_ids$month[raster_ids$month <= 9])

viirs_subet_ids <- raster_ids$id[(raster_ids$year_month >= VIIRS_YEAR_MONTH_START) & (raster_ids$year_month <= VIIRS_YEAR_MONTH_END)]

# Raster Brick of Select VIIRS Year/Months
viirs_rad_brick <- stack(file.path(project_file_path, "Data", "RawData","viirs_monthly","kenya_viirs_2012_to_2018_avg_rad.tif")) %>% subset(viirs_subet_ids) %>% brick() %>% crop(ken_adm0) %>% mask(ken_adm0)
viirs_cf_brick <- stack(file.path(project_file_path, "Data", "RawData","viirs_monthly","kenya_viirs_2012_to_2018_cf_cvg.tif")) %>% subset(viirs_subet_ids) %>% brick() %>% crop(ken_adm0) %>% mask(ken_adm0)

# VIIRS pixel needs to have at least 1 cloud free day; otherwise, set NA radiance
for(id in viirs_subet_ids){
  viirs_rad_brick[[paste0("kenya_viirs_2012_to_2018_avg_rad.",id)]][][viirs_cf_brick[[paste0("kenya_viirs_2012_to_2018_cf_cvg.",id)]][] <= VIIRS_NUMBER_CLOUDFREE_PIXELS] <- NA
}

# Collapse raster brick
s <- Sys.time()
viirs_rad <- calc(viirs_rad_brick, fun = function(x) {quantile(x, probs = c(.9),na.rm=TRUE)})
e <- Sys.time()
e-s

# Add raw viirs radiance to viirs-urban dataframe
viirs_urban_df$viirs_rad_raw <- viirs_rad[][!is.na(viirs_rad[])]

# Illustrate Trade Offs of NTL Threshold ---------------------------------------

# Summarize VIIRS Rad values across different %s of urban
summarize_viirs_urban <- function(threshold){
  summarize_viirs_urban_per_threshold <- function(per_urban, threshold){
    
    df_out <- per_urban %>% as.data.frame %>% dplyr::rename(prop_urban = ".")
    df_out$ntl_threshold <- threshold
    
    viirs_values_nourban <- viirs_urban_df$viirs_rad_raw[(viirs_urban_df$africover_urban == 0)]
    df_out$number_nourban <- length(viirs_values_nourban)
    df_out$number_classify_nourban_as_urban <- sum(viirs_values_nourban > threshold)
    df_out$prop_nonurban_class_urban <- df_out$number_classify_nourban_as_urban/df_out$number_nourban
    
    viirs_values_urban <- viirs_urban_df$viirs_rad_raw[(viirs_urban_df$africover_urban >= per_urban)]
    df_out$number_urban <- length(viirs_values_urban)
    df_out$number_classify_urban <- sum(viirs_values_urban > threshold)
    df_out$prop_urban_drop <- 1 - (df_out$number_classify_urban/df_out$number_urban)

    return(df_out)
  }
  
  percent_urban_list <- seq(from=0.01,to=.9,by=.01)
  df_out <- lapply(percent_urban_list, summarize_viirs_urban_per_threshold, threshold) %>% bind_rows
  return(df_out)
}

ntl_thresholds <- seq(from=.1, to=1, by=.1)
viirs_rad_urban_percent_df <- lapply(ntl_thresholds, summarize_viirs_urban) %>% bind_rows
viirs_rad_urban_percent_df$ntl_threshold_factor <- factor(viirs_rad_urban_percent_df$ntl_threshold)

ggplot(viirs_rad_urban_percent_df[viirs_rad_urban_percent_df$ntl_threshold %in% c(0.2,0.3,0.5,0.8,1),]) + 
  geom_line(aes(x=prop_urban, y=prop_urban_drop, color=ntl_threshold_factor), size=1) + 
  labs(x="Proportion Cell Urban, According to Africover [Equal or greater]",
       y="Percent Urban Cells Drop",
       color="") +
  geom_vline(xintercept=.1,col="red") + 
  theme_minimal()

ggplot(viirs_rad_urban_percent_df[as.character(viirs_rad_urban_percent_df$prop_urban) == "0.1",]) + 
  geom_line(aes(x=ntl_threshold, y=number_classify_nourban_as_urban), size=1) + 
  labs(x="NTL Threshold",
       y="Number VIIRS Cells with 0% Africover-Urban Classified as Urban") +
  theme_minimal()
  
ggplot(viirs_rad_urban_percent_df[as.character(viirs_rad_urban_percent_df$ntl_threshold) == "1",]) + 
  geom_line(aes(x=prop_urban, y=number_urban), size=1) + 
  labs(x="Proportion Cell Urban, According to Africover [Equal or greater]",
       y="Number Cells") +
  theme_minimal()

# Set VIIRS Noise to 0 ---------------------------------------------------------
# Determine maximum NTL near Nairobi
viirs_nairobi <- viirs_rad %>% crop(nairobi_boundary) %>% mask(nairobi_boundary)
viirs_dnb_max_nairobi <- viirs_nairobi[] %>% max(na.rm=T)

# Cut offs
viirs_rad[][viirs_rad[] <= VIIRS_DNB_PIXEL_MIN] <- 0
viirs_rad[][viirs_rad[] >= viirs_dnb_max_nairobi] <- 0

# Potential Urban Clusters -----------------------------------------------------
calc_raster_gradient_magnitude <- function(r){
  g <- rast.grad(r)
  r_gm <- sqrt(g$rast.grad.x^2 + g$rast.grad.y^2)
  return(r_gm)
}

# Segmentation
viirs_rad_cluster <- viirs_rad %>% 
  calc_raster_gradient_magnitude %>%
  as.cimg %>% 
  distmap %>% 
  watershed(tolerance=0.3,ext=5)
#plot(colorLabels(getFrame(viirs_rad_segment, 1)))

# Raster of Segmentation Values
viirs_rad_cluster_r <- viirs_rad # Make raster to replace values of
viirs_rad_cluster_r[] <- viirs_rad_cluster[] %>% as.list %>% unlist

# Extract Raster Segment IDs to VIIRS DF Grid
viirs_urban_sdf <- viirs_urban_df@data
coordinates(viirs_urban_sdf) <- ~longitude+latitude
crs(viirs_urban_sdf) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

viirs_urban_sdf$cluster_id <- extract(viirs_rad_cluster_r, viirs_urban_sdf)

# Create Urban-Cluster Level Datasert ------------------------------------------
viirs_urban_df <- viirs_urban_sdf@data
viirs_urban_df$Ncells <- 1
viirs_urban_df$seg_size <- 750

create_cluster_level_dataset <- function(cluster_id, dt){
  dt_i <- dt[dt$cluster_id == cluster_id,]
  proportion_urban_truth <- mean(dt_i$africover_urban)
  
  potential_thresholds <- seq(0,2,by=.05)
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

data_mean_sum <- pbmclapply(sort(unique(viirs_level_dt$cluster_id)), create_cluster_level_dataset, viirs_urban_df, mc.cores=1) %>% bind_rows

# Remove Obviously Non-Urban Clusters ------------------------------------------
data_mean_sum <- data_mean_sum[data_mean_sum$cluster_id != 0,]
data_mean_sum <- data_mean_sum[data_mean_sum$proportion_urban_truth > 0,]

# Method 1 ---------------------------------------------------------------------
data_mean_sum$y <- log((data_mean_sum$viirs_rad_max - data_mean_sum$viirs_rad_min)/(data_mean_sum$ntl_threshold - data_mean_sum$viirs_rad_min)-1)
data_mean_sum$ln_viirs_rad_mean <- -log(data_mean_sum$viirs_rad_mean)
data_mean_sum$ln_Ncells <- -log(data_mean_sum$Ncells)
data_mean_sum$neg_one <- -1

lm1 <- lm(y~ln_viirs_rad_mean+ln_Ncells+neg_one-1, data=data_mean_sum)

n <- coefficients(lm1)[1] %>% as.numeric
a <- coefficients(lm1)[1] %>% as.numeric
b <- coefficients(lm1)[1] %>% as.numeric

data_mean_sum$opt_thresh_predict <- 1/(1+exp(1)^(-(a*data_mean_sum$ln_viirs_rad_mean + b*data_mean_sum$ln_Ncells + n))) * (data_mean_sum$viirs_rad_max-data_mean_sum$viirs_rad_min) + data_mean_sum$viirs_rad_min


plot(data_mean_sum$opt_thresh_predict, data_mean_sum$ntl_threshold)

data_mean_sum$viirs_rad_max
data_mean_sum$viirs_rad_min

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











































# PURGATORY --------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Collapse raster brick
viirs_rad_mean <- calc(viirs_rad_brick, fun = mean, na.rm = T)
viirs_rad_median <- calc(viirs_rad_brick, fun = median, na.rm = T)

viirs_band1_shp_africover_urban$viirs_rad_mean <- viirs_rad_mean[][!is.na(viirs_rad_mean[])]
viirs_band1_shp_africover_urban$viirs_rad_median <- viirs_rad_median[][!is.na(viirs_rad_median[])]

# Cut off to say things aren't urban
viirs_band1_shp_africover_urban$viirs_rad_median[viirs_band1_shp_africover_urban$africover_urban > .1] %>% summary
quantile(viirs_band1_shp_africover_urban$viirs_rad_median[viirs_band1_shp_africover_urban$africover_urban > .1], .1)

# Check cut-offs ---------------------------------------------------------------
# If NTL is greater than THRESH (eg 0.5); number non-urban cells.
# Classify non-urban as 0% urban

ntl_thresh <- 0.5
africover_thresh <- 0
sum(viirs_band1_shp_africover_urban$africover_urban > africover_thresh)
num_incorrectly_class_urban <- sum((viirs_band1_shp_africover_urban$viirs_rad_median > ntl_thresh) & (viirs_band1_shp_africover_urban$africover_urban == africover_thresh))
num_correctly_class_urban <- sum((viirs_band1_shp_africover_urban$viirs_rad_median > ntl_thresh) & (viirs_band1_shp_africover_urban$africover_urban >= africover_thresh))

# Algorithm: Cutoffs -----------------------------------------------------------
africover_urban_percent_cutoff <- 0.1
viirs_thresh <- 0.5
sum_type <- "median"
df <- viirs_band1_shp_africover_urban
subset_info <- ""

det_accuracry_alg_cutoff <- function(viirs_thresh, africover_urban_percent_cutoff, sum_type, subset_info, df){
  df <- df@data
  df$is_urban <- as.numeric(df$africover_urban > africover_urban_percent_cutoff)
  
  num_cells <- nrow(df)
  num_urban_cells <- sum(df$is_urban == 1)
  num_nonurban_cells <- sum(df$is_urban == 0)
  num_correctly_classify_urban <- sum((df$is_urban == 1) & (df[[paste0("viirs_rad_",sum_type)]] >= viirs_thresh))     
  num_incorrectly_classify_urban <- sum((df$is_urban == 0) & (df[[paste0("viirs_rad_",sum_type)]] >= viirs_thresh))     
  num_correctly_classify_nonurban <- sum((df$is_urban == 0) & (df[[paste0("viirs_rad_",sum_type)]] < viirs_thresh))     
  of_urban_prop_correct <- num_correctly_classify_urban/num_urban_cells
  of_nonurban_prop_correct <- num_correctly_classify_nonurban/num_nonurban_cells
  
  urban_cells <- df[df$is_urban == 1,]
  incorrectly_classified_as_urban_cells <- df[(df$is_urban == 0) & (df[[paste0("viirs_rad_",sum_type)]] >= viirs_thresh),]
  
  coordinates(urban_cells) <- ~longitude+latitude
  urban_cells <- urban_cells %>% gBuffer(width=.01/111.12, byid=T) # 10 meters
  urban_cells <- gUnaryUnion(urban_cells)
  
  coordinates(incorrectly_classified_as_urban_cells) <- ~longitude+latitude
  
  #distances <- as.numeric(gDistance(incorrectly_classified_as_urban_cells,urban_cells,byid=T)) * 111.12
  #distance_max <- max(distances)
  #distance_95p <- quantile(distances,.95)
  #distance_90p <- quantile(distances,.90)
  #distance_80p <- quantile(distances,.80)
  #distance_50p <- quantile(distances,.50)
  #distance_25p <- quantile(distances,.25)
  
  #df_out <- cbind(num_cells, num_urban_cells, num_nonurban_cells, num_correctly_classify_urban,
  #                num_incorrectly_classify_urban, num_correctly_classify_nonurban, of_urban_prop_correct,
  #                of_nonurban_prop_correct, africover_urban_percent_cutoff, viirs_thresh, sum_type, subset_info,
  #                distance_max,distance_95p,distance_90p,distance_80p,distance_50p,distance_25p) %>% as.data.frame
  
  df_out <- cbind(num_cells, num_urban_cells, num_nonurban_cells, num_correctly_classify_urban,
                  num_incorrectly_classify_urban, num_correctly_classify_nonurban, of_urban_prop_correct,
                  of_nonurban_prop_correct, africover_urban_percent_cutoff, viirs_thresh, sum_type, subset_info) %>% as.data.frame
  
  return(df_out)
}

results_cutoff_df <- as.data.frame(matrix(nrow=0,ncol=0))
for(ntl_cut_off in c(0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3)){
  results_cutoff_df_i <- det_accuracry_alg_cutoff(ntl_cut_off, .1, "median", "none", viirs_band1_shp_africover_urban[viirs_band1_shp_africover_urban$distance_nairobi > 30,])
  results_cutoff_df <- bind_rows(results_cutoff_df, results_cutoff_df_i)
  
  print(ntl_cut_off)
}

# Algorithm --------------------------------------------------------------------
# Gradient Magitude
calc_raster_gradient_mag <- function(r){
  g <- rast.grad(r)
  r_gm <- sqrt(g$rast.grad.x^2 + g$rast.grad.y^2)
  return(r_gm)
}

# Cut off: Areas below here cannot be urban.
viirs_rad_mean[][viirs_rad_mean[] <= 0.2] <- 0
viirs_rad_median[][viirs_rad_median[] <= 0.2] <- 0

# Plot Areas Above Threshold
if(F){
viirs_rad_mean_coord <- coordinates(viirs_rad_mean) %>% as.data.frame
viirs_rad_mean_coord$pos_ntl <- (viirs_rad_mean[] > 0) %>% as.numeric
viirs_rad_mean_coord <- viirs_rad_mean_coord[!is.na(viirs_rad_mean_coord$pos_ntl),]
coordinates(viirs_rad_mean_coord) <- ~x+y
plot(viirs_rad_mean_coord[viirs_rad_mean_coord$pos_ntl==1,])
}

# Calculate Gradient Magniture
viirs_rad_mean_grad_mad <- calc_raster_gradient_mag(viirs_rad_mean)
viirs_rad_median_grad_mad <- calc_raster_gradient_mag(viirs_rad_median)

# Segmentation
viirs_rad_mean_grad_mad_seg <- viirs_rad_mean_grad_mad %>% as.cimg %>% distmap %>% watershed(tolerance=0.3,ext=5)
viirs_rad_median_grad_mad_seg <- viirs_rad_median_grad_mad %>% as.cimg %>% distmap %>% watershed(tolerance=0.3,ext=5)
# plot(erode(a, kern))

# Plot
png(file.path(figures_file_path, "Urbanization","urban_potential_segments_avgviirs.png"),height=600,width=600)
plot(colorLabels(getFrame(viirs_rad_mean_grad_mad_seg, 1)))
dev.off()

png(file.path(figures_file_path,  "Urbanization","urban_potential_segments_medianviirs.png"),height=600,width=600)
plot(colorLabels(getFrame(viirs_rad_median_grad_mad_seg, 1)))
dev.off()

# To raster
viirs_rad_mean_segval <- viirs_rad_mean
viirs_rad_mean_segval[] <- viirs_rad_mean_grad_mad_seg[] %>% as.list %>% unlist

viirs_rad_median_segval <- viirs_rad_median
viirs_rad_median_segval[] <- viirs_rad_median_grad_mad_seg[] %>% as.list %>% unlist

# Extract Raster Segment Values to Grid
viirs_band1_shp_africover_urban_points <- viirs_band1_shp_africover_urban@data
coordinates(viirs_band1_shp_africover_urban_points) <- ~longitude+latitude
crs(viirs_band1_shp_africover_urban_points) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

viirs_band1_shp_africover_urban_points$segment_value_mean <- extract(viirs_rad_mean_segval, viirs_band1_shp_africover_urban_points)
viirs_band1_shp_africover_urban_points$segment_value_median <- extract(viirs_rad_median_segval, viirs_band1_shp_africover_urban_points)

# Determine Optimal Threshold --------------------------------------------------
data.dt <- as.data.table(viirs_band1_shp_africover_urban_points)
data.dt$Ncells <- 1
data.dt$seg_size <- 750

dt <- data.dt

determine_optimal_threshold_per_cluster <- function(cluster_id, dt, viirs_type){
  dt_i <- dt[dt$segment_value_mean %in% cluster_id,]
  proportion_urban_truth <- mean(dt_i$africover_urban)

  potential_thresholds <- seq(0,3,by=.01)
  if(viirs_type == "mean") proportion_urban_thresh <- lapply(potential_thresholds, function(thresh) mean(dt_i$viirs_rad_mean < thresh)) %>% unlist
  if(viirs_type == "median") proportion_urban_thresh <- lapply(potential_thresholds, function(thresh) mean(dt_i$viirs_rad_median < thresh)) %>% unlist
  
  threshold_urban_df <- cbind(potential_thresholds, proportion_urban_thresh)
  optm_thresh_df <- threshold_urban_df[which.min(abs(threshold_urban_df[,2] - proportion_urban_truth)),] %>% t %>% as.data.frame
  
  # Add other variables
  optm_thresh_df$cluster_id <- cluster_id
  optm_thresh_df$Ncells <- sum(dt_i$Ncells)
  optm_thresh_df$seg_size <- sum(dt_i$seg_size)
  optm_thresh_df$proportion_urban_truth <- proportion_urban_truth
  if(viirs_type == "mean") optm_thresh_df$viirs_rad_mean <- mean(dt_i$viirs_rad_mean)
  if(viirs_type == "mean") optm_thresh_df$viirs_rad_min <- min(dt_i$viirs_rad_mean)
  if(viirs_type == "mean") optm_thresh_df$viirs_rad_max <- max(dt_i$viirs_rad_mean)
  
  if(viirs_type == "median") optm_thresh_df$viirs_rad_mean <- mean(dt_i$viirs_rad_median)
  if(viirs_type == "median") optm_thresh_df$viirs_rad_min <- min(dt_i$viirs_rad_median)
  if(viirs_type == "median") optm_thresh_df$viirs_rad_max <- max(dt_i$viirs_rad_median)
  
  return(optm_thresh_df)
}

data_mean_sum <- lapply(sort(unique(data.dt$segment_value_mean)), determine_optimal_threshold_per_cluster, data.dt, "mean") %>% bind_rows
data_median_sum <- lapply(sort(unique(data.dt$segment_value_median)), determine_optimal_threshold_per_cluster, data.dt, "median") %>% bind_rows

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
