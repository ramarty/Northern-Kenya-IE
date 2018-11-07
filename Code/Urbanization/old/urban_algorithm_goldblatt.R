# Urban 
# Northern Kenya IE

# Setup ------------------------------------------------------------------------
if(Sys.info()["user"] == "WB521633") source("C:/Users/wb521633/Dropbox/World Bank/IEs/Northern Kenya IE/Code/northern_kenya_ie_master.R")

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
viirs_rad_mean <- calc(viirs_rad_brick, fun = mean, na.rm = T)
viirs_rad_median <- calc(viirs_rad_brick, fun = median, na.rm = T)

