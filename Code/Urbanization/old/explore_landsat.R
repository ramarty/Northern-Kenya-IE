
# Setup ------------------------------------------------------------------------
if(Sys.info()[["user"]] == "robmarty") project_file_path <- "~/Dropbox/World Bank/IEs/Northern Kenya IE/"
if(Sys.info()[["user"]] == "WB521633") project_file_path <- "C:/Users/wb521633/Dropbox/World Bank/IEs/Northern Kenya IE/"

library(raster)
library(rgdal)
library(dplyr)

a <- raster(file.path(project_file_path, "Data", "RawData", "Landsat", "kenya_road_1_20km_landsat_2016_median-0000000000-0000000000.tif"))
aa <- a[][!is.na(a[])]
aa %>% length
aa %>% head



