# Create Hexagons Landsat 
# Northern Kenya IE

# Setup ------------------------------------------------------------------------
if(Sys.info()["user"] == "WB521633") source("C:/Users/wb521633/Dropbox/World Bank/IEs/Northern Kenya IE/Code/northern_kenya_ie_master.R")

library(raster)
library(rgdal)

buffer_size <- 20 #km

# Load Data --------------------------------------------------------------------
#### Study Area
setwd(file.path(intermediate_data_file_path, "Treatment Roads Buffer"))
roads_study_area <- readOGR(dsn=".", layer=paste0("kenya_treat_roads_",buffer_size,"km_buff"))

#### Kenya ADM
setwd(file.path(raw_data_file_path,"GADM"))
ken_adm0 <- getData('GADM', country='KEN', level=0)

#### Africover
africover <- raster(file.path(intermediate_data_file_path, "Africover", "africover_kenya.tif"))

#### Landsat
africover <- raster(file.path(raw_data_file_path, "Landsat", "kenya_road_1_20km_landsat_2016_median-0000000000-0000000000.tif"))



