# Create Hexagons Landsat 
# Northern Kenya IE

# Setup ------------------------------------------------------------------------
if(Sys.info()["user"] == "WB521633") source("C:/Users/wb521633/Dropbox/World Bank/IEs/Northern Kenya IE/Code/northern_kenya_ie_master.R")

library(raster)

buffer_size <- 20 #km

# Load Data --------------------------------------------------------------------
setwd(file.path(intermediate_data_file_path, "Roads"))

a1 <- readOGR(dsn=".", layer="a1")
a2 <- readOGR(dsn=".", layer="a2")
b9 <- readOGR(dsn=".", layer="b9")

a1_buff <- gBuffer(a1, width=buffer_size/111.12) 
a2_buff <- gBuffer(a2, width=buffer_size/111.12) 
b9_buff <- gBuffer(b9, width=buffer_size/111.12) 

a1_buff$id <- 1
a2_buff$id <- 2
b9_buff$id <- 3

a1_buff$road <- "a1"
a2_buff$road <- "a2"
b9_buff$road <- "b9"

roads_buff <- rbind(a1_buff, a2_buff, b9_buff)

setwd(file.path(intermediate_data_file_path, "Treatment Roads Buffer"))
writeOGR(obj=roads_buff, dsn=".", layer=paste0("kenya_treat_roads_",buffer_size,"km_buff"),driver="ESRI Shapefile", overwrite_layer=T)
#writeOGR(obj=roads_buff, dsn="kenya_treat_roads_20km_buff_geoj", layer="kenya_treat_roads_20km_buff_geoj",driver="GeoJSON", overwrite_layer=T)

