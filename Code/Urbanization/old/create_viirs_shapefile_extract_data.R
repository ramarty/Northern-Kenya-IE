# Create Urban Rasters from VIIRS
# Northern Kenya IE

# Setup ------------------------------------------------------------------------
if(Sys.info()["user"] == "WB521633") source("C:/Users/wb521633/Dropbox/World Bank/IEs/Northern Kenya IE/Code/northern_kenya_ie_master.R")

import_raw_africover <- FALSE
create_viirs_shapefile_raw <- FALSE

# Load Data --------------------------------------------------------------------
# Nighttime Lights
viirs_band1 <- raster(file.path(project_file_path, "Data", "RawData","viirs_monthly","kenya_viirs_2012_to_2018_avg_rad.tif"), band=1)

# Roads Data
setwd(file.path(project_file_path, "Data", "IntermediateData", "Roads"))
road_a1 <- readOGR(dsn=".", layer="a1")
road_a2 <- readOGR(dsn=".", layer="a2")
road_b9 <- readOGR(dsn=".", layer="b9")
road_a109 <- readOGR(dsn=".", layer="a109")

buff_width <- 50
road_a1_buff <- gBuffer(road_a1, width=buff_width/111.12, byid=F)
road_a2_buff <- gBuffer(road_a2, width=buff_width/111.12, byid=F)
road_b9_buff <- gBuffer(road_b9, width=buff_width/111.12, byid=F)
road_a109_buff <- gBuffer(road_a109, width=buff_width/111.12, byid=F)

# GADM
ken_adm0 <- getData('GADM', country='KEN', level=0)

# Nairobi
nairobi <- c(-1.283333, 36.816667) %>% t %>% as.data.frame %>% dplyr::rename(lat=V1) %>% dplyr::rename(lon=V2)

# Africover
if(import_raw_africover){
  convert_africover_urban <- function(r){
    r[] <- as.numeric(r[] %in% c(8))
    return(r)
  }
  
  africover <- raster(file.path(file_path, "Data", "RawData", "Africover", "ESACCI-LC-L4-LC10-Map-20m-P1Y-2016-v1.0.tif"), band=1)
  africover <- crop(africover, extent(viirs_band1))
  africover_urban <- calc(africover, fun=convert_africover_urban)
  
  writeRaster(africover, file.path(file_path, "Data", "IntermediateData", "Africover", "africover_kenya.tif"))
  writeRaster(africover_urban, file.path(file_path, "Data", "IntermediateData", "Africover", "africover_kenya_urban.tif"))
} else{
  africover_urban <- raster(file.path(project_file_path, "Data", "IntermediateData", "Africover", "africover_kenya_urban.tif"))
}

# VIIRS-Level Shapefile --------------------------------------------------------
if(create_viirs_shapefile_raw){
  viirs_band1 <- viirs_band1 %>% crop(ken_adm0) %>% mask(ken_adm0)
  viirs_band1_shp <- rasterToPolygons(viirs_band1, na.rm=T)
  save(viirs_band1_shp, file=file.path(intermediate_data_file_path, "viirs_shapefile", "viirs_band1_shp.Rda"))
} else{
  load(file.path(intermediate_data_file_path, "viirs_shapefile", "viirs_band1_shp.Rda"))
}

# Split Africover into Multiple Rasters ----------------------------------------
viirs_band1_shp$africover_urban <- NA

chunk_size <- 50000
start_ids <- seq(1, nrow(viirs_band1_shp), by=chunk_size)

for(start_i in start_ids){
  end_i <- min(start_i + chunk_size - 1, nrow(viirs_band1_shp))

  num_range <- start_i:end_i
  viirs_band1_shp_i <- viirs_band1_shp[num_range,]
  africover_urban_i <- africover_urban %>% crop(viirs_band1_shp_i)
  
  viirs_band1_shp$africover_urban[num_range] <- as.numeric(velox(africover_urban_i)$extract(sp=viirs_band1_shp_i, fun=function(x){mean(x, na.rm=TRUE)}))

  print(start_i)
}

# Calculate Distances ----------------------------------------------------------
viirs_band1_shp$longitude <- coordinates(viirs_band1_shp)[,1] %>% as.numeric
viirs_band1_shp$latitude <- coordinates(viirs_band1_shp)[,2] %>% as.numeric

viirs_band1_shp$distance_nairobi <- sqrt((viirs_band1_shp$longitude - nairobi$lon)^2 + (viirs_band1_shp$latitude - nairobi$lat)^2) * 111.12

# Distance to roads cutoffs
viirs_band1_shp_points <- viirs_band1_shp@data
coordinates(viirs_band1_shp_points) <- ~longitude+latitude
crs(viirs_band1_shp_points) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

for(buff_width in c(5,10,20,50)){
  road_a1_buff <- gBuffer(road_a1, width=buff_width/111.12, byid=F)
  road_a2_buff <- gBuffer(road_a2, width=buff_width/111.12, byid=F)
  road_b9_buff <- gBuffer(road_b9, width=buff_width/111.12, byid=F)
  road_a109_buff <- gBuffer(road_a109, width=buff_width/111.12, byid=F)
  
  road_a1_buff$in_buff <- 1
  road_a2_buff$in_buff <- 1
  road_b9_buff$in_buff <- 1
  road_a109_buff$in_buff <- 1
  
  viirs_band1_shp[[paste0("within_",buff_width,"km_a1")]] <- over(viirs_band1_shp_points, road_a1_buff)$in_buff
  viirs_band1_shp[[paste0("within_",buff_width,"km_a2")]] <- over(viirs_band1_shp_points, road_a2_buff)$in_buff
  viirs_band1_shp[[paste0("within_",buff_width,"km_b9")]] <- over(viirs_band1_shp_points, road_b9_buff)$in_buff
  viirs_band1_shp[[paste0("within_",buff_width,"km_a109")]] <- over(viirs_band1_shp_points, road_a109_buff)$in_buff
  
  print(buff_width)
}

# Export -----------------------------------------------------------------------
viirs_band1_shp_africover_urban <- viirs_band1_shp
save(viirs_band1_shp_africover_urban, file=file.path(intermediate_data_file_path, "viirs_shapefile", "viirs_band1_shp_africover_urban.Rda"))












# Extract land cover data to VIIRS resolution ----------------------------------
#viirs_band1_shp <- rasterToPolygons(viirs_band1)

viirs_band1_coords <- coordinates(viirs_band1) %>% as.data.frame
names(viirs_band1_coords) <- c("lat", "lon")
viirs_band1_coords$id <- 1:nrow(viirs_band1_coords)
coordinates(viirs_band1_coords) <- ~lat+lon

viirs_over_urban <- extract(africover_urban, viirs_band1_coords)

# NTL to Urban -----------------------------------------------------------------
viirs_band1[] <- log(viirs_band1[] + 1)

viirs_band1_bin <- viirs_band1
viirs_band1_bin[] <- as.numeric(viirs_band1[] >= .5)

table(viirs_over_urban, viirs_band1_bin[])
table(viirs_over_urban, crowns_r_bin[])

viirs_band1 <- raster(file.path(file_path, "Data", "RawData","viirs_monthly","kenya_viirs_2012_to_2018_avg_rad.tif"), band=1)
viirs_band1[] <- log(viirs_band1[] + 1)

lin <- function(x){x * 0.0005 + 0.6}
ttops <- vwf(CHM = viirs_band1, winFun = lin, minHeight = 0.5, maxWinDiameter=NULL)
crowns_r <- mcws(treetops = ttops, CHM = viirs_band1, minHeight = .5, verbose = FALSE)
crowns_r_bin <- crowns_r
crowns_r_bin[] <- as.numeric(crowns_r[] > 0)

crowns_p <- mcws(treetops = ttops, CHM = viirs_band1, minHeight = .5, verbose = FALSE, format="polygons")


# Look at certain cities
wajir <- c(1.749975, 40.055116) %>% t %>% as.data.frame
dela <- c(2.3030871879724506, 39.70882137112494) %>% t %>% as.data.frame

coordinates(wajir) <- ~V2+V1
coordinates(dela) <- ~V2+V1

crs(wajir) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
crs(dela) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

extent(gBuffer(wajir, 0.1, byid = T))
wajir <- extent(gBuffer(wajir, .1, byid=T)) - .01

a <- crop(africover_urban, extent(gBuffer(dela, 0.00001, byid = T)))




globcover.2012.urban


# Remove plot margins (optional)
par(mar = rep(0.5, 4))

# Plot CHM (extra optional arguments remove labels and tick marks from the plot)
plot(kootenayCHM, xlab = "", ylab = "", xaxt='n', yaxt = 'n')
lin <- function(x){x * 0.05 + 0.6}
ttops <- vwf(CHM = kootenayCHM, winFun = lin, minHeight = 2)

# Create crown map
crowns <- mcws(treetops = ttops, CHM = kootenayCHM, minHeight = 1.5, verbose = FALSE)

# Plot crowns
plot(crowns, col = sample(rainbow(50), length(unique(crowns[])), replace = TRUE), legend = FALSE, xlab = "", ylab = "", xaxt='n', yaxt = 'n')


ntl <- raster(file.path(file_path, "Data", "RawData","viirs_monthly","kenya_viirs_2012_to_2018_avg_rad.tif"), band=1)
ntl[] <- log(ntl[]+1)
ntl_bin <- ntl
ntl_bin[] <- as.numeric(ntl_bin[] > 0.5)
lin <- function(x){x * 0.00002 + 0.01}
ttops <- vwf(CHM = ntl, winFun = lin, minHeight = 2)
