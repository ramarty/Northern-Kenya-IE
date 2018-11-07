# Northern Kenya IE Master

# Filepaths --------------------------------------------------------------------
if(Sys.info()["user"] == "WB521633") project_file_path <- "C:/Users/wb521633/Dropbox/World Bank/IEs/Northern Kenya IE"
if(Sys.info()["user"] == "robmarty") project_file_path <- "~/Dropbox/World Bank/IEs/Northern Kenya IE"
if(Sys.info()["user"] == "r521633") project_file_path <- "/home/wb521633/IEs/Northern-Kenya-IE"

raw_data_file_path <- file.path(project_file_path, "Data", "RawData")
intermediate_data_file_path <- file.path(project_file_path, "Data", "IntermediateData")
final_data_file_path <- file.path(project_file_path, "Data", "FinalData")
figures_file_path <- file.path(project_file_path, "Figures")

# Packages ---------------------------------------------------------------------
library(rgdal)
library(raster)
library(ggplot2)
library(sp)
library(rgeos)
library(dplyr)
library(ForestTools)
library(sf)
library(velox)
library(spex)
library(ctmcmove)
library(imager)
library(lidR)
library(rsMove)
library(data.table)
library(parallel)
library(pbmcapply)

#source("https://bioconductor.org/biocLite.R")
#biocLite("EBImage")
library(EBImage)
