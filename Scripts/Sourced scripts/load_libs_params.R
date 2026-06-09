# Load packages
library(tidyverse)
library(odbc)
library(DBI)
library(getPass)
library(keyring)
library(lifecycle)
library(data.table)
library(crabpack)
library(sf)
library(gstat)
library(rnaturalearth)
library(raster)
library(concaveman)
library(png)
library(mgcv)
library(dplyr)
library(mgcv)
library(ggplot2)
library("rnaturalearth")
library(patchwork)
library(gratia)
library(MuMIn)
library(DHARMa)
library(mgcViz)
library(akgfmaps)
library(INLA)
library(sdmTMB)
library(glmmTMB)
library(broom)
library(gstat)
library(devtools)
#install_github("vast-lib/tinyVAST", dependencies = TRUE)
library(tinyVAST)
library(ecmwfr)
library(tidync)

# # Set years
# current.year <- 2025
# years <- c(1989:current.year)

# Specify directories
dir <- "Y:/KOD_Research/Ryznar/Crab functional maturity"

data_dir <- "Y:/KOD_Survey/EBS Shelf/Data_Processing/Data/" # for survey data

remote_dir <- "Y:/KOD_Research/Ryznar/Crab functional maturity/"

# CRS for spatial stufff
region_layers <- akgfmaps::get_base_layers("sebs")

map.crs <- region_layers$crs

in.crs = "+proj=longlat +datum=NAD83"

crs.latlon <- "epsg:4326" 

# Read in spatial layers
source("Y:/KOD_Survey/EBS Shelf/Spatial crab/load.spatialdata.R")


