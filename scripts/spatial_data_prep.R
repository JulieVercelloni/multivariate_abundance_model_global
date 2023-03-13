### Preparing human densities spatial layer for modelling 
#Reference dataset:
#Center for International Earth Science Information Network - CIESIN - Columbia University. 2018. 
#Gridded Population of the World, Version 4 (GPWv4): Population Density Adjusted to Match 2015 Revision UN WPP Country Totals, Revision 11. Palisades, 
#New York: NASA Socioeconomic Data and Applications Center (SEDAC). https://doi.org/10.7927/H4F47M65. Accessed 06 December 2019.

# Previous manipulations aimed to re-sample, re-project and crop the raw raster to match with the other global disturbances used in the study.  

rm(list=ls())

source("R/packages.R")
source("R/functions.R")

# Read spatial layer 
density_84 <- raster("data/density_resampled_proj.tif")

#### Read the reef data
# ######################################################### Transform use 
proy_merc<-CRS("+proj=merc +lon_0=150 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")

# # Make it spatial
Rec_spat_sf <- read.csv('data/reef_data.csv')%>%
  st_as_sf(coords = c('lng', 'lat'), crs = st_crs(4326)) %>%
  st_transform(crs = proy_merc) 

Rec_spat_sp <- as_Spatial(Rec_spat_sf)

# Extract values of human densities that were the closest from the observations 

new.spdf1 <- extract.mod(density_84, Rec_spat_sp)

## Check number of NAs

sum(is.na(new.spdf1$extract.col)) # need to be = 0

# Transform back to long and lat 

reef_data_with_density_sf <- Rec_spat_sf%>%
  dplyr::select(reef_name) %>%
  cbind(new.spdf1@data$extract.col) %>%
  dplyr::rename(Density = new.spdf1.data.extract.col) %>% distinct() %>%
  st_transform(crs = 4326) 
