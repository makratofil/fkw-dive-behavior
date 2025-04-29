## geoprocess_for_seafloor_depth.R: for analysis of dive depths relative to 
## seafloor depths, geoprocess multiple imputation pseudotracks 

## 03b: multiple imputation

## for this analysis: only doing seafloor depth and time variables

## Author: Michaela A. Kratofil, Oregon State University, Cascadia Research
## Updated: 23 Feb 2025

## -------------------------------------------------------------------------- ##

## load packages 
library(lubridate)
library(oce)
library(lunar)
library(suntools)
library(sf)
library(stars)
library(purrr)
library(dplyr)
library(here)

## read in all behavior log pseudotracks and prep for geoprocessing ## ------ ##
files <- list.files(path = here("pipeline","crawl_pseudotracks_multiple_imputation"),
                    pattern = "2025Feb23", full.names = T, recursive = F)
files

# function to read in and format the files 
fread <- function(x){
  
  d <- read.csv(x, header = T, stringsAsFactors = F)
  d$start_utc <- as.POSIXct(d$start_utc, tz = "UTC")
  d$start_hst <- as.POSIXct(d$start_hst, tz = "Pacific/Honolulu")
  
  return(d)
}

# apply the function
dfs <- lapply(files, fread)
tags <- bind_rows(dfs)

# filter "p" locations 
tags2 <- filter(tags, locType == "p")

## raster-based variables ## ------------------------------------------------- ##

# source the geoprocessing helper functions
source(here("code","data_processing","oce_sun_moon_angle_devel_functions.R"))

# read in raster files (three different resolutions). read in from outside 
# directory so don't have to unnecessarily copy large files into this 
# directory

rasts <- "D:/04_OSU/06_Thesis/05_Analyses/pc_tag_processing/data/rasters/"

depthL <- stars::read_stars(paste0(rasts, "GebcodepthWUTM4N.nc"))

depthM <- stars::read_stars(paste0(rasts, "FalkorDepthUTM4N.nc"))

depthH <- stars::read_stars(paste0(rasts, "multibeamUTM4N.nc"))

aspectL <- stars::read_stars(paste0(rasts, "GebcoAspectWUTM4N.nc"))

aspectM <- stars::read_stars(paste0(rasts, "FalkorAspectUTM4N.nc"))

aspectH <- stars::read_stars(paste0(rasts, "MultibeamAspectUTM4N.nc"))

slopeL <- stars::read_stars(paste0(rasts, "GebcoSlopeWUTM4N.nc"))

slopeM <- stars::read_stars(paste0(rasts, "FalkorSlopeUTM4N.nc"))

slopeH <- stars::read_stars(paste0(rasts, "MultibeamSlopeUTM4N.nc"))


# make tags dataframe into a spatial object 
tags_sf <- st_as_sf(tags2, coords = c("lon","lat"), crs = 4326) %>%
  st_transform(crs = st_crs(depthH))

tags_sf <- tags_sf %>%
  mutate(
    # bathymetry variables
    depthH = as.numeric(stars::st_extract(depthH, ., method = "bilinear")[[1]]),
    slopeH = as.numeric(stars::st_extract(slopeH, ., method = "bilinear")[[1]]),
    aspectH = as.numeric(stars::st_extract(aspectH, ., method = "bilinear")[[1]]),
    depthM = as.numeric(stars::st_extract(depthM, ., method = "bilinear")[[1]]),
    slopeM = as.numeric(stars::st_extract(slopeM, ., method = "bilinear")[[1]]),
    aspectM = as.numeric(stars::st_extract(aspectM, ., method = "bilinear")[[1]]),
    depthL = as.numeric(stars::st_extract(depthL, ., method = "bilinear")[[1]]),
    slopeL = as.numeric(stars::st_extract(slopeL, ., method = "bilinear")[[1]]),
    aspectL = as.numeric(stars::st_extract(aspectL, ., method = "bilinear")[[1]])
  ) %>%
  # for temporal variables, get back into geographic coords
  st_transform(crs = 4326) %>%
  mutate(
    sunrise = as.POSIXct(suntools::sunriset(., start_hst, direction = "sunrise", POSIXct.out=T)$time),
    sunset = as.POSIXct(suntools::sunriset(., start_hst, direction = "sunset", POSIXct.out=T)$time),
    solar_noon = as.POSIXct(suntools::solarnoon(., start_hst, POSIXct.out=T)$time),
    civil_dawn = as.POSIXct(suntools::crepuscule(., start_hst, solarDep = 6, direction="dawn", POSIXct.out=T)$time),
    end_dawn = as.POSIXct(suntools::crepuscule(., start_hst, solarDep = -6, direction="dawn", POSIXct.out=T)$time),
    civil_dusk = as.POSIXct(suntools::crepuscule(., start_hst, solarDep = 6, direction="dusk", POSIXct.out=T)$time),
    start_dusk = as.POSIXct(suntools::crepuscule(., start_hst, solarDep = -6, direction="dusk", POSIXct.out=T)$time)
  )

summary(tags_sf)
str(tags_sf)

# get data back into regular dataframe with lat/lon, then add sun and moon
# angle and azimuth information (non-sf operation)
coords <- as.data.frame(st_coordinates(tags_sf))
pts_df <- bind_cols(tags_sf, coords) %>%
  rename(
    lon = X,
    lat = Y
  ) %>%
  st_drop_geometry() %>%
  mutate(
    ref = rep(T, length(start_utc)),
    sun_azimuth = sunAngle(start_utc, lon, lat, ref)$azimuth,
    sun_altitude = sunAngle(start_utc, lon, lat, ref)$altitude,
    moon_azimuth = moonAngle(start_utc, lon, lat, ref)$azimuth,
    moon_altitude = moonAngle(start_utc, lon, lat, ref)$altitude,
    moon_ill_fraction = moonAngle(start_utc, lon, lat, ref)$illuminatedFraction,
    moon_phase = lunar.phase(start_utc, name = T),
    tod = case_when(
      start_hst < civil_dawn | start_hst > civil_dusk ~ "night",
      start_hst > end_dawn & start_hst < start_dusk ~ "day",
      start_hst >= civil_dawn & start_hst <= end_dawn ~ "dawn",
      start_hst >= start_dusk & start_hst <= civil_dusk ~ "dusk"
    )
  )

## save the file ## ---------------------------------------------------------- ##
write.csv(pts_df, here("pipeline","geoprocessed","all_behavlog_pseudotracks_20imputes_rerouted20mIso_geoprocessed_seafloor_2025Feb23.csv"), row.names = F)
saveRDS(pts_df, here("pipeline","geoprocessed","all_behavlog_pseudotracks_20imputes_rerouted20mIso_geoprocessed_seafloor_2025Feb23.rds"))
