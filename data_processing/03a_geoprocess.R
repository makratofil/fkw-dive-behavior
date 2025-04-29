## geoprocess.R: populate behavior log pseudotracks with temporal and spatial
## variables for further analysis.

## 03a: regular pseudotracks (i.e., single track, not multiple imputation)

## Author: Michaela A. Kratofil, Oregon State University, Cascadia Research
## Updated: 18 Feb 2025

## ---------------------------------------------------------------------------- ##

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

## read in all behavior log pseudotracks and prep for geoprocessing ## -------- ##
files <- list.files(path = here("pipeline","crawl_pseudotracks"),
                    pattern = "2025Feb18.csv", full.names = T, recursive = F)
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
str(tags)

## raster-based variables ## ------------------------------------------------- ##

# source the geoprocessing helper functions (specifically, devel version of 
## oce functions)
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
tags_sf <- st_as_sf(tags, coords = c("lon","lat"), crs = 4326) %>%
  st_transform(crs = st_crs(depthH))

tags_sf <- tags_sf %>%
  mutate(
    # bathymetric variables
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

## add chlorophyll-a and mixed layer depth ## -------------------------------- ##

## create list of netCDF files to extract 
ch.files <- list.files(here("data","copernicus_oceancolor")) # chla
ch.files
mld.files <- list.files(here("data","copernicus_global_ocean_physics_reanalysis_mld")) # mld
mld.files

# read in each chla file for each deployment year 
chla2010 <- read_ncdf(here("data","copernicus_oceancolor", ch.files[1]), eps = 1e-3)
chla2012 <- read_ncdf(here("data","copernicus_oceancolor", ch.files[2]), eps = 1e-3)
chla2013 <- read_ncdf(here("data","copernicus_oceancolor", ch.files[3]), eps = 1e-3)
chla2015 <- read_ncdf(here("data","copernicus_oceancolor", ch.files[4]), eps = 1e-3)
chla2017 <- read_ncdf(here("data","copernicus_oceancolor", ch.files[5]), eps = 1e-3)
chla2021 <- read_ncdf(here("data","copernicus_oceancolor", ch.files[6]), eps = 1e-3)
chla2023 <- read_ncdf(here("data","copernicus_oceancolor", ch.files[7]), eps = 1e-3)
chla2024 <- read_ncdf(here("data","copernicus_oceancolor", ch.files[8]), eps = 1e-3)

# read in each mld file for each deployment year
mld2010 <- read_ncdf(here("data","copernicus_global_ocean_physics_reanalysis_mld", mld.files[1]), eps = 1e-3)
mld2012 <- read_ncdf(here("data","copernicus_global_ocean_physics_reanalysis_mld", mld.files[2]), eps = 1e-3)
mld2013 <- read_ncdf(here("data","copernicus_global_ocean_physics_reanalysis_mld", mld.files[3]), eps = 1e-3)
mld2015 <- read_ncdf(here("data","copernicus_global_ocean_physics_reanalysis_mld", mld.files[4]), eps = 1e-3)
mld2017 <- read_ncdf(here("data","copernicus_global_ocean_physics_reanalysis_mld", mld.files[5]), eps = 1e-3)
mld2021 <- read_ncdf(here("data","copernicus_global_ocean_physics_reanalysis_mld", mld.files[6]), eps = 1e-3)
mld2023 <- read_ncdf(here("data","copernicus_global_ocean_physics_reanalysis_mld", mld.files[7]), eps = 1e-3)
mld2024 <- read_ncdf(here("data","copernicus_global_ocean_physics_reanalysis_mld", mld.files[8]), eps = 1e-3)
st_crs(mld2024) == st_crs(chla2024)

# nest the data by deployment year, so we can extract data per file  
tags_sf$dmy <- as.character(date(tags_sf$start_utc))
dmy_nest <- tags_sf %>%
  st_transform(crs = st_crs(chla2010)) %>% # projection the same as the files (same for chla and mld)
  group_by(dmy) %>%
  tidyr::nest()

# create a function for identifying corresponding ocean vars files and extracting values at each
# of the points by date 
get_vars <- function(data){
  
  # for testing
  #data <- dmy_nest[[2]][[2]]
  
  # get the deployment year 
  dep_yr <- as.character(first(year(data$start_utc)))
  
  # get the dmy 
  dmy <- as.character(date(first(data$start_utc)))
  
  # assign the dataset based on the deployment year 
  if(dep_yr == "2010" | dep_yr == "2011"){
    chla.file <- chla2010
    mld.file <- mld2010
  } else if(dep_yr == "2012"){
    chla.file <- chla2012
    mld.file <- mld2012
  } else if(dep_yr == "2013"){
    chla.file <- chla2013
    mld.file <- mld2013
  } else if(dep_yr == "2015"){
    chla.file <- chla2015
    mld.file <- mld2015
  } else if(dep_yr == "2017"){
    chla.file <- chla2017
    mld.file <- mld2017
  } else if(dep_yr == "2021"){
    chla.file <- chla2021
    mld.file <- mld2021
  } else if(dep_yr == "2023"){
    chla.file <- chla2023
    mld.file <- mld2023
  } else if(dep_yr == "2024"){
    chla.file <- chla2024
    mld.file <- mld2024
  }
  
  # filter the chla file for the dmy
  chla_day <- filter(chla.file, stringr::str_detect(time, as.character(dmy)))
  
  # filter the mld file for the dmy
  mld_day <- filter(mld.file, stringr::str_detect(time, as.character(dmy)))
  
  # extract chla and mld values for each point within data
  data_new <- data %>%
    mutate(
      chla = as.numeric(stars::st_extract(chla_day, .)$CHL),
      chla_unc = as.numeric(st_extract(chla_day, .)$CHL_uncertainty),
      chla_flag = as.numeric(st_extract(chla_day, .)$flags),
      mld = as.numeric(stars::st_extract(mld_day, .)$mlotst)
    )
  
  # return the new dataframe
  return(data_new)
  
}

# apply the function
dmy_nest <- dmy_nest %>%
  mutate(
    oc_data = purrr::map(data, ~get_vars(.x))
  )

# unnest the data 
dmy_df <- dmy_nest %>%
  dplyr::select(dmy, oc_data) %>%
  tidyr::unnest() %>%
  sf::st_sf()

# get data back into regular dataframe with lat/lon, then add sun and moon
# angle and azimuth information (non-sf operation)
coords <- as.data.frame(st_coordinates(dmy_df))
pts_df <- bind_cols(dmy_df, coords) %>%
  rename(
    lon = X,
    lat = Y
  ) %>%
  arrange(DeployID, start_utc) %>%
  st_drop_geometry() %>%
  dplyr::select(-dmy) %>%
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

# review the file 
summary(pts_df)

## save the file ## ---------------------------------------------------------- ##

# for geoprocessed to by split by time of day:
write.csv(pts_df, here("pipeline","geoprocessed","all_behavlog_pseudotracks_rerouted20mIso_geoprocessed_2025Feb18.csv"), row.names = F)
saveRDS(pts_df, here("pipeline","geoprocessed","all_behavlog_pseudotracks_rerouted20mIso_geoprocessed_2025Feb18.rds"))

