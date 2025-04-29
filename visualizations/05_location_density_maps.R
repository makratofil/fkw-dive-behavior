## 05_location_density_maps.R: make maps of location (from satellite tags)
## density to visually compare with the density/distribution of dives

## Author: Michaela A. Kratofil, Oregon State University, Cascadia Research
## Updated: 29 Apr 2025

## --------------------------------------------------------------------------- ##

## load packages 
library(dplyr)
library(lubridate)
library(here)
library(sf)
library(ggplot2)
library(ggspatial)

## read in crawl locations ## ------------------------------------------------ ##
crawl <- readRDS(here("pipeline","PcTags_SPLASH_thru092_Kalman_Crawl1hr_rerouted20mIso_wSegments_2024Dec29.rds"))

# get original time points (i.e., original Argos/GPS locations but smoothed by crawl)
o <- filter(crawl, locType == "o")

## read in behavior log to get dives ## -------------------------------------- ##
beh <- readRDS(here("pipeline","clean_data_for_analysis",
                    "all_behavlog_pseudotracks_rerouted20mIso_geoprocessed_split_tod_2025Feb18.rds"))

# filter out dives and make into sf object for mapping 
beh_sf <- beh %>%
  filter(., !is.na(lon)) %>%
  st_as_sf(., coords = c("lon","lat"), crs = 4326) %>%
  mutate(timestamp = start_utc,
         animal = DeployID)

dives_sf <- beh %>%
  filter(., What == "Dive")%>%
  st_as_sf(., coords = c("lon","lat"), crs = 4326)

# get the last timestamp for each animal
endt <- beh %>%
  group_by(DeployID) %>%
  summarise(
    last = last(end_utc)
  ) 
endt$animal <- endt$DeployID

# remove locations from argos track that occurred after the last behavior log record
o <- left_join(o, endt, by = "animal")

o_sub <- o %>%
  group_by(animal) %>%
  filter(timestamp <= last) %>%
  ungroup()

o_sub %>%
  group_by(animal) %>%
  tally()

# add population column 
o_sub <- o_sub %>%
  mutate(
    population = case_when(
      animal %in% c("PcTag026","PcTag028","PcTag030","PcTag032","PcTag055",
                      "PcTag074") ~ "MHI",
      animal %in% c("PcTag035","PcTag037","PcTag049") ~ "NWHI",
      animal %in% c("PcTag090","PcTag092","PcTagP09") ~ "Open-ocean"
    )
  )

## get map objects ## -------------------------------------------------------- ##

# read in coastline for mapping 
coast <- st_read(here("data","shapefiles","FisheriesIslands.shp"))

## make MHI grid and map ## -------------------------------------------------- ##
mhi <- filter(o_sub, population == "MHI") %>%
  st_transform(., crs = 26904)

# make hexagonal grid with target area of X squared kilometers 
cell_area <- units::as_units(250, "km^2")
grid_spacing <- sqrt(2*cell_area/sqrt(3))

hex_grid <- sf::st_make_grid(
  st_as_sfc(st_bbox(mhi)),
  cellsize = grid_spacing,
  square = F,
  crs = st_crs(mhi)) %>% st_sf() %>%
  tibble::rowid_to_column(var = "hexid")

# now need to join points with hex grid and count the number of points per grid cell
hexbins <- st_join(mhi, hex_grid, join = st_intersects) %>% 
  sf::st_set_geometry(NULL) %>% 
  group_by(hexid) %>% 
  summarise(n = n()) %>% 
  right_join(hex_grid, by = "hexid") %>% 
  sf::st_sf()

# calculate the proportion of locations in each bin
ntag <- nrow(mhi)
hexbins <- hexbins %>%
  mutate(
    prop = n/ntag
  )

# remove bins without any location data
hexbins_sub <- filter(hexbins, !is.na(n))

# map
ggplot() +

  # location density grid
  geom_sf(data = hexbins_sub, color = "gray79", aes(fill = prop)) +
  scale_fill_gradient(low = "white", high = "#015b58") +
  
  # coastline
  geom_sf(data = coast, color = NA, fill = "gray19") +
  
  # aesthetics
  coord_sf(crs = 4326,
           xlim = c(-160.5, -154.7),
           ylim = c(18.7, 22.4)) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black", size = 11),
    legend.text = element_text(color = "black", size = 11)
  ) +
  ggspatial::annotation_scale(text_cex = 0.9, text_col = "black") +
  ggspatial::annotation_north_arrow(location = "tr")

ggsave(here("outputs","mhi_original_location_density_map_250kmhexgrid_2025Mar17_v2.png"),
       width = 9, height = 5, units = "in")

## make NWHI grid and map ## -------------------------------------------------- ##
nwhi <- filter(o_sub, population == "NWHI") %>%
  st_transform(., crs = 26904)

# make hexagonal grid with target area of X squared kilometers 
cell_area <- units::as_units(250, "km^2")
grid_spacing <- sqrt(2*cell_area/sqrt(3))

hex_grid <- sf::st_make_grid(
  st_as_sfc(st_bbox(nwhi)),
  cellsize = grid_spacing,
  square = F,
  crs = st_crs(nwhi)) %>% st_sf() %>%
  tibble::rowid_to_column(var = "hexid")

# now need to join points with hex grid and count the number of points per grid cell
hexbins <- st_join(nwhi, hex_grid, join = st_intersects) %>% 
  sf::st_set_geometry(NULL) %>% 
  group_by(hexid) %>% 
  summarise(n = n()) %>% 
  right_join(hex_grid, by = "hexid") %>% 
  sf::st_sf()

# calculate the proportion of locations in each bin
ntag <- nrow(nwhi)
hexbins <- hexbins %>%
  mutate(
    prop = n/ntag
  )

# remove bins without any location data
hexbins_sub <- filter(hexbins, !is.na(n))

# map
ggplot() +

  # location density grid
  geom_sf(data = hexbins_sub, color = "gray79", aes(fill = prop)) +
  scale_fill_gradient(low = "white", high = "#89689d") +
  
  # coastline
  geom_sf(data = coast, color = NA, fill = "gray19") +
  
  # aesthetics
  coord_sf(crs = 4326,
           xlim = c(-164, -159),
           ylim = c(21, 24.4)) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black", size = 11),
    legend.text = element_text(color = "black", size = 11)
  ) +
  ggspatial::annotation_scale(text_cex = 0.9, text_col = "black") +
  ggspatial::annotation_north_arrow(location = "tr")

ggsave(here("outputs","nwhi_original_location_density_map_250kmhexgrid_2025Mar17_v2.png"),
       width = 7, height = 5, units = "in")

## make pelagic grid and map ## ---------------------------------------------- ##
pel <- filter(o_sub, population == "Open-ocean") %>%
  st_transform(., crs = 26904)

# make hexagonal grid with target area of X squared kilometers 
cell_area <- units::as_units(250, "km^2")
grid_spacing <- sqrt(2*cell_area/sqrt(3))

hex_grid <- sf::st_make_grid(
  st_as_sfc(st_bbox(pel)),
  cellsize = grid_spacing,
  square = F,
  crs = st_crs(pel)) %>% st_sf() %>%
  tibble::rowid_to_column(var = "hexid")

# now need to join points with hex grid and count the number of points per grid cell
hexbins <- st_join(pel, hex_grid, join = st_intersects) %>% 
  sf::st_set_geometry(NULL) %>% 
  group_by(hexid) %>% 
  summarise(n = n()) %>% 
  right_join(hex_grid, by = "hexid") %>% 
  sf::st_sf()

# calculate the proportion of locations in each bin
ntag <- nrow(pel)
hexbins <- hexbins %>%
  mutate(
    prop = n/ntag
  )

# remove bins without any location data
hexbins_sub <- filter(hexbins, !is.na(n))

# map
ggplot() +

  # location density grid
  geom_sf(data = hexbins_sub, color = "gray79", aes(fill = prop)) +
  scale_fill_gradient(low = "white", high = "#e69b99") +
  
  # coastline
  geom_sf(data = coast, color = NA, fill = "gray19") +
  
  # aesthetics
  coord_sf(crs = 4326,
           xlim = c(-167, -156),
           ylim = c(14.7, 20)) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black", size = 11),
    legend.text = element_text(color = "black", size = 11)
  ) +
  ggspatial::annotation_scale(text_cex = 0.9, text_col = "black") +
  ggspatial::annotation_north_arrow(location = "br")

ggsave(here("outputs","pelagic_original_location_density_map_250kmhexgrid_2025Mar17_v2.png"),
       width = 9, height = 4, units = "in")

