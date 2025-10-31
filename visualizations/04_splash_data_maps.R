## 04_splash_data_maps.R: make maps of satellite-linked depth transmitting tag
## (SPLASH) data for manuscript 

## Author: Michaela A. Kratofil, Oregon State University, Cascadia Research
## Updated: 15 Oct 2025

## --------------------------------------------------------------------------- ##

## load packages 
library(dplyr)
library(lubridate)
library(sf)
library(ggplot2)
library(ggspatial)
library(here)

## read in and format files for mapping ## ----------------------------------- ##

# behavior log data
beh <- readRDS(here("pipeline","clean_data_for_analysis",
                    "all_behavlog_pseudotracks_rerouted20mIso_geoprocessed_split_tod_2025Oct10.rds"))

# add a population column
beh <- beh %>%
  mutate(
    population = case_when(
      DeployID %in% c("PcTag026","PcTag028","PcTag030","PcTag032","PcTag055",
                      "PcTag074","PcTag095","PcTag099") ~ "MHI",
      DeployID %in% c("PcTag035","PcTag037","PcTag049","PcTag096","PcTag097") ~ "NWHI",
      DeployID %in% c("PcTag090","PcTag092","PcTagP09") ~ "Open-ocean"
    )
  )


# make spatial
dives <- filter(beh, What == "Dive")
dives_sf <- st_as_sf(dives, coords = c("lon","lat"), crs = 4326)

# get tag deployment locations 
deps <- read.csv(here("data","location","PcTag001-099_DouglasFiltered_r20d3lc2_ArgosGPS_2025OCTv1.csv")) %>%
  filter(animal %in% c("PcTag026","PcTag028","PcTag030","PcTag032",
                       "PcTag035","PcTag037","PcTag049","PcTag055",
                       "PcTag074","PcTag090","PcTag092","PcTagP09",
                       "PcTag095","PcTag096","PcTag097","PcTag099")) %>%
  filter(., LC == "DP") %>%
  st_as_sf(., coords = c("longitud","latitude"), crs = 4326) %>%
  mutate(
    population = case_when(
      animal %in% c("PcTag026","PcTag028","PcTag030","PcTag032","PcTag055",
                      "PcTag074","PcTag095","PcTag099") ~ "MHI",
      animal %in% c("PcTag035","PcTag037","PcTag049","PcTag096","PcTag097") ~ "NWHI",
      animal %in% c("PcTag090","PcTag092","PcTagP09") ~ "Open-ocean"
    )
  )

# mote locations 
motes <- read.csv(here("pipeline","fkw_splash_mote_ids_by_deploy_id.csv")) %>%
  filter(., !is.na(Mote.Id)) %>%
  group_by(Mote.Id) %>%
  slice(1) %>%
  st_as_sf(., coords = c("Longitude","Latitude"), crs = 4326)

## mapping objects ## -------------------------------------------------------- ##

# Ocean basemap 
esri_ocean <- paste0('https://services.arcgisonline.com/arcgis/rest/services/',
                     'Ocean/World_Ocean_Base/MapServer/tile/${z}/${y}/${x}.jpeg')

# islands
coast <- st_read(here("data","shapefiles","FisheriesIslands.shp")) %>%
  st_transform(crs = 3857)

## dive locations for NWHI and MHI region ## --------------------------------- ## 
ggplot() +
  
  # ESRI ocean basemap
  annotation_map_tile(type = esri_ocean, zoomin = 1, progress = "none") +
  
  # dive locations
  geom_sf(data = dives_sf, aes(color = population), alpha = 0.8, size = .7) +
  
  # islands/coastline
  geom_sf(data = coast, fill = "gray19", color = "gray19") +
  
  # deployment locations
  geom_sf(data = deps, shape = 21, aes(color = population), fill = "white", size = 1) +
  
  # color of dive locations
  scale_color_manual(values = c("#015b58","#89689d","#e69b99")) +
  
  # MOTE locations
  geom_sf(data = motes, shape = 8, color = "white", size = 1.5) +
  
  # aesthetics
  coord_sf(crs = 4326,
           xlim = c(-163.7,-154.5),
           ylim = c(18.6,24)) +
  theme_bw() +
  labs(fill = "") +
  theme(
    axis.text = element_text(color = "black", size = 10),
    legend.position = "top",
    axis.ticks.length=unit(-0.1, "cm"),
    legend.text = element_text(size = 10, color = "black"),
    legend.background = element_rect(fill = "white", color = "black")
  ) +
  ggspatial::annotation_north_arrow(location = "tr") +
  ggspatial::annotation_scale(location = "br")

ggsave(here("outputs","misc_plots","map_dive_locs_colored_bypop_mhi_nwhi_region_2025Oct10.png"),
       width = 8, height = 6, units = "in")

## PcTagP09 map ## ----------------------------------------------------------- ##
ggplot() +
  # ESRI ocean basemap
  annotation_map_tile(type = esri_ocean, zoomin = 1, progress = "none") +
  
  # coastline
  geom_sf(data = coast, fill = "gray19", color = "gray19") +
  
  # dive locations
  geom_sf(data = dives_sf, aes(fill = population), color = "gray19", shape = 21, alpha = 0.8) +
  
  # deployment location
  geom_sf(data = deps, shape = 21, aes(color = population), fill = "white", size = 1.5) +
  
  # aesthetics
  scale_fill_manual(values = c("#015b58","#89689d","#e69b99")) +
  scale_color_manual(values = c("#015b58","#89689d","#e69b99")) +
  coord_sf(crs = 4326,
           xlim = c(-167,-160.6),
           ylim = c(15,18)) +
  scale_y_continuous(breaks = c(15.5, 16.5, 17.5)) +
  scale_x_continuous(breaks = c(-161,-163,-165,-167)) +
  theme_bw() +
  labs(fill = "") +
  theme(
    axis.text = element_text(color = "black", size = 10),
    panel.background = element_rect(fill = "transparent",
                                    colour = NA_character_), # necessary to avoid drawing panel outline
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    plot.background = element_rect(fill = "transparent",
                                   colour = NA_character_), # necessa
    legend.position = "none",
    axis.ticks.length=unit(-0.1, "cm"),
    legend.text = element_text(size = 10, color = "black"),
    legend.background = element_rect(fill = "white", color = "black")
  ) +
  #ggspatial::annotation_north_arrow(location = "tr") +
  ggspatial::annotation_scale()

ggsave(here("outputs","misc_plots","map_dive_locs_colored_bypop_PcTagP09_region_2025Jan12v3.png"),
       width = 6, height = 3, units = "in", bg = "transparent")

## make gridded data summary maps ## ----------------------------------------- ##

# check the west most extent for new NWHI tags
n <- dives %>%
  filter(DeployID %in% c("PcTag096","PcTag097")) %>%
  summarise(
    minlon = min(lon),
    maxlon = max(lon)
  )

# first transform to UTM
dives_sf <- st_transform(dives_sf, crs = 3857)

# make hexagonal grid with target area of X squared kilometers 
cell_area <- units::as_units(750, "km^2")
grid_spacing <- sqrt(2*cell_area/sqrt(3))

hex_grid <- sf::st_make_grid(
  st_as_sfc(st_bbox(dives_sf)),
  cellsize = grid_spacing,
  square = F,
  crs = st_crs(dives_sf)) %>% st_sf() %>%
  tibble::rowid_to_column(var = "hexid")

# now need to join points with hex grid and count the number of unique dives 
# in the cell
hexbins <- st_join(dives_sf, hex_grid, join = st_intersects) %>% 
  sf::st_set_geometry(NULL) %>% 
  group_by(hexid) %>% 
  summarise(n = n(),
            mean_depth = mean(depth_avg50),
            med_depth = median(depth_avg50),
            cv_depth = sd(depth_avg50)/mean_depth,
            max_depth = max(depth_avg50),
            mean_dur = mean(dur_mins),
            med_dur = median(dur_mins),
            cv_dur = sd(dur_mins)/mean_dur) %>% 
  right_join(hex_grid, by = "hexid") %>% 
  sf::st_sf()

# remove those without any dives
hexbins <- filter(hexbins, !is.na(n))
hexbins$prop <- hexbins$n/1508

# filter out bins that only have 1 dive
hexbins_sub <- filter(hexbins, n > 1) %>%
  st_transform(crs = 4326) # transform back to 4326 (easy to plot, and crop)

# color palettes
pal <- PNWColors::pnw_palette("Shuksan",100)
pal2 <- PNWColors::pnw_palette("Moth",100)

# remove grids with PcTagP09 (focus on where the most data is)
hexbins_sub2 <- st_crop(hexbins_sub, xmin = -164, xmax = -154.5,
                        ymin = 18, ymax = 24)

# save for later 
save(coast, hexbins_sub2, file = here("pipeline","data_for_nwhi_mhi_gridded_dive_maps_750km2_2025Oct10.RData"))

# CV dive depth
ggplot() +
  # hexbins
  geom_sf(data = hexbins_sub2, aes(fill = cv_depth*100), color = "gray64", alpha = 0.9) +
  #geom_sf(data = hexbins_sub2, aes(fill = cv_depth*100), color = NA, alpha = 0.6) +
  
  # coastline
  #geom_sf(data = coast, fill = "gray19", color = "gray19") +
  geom_sf(data = coast, fill = "gray19", color = "gray50") +
  #geom_sf(data = coast, fill = "gray90", color = "gray19") +
  
  # aesthetics
  scale_fill_gradientn(colours = rev(pal2), na.value = NA) +
  coord_sf(crs = 4326,
           xlim = c(-163.8,-154.5),
           ylim = c(18,24)) +
  theme_bw() +
  labs(fill = stringr::str_wrap("CV dive depth (%)",10)) +
  theme(
    axis.text = element_text(color = "black", size = 10),
    legend.position = c(.1,.3),
    axis.ticks.length=unit(-0.1, "cm"),
    legend.text = element_text(size = 10, color = "black"),
    legend.background = element_rect(fill = "white", color = "black")
  ) +
  ggspatial::annotation_north_arrow(location = "tr") +
  ggspatial::annotation_scale(text_cex = 0.9)

ggsave(here("outputs","misc_plots","map_dive_grid_cvdepth_nwhi-mhi_region_2025Oct10v5.png"),
       width = 7.5, height = 5, units = "in")

ggsave(here("outputs","misc_plots","map_dive_grid_cvdepth_nwhi-mhi_region_2025Oct10v5.pdf"),
       width = 7.5, height = 5, units = "in")

# mean dive depth 
ggplot() +
  # hexbins
  geom_sf(data = hexbins_sub2, aes(fill = -mean_depth), color = "gray64", alpha = 0.9) +
  #geom_sf(data = hexbins_sub2, aes(fill = -mean_depth), color = NA, alpha = 0.6) +
  
  # coastline
  geom_sf(data = coast, fill = "gray19", color = "gray50") +
  #geom_sf(data = coast, fill = "gray90", color = "gray19") +
  
  # aesthetics
  scale_fill_viridis_c(option = "G", na.value = NA, direction = 1) +
  coord_sf(crs = 4326,
           xlim = c(-163.8,-154.5),
           ylim = c(18,24)) +
  theme_bw() +
  labs(fill = stringr::str_wrap("Mean dive depth (m)", 10)) +
  theme(
    axis.text = element_text(color = "black", size = 10),
    legend.position = c(.1,.3),
    axis.ticks.length=unit(-0.1, "cm"),
    legend.text = element_text(size = 10, color = "black"),
    legend.background = element_rect(fill = "white", color = "black")
  ) +
  ggspatial::annotation_north_arrow(location = "tr") +
  ggspatial::annotation_scale(text_cex = 0.9)

ggsave(here("outputs","misc_plots","map_dive_grid_meandepth_nwhi-mhi_region_2025Oct10v5.png"),
       width = 7.5, height = 5, units = "in")

ggsave(here("outputs","misc_plots","map_dive_grid_meandepth_nwhi-mhi_region_2025Oct10v5.pdf"),
       width = 7.5, height = 5, units = "in")

# max dive depth
ggplot() +
  # hexbins
  geom_sf(data = hexbins_sub2, aes(fill = -max_depth), color = "gray64", alpha = 0.9) +
  #geom_sf(data = hexbins_sub2, aes(fill = -max_depth), color = NA, alpha = 0.6) +
  
  # coastline
  geom_sf(data = coast, fill = "gray19", color = "gray50") +
  #geom_sf(data = coast, fill = "gray90", color = "gray19") +
  #geom_sf(data = coast, fill = "gray69", color = "gray19") +
  
  # aesthetics
  scale_fill_viridis_c(option = "G", na.value = NA, direction = 1) +
  coord_sf(crs = 4326,
           xlim = c(-163.8,-154.5),
           ylim = c(18,24)) +
  theme_bw() +
  labs(fill = stringr::str_wrap("Max dive depth (m)", 10)) +
  theme(
    axis.text = element_text(color = "black", size = 10),
    #legend.position = "top",
    legend.position = c(.1,.3),
    axis.ticks.length=unit(-0.1, "cm"),
    legend.text = element_text(size = 10, color = "black"),
    legend.background = element_rect(fill = "white", color = "black")
  ) +
  ggspatial::annotation_north_arrow(location = "tr") +
  ggspatial::annotation_scale(text_cex = 0.9)

ggsave(here("outputs","misc_plots","map_dive_grid_maxdepth_nwhi-mhi_region_2025Oct10v5.png"),
       width = 7.5, height = 5, units = "in")

ggsave(here("outputs","misc_plots","map_dive_grid_maxdepth_nwhi-mhi_region_2025Oct10v5.pdf"),
       width = 7.5, height = 5, units = "in")

# number of dives
ggplot() +
  # hexbins
  geom_sf(data = hexbins_sub2, aes(fill = n), color = "gray64") +
  #geom_sf(data = hexbins_sub2, aes(fill = n), color = NA, alpha = 0.99) +
  
  # coastline
  geom_sf(data = coast, fill = "gray19", color = "gray50") +
  #geom_sf(data = coast, fill = "gray90", color = "gray19") +
  #geom_sf(data = coast, fill = "gray69", color = "gray19") +
  
  # aesthetics
  scale_fill_gradientn(colours = rev(pal), na.value = NA) +
  coord_sf(crs = 4326,
           xlim = c(-163.8,-154.5),
           ylim = c(18,24)) +
  theme_bw() +
  labs(fill = stringr::str_wrap("# dives", 10)) +
  theme(
    axis.text = element_text(color = "black", size = 10),
    legend.position = c(.1,.3),
    axis.ticks.length=unit(-0.1, "cm"),
    legend.text = element_text(size = 10, color = "black"),
    legend.background = element_rect(fill = "white", color = "black")
  ) +
  ggspatial::annotation_north_arrow(location = "tr") +
  ggspatial::annotation_scale(text_cex = 0.9)

ggsave(here("outputs","misc_plots","map_dive_grid_ndives_nwhi-mhi_region_2025Oct10v5.png"),
       width = 7.5, height = 5, units = "in")

ggsave(here("outputs","misc_plots","map_dive_grid_ndives_nwhi-mhi_region_2025Oct10v5.pdf"),
       width = 7.5, height = 5, units = "in")
