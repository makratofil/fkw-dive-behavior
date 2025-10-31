## 07_near_seafloor_dive_figures.R: figures of near-seafloor diving behavior 
## analysis products 

## Author: Michaela A. Kratofil, Oregon State University, Cascadia Research
## Updated: 19 Oct 2025

## ---------------------------------------------------------------------------- ##

## load packages 
library(dplyr)
library(here)
library(ggplot2)
library(sf)
library(ggspatial)
library(ggridges)
library(ggpubr)

## load data from near seafloor dive analysis script ## ----------------------- ##
load(here("pipeline","data_objects_for_seafloor_dives_analysis_2025Oct19.RData"))

## make plot of standard deviation values for all individuals, all dives ## --- ##

# make deployID an ordered factor
dive_sf_sum$DeployID <- factor(dive_sf_sum$DeployID, levels = c("PcTag026",
                                                                "PcTag028",
                                                                "PcTag030",
                                                                "PcTag032",
                                                                "PcTag055",
                                                                "PcTag074",
                                                                "PcTag095",
                                                                "PcTag099",
                                                                "PcTag035",
                                                                "PcTag037",
                                                                "PcTag049",
                                                                "PcTag096",
                                                                "PcTag097"))

# add population back
dive_sf_sum <- dive_sf_sum %>%
  mutate(
    population = case_when(
      DeployID %in% c("PcTag026","PcTag028","PcTag030","PcTag032","PcTag055",
                      "PcTag074","PcTag095","PcTag099") ~ "MHI",
      DeployID %in% c("PcTag035","PcTag037","PcTag049","PcTag096","PcTag097") ~ "NWHI"
    )
  )

# this is figure s7
ggplot(dive_sf_sum, aes(x = DeployID, y = sd_sf, fill = population)) +
  geom_boxplot(alpha = 0.8) +
  xlab("") +
  ylab("SD of seafloor depth (m)") +
  scale_fill_manual(values = c("#015b58","#89689d")) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black", size = 10),
    axis.text.x = element_text(angle = 45, vjust = .9, hjust = .95),
    axis.title = element_text(color = "black", size = 12),
    legend.position = "none"
  )

ggsave(here("outputs","seafloor_depth_plots","SD_seafloor_depth_boxplot_all_2025Oct19.png"),
       width = 6, height = 5, units = "in")

## make boxplots of dive metrics by tod ## --------------------------------------- ##

# dive depth
ggplot(dives_100sd_abs, aes(x = tod, y = dive_depth, fill = tod)) +
  geom_boxplot(outliers = F) +
  scale_fill_manual(values = c("#f8e3d1","#fefbe9","#f5f5f5","#81a9ad")) +
  geom_jitter(data = dives_100sd_abs, width=0.15, size=2, shape = 21, color = "gray39", aes(x = tod, fill = tod, y = dive_depth)) +
  labs(fill = "") +
  xlab("") +
  ylab("Dive depth (m)") +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black", size = 10),
    legend.position = "none",
    #axis.text.x = element_text(angle = 45, vjust = .9, hjust = .95),
    axis.title = element_text(color = "black", size = 12, face = "bold")
  )

ggsave(here("outputs","seafloor_depth_plots","probable_seafloor_dives_divedepth_tod_boxplot_all_2025Oct19.png"),
       width = 5, height = 3, units = "in", bg = "transparent")

# dive duration
ggplot(dives_100sd_abs, aes(x = tod, y = dive_dur, fill = tod)) +
  geom_boxplot(outliers = F) +
  scale_fill_manual(values = c("#f8e3d1","#fefbe9","#f5f5f5","#81a9ad")) +
  geom_jitter(data = dives_100sd_abs, width=0.15, size=2, shape = 21, color = "gray39", aes(x = tod, fill = tod, y = dive_dur)) +
  labs(fill = "") +
  xlab("") +
  ylab("Dive duration (min)") +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black", size = 10),
    legend.position = "none",
    axis.title = element_text(color = "black", size = 12, face = "bold")
  )

ggsave(here("outputs","seafloor_depth_plots","probable_seafloor_dives_divedur_tod_boxplot_all_2025Oct19.png"),
       width = 5, height = 3, units = "in", bg = "transparent")

# dive rate
ggplot(merge_tod2, aes(x = tod, y = sf_rate, fill = tod)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c("#f8e3d1","#fefbe9","#f5f5f5","#81a9ad")) +
  geom_jitter(data = merge_tod2, width=0.15, size=2, shape = 21, color = "gray39", aes(x = tod, fill = tod, y = sf_rate)) +
  scale_color_manual(values = c("#f8e3d1","#fefbe9","#f5f5f5","#81a9ad")) +
  xlab("") +
  ylab("Dive rate (# dives/h)") +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black", size = 10),
    legend.position = "none",
    axis.title = element_text(color = "black", size = 12, face = "bold")
  )

ggsave(here("outputs","seafloor_depth_plots","dive_rates_seafloor_dives_by_tod_boxplot_2025Oct19.png"),
       width = 5, height = 3, units = "in")


## map dives with a high probability of being close to the seafloor ## -------- ##

# read in the regular pseudotrack to get mean estimated position 
pseu <- readRDS(here("pipeline","clean_data_for_analysis",
                     "all_behavlog_pseudotracks_rerouted20mIso_geoprocessed_split_tod_2025Oct17.rds"))

# assign dive ID for each tag to use to filter the other dies
dives_pseu <- filter(pseu, What == "Dive") %>%
  group_by(DeployID) %>%
  mutate(
    dive_id2 = paste0(DeployID,"-",row_number())
  ) %>%
  ungroup()

# get sf points of probable seafloor dives for mapping
dd_sd100_sf <- dives_pseu %>%
  filter(., dive_id2 %in% dives_100sd_abs$dive_id2) %>%
  st_as_sf(., coords = c("lon","lat"), crs = 4326)

# make time of day an ordered factor 
dd_sd100_sf$tod <- factor(dd_sd100_sf$tod, levels = c("dawn","day","dusk","night"))

# get the seafloor depth summaries 
dd_sd100_sf <- left_join(dd_sd100_sf, dives_100sd_abs[,c(1,2,7:15)], by = c("DeployID","dive_id2"))

# quick map of the deep dives 
deep <- filter(dd_sd100_sf, depth_avg50 > 1000)
mapview::mapview(deep)

## map objects ## ------------------------------------------------------------- ##

# custom ggplot theme 
theme_map <- function() {
  theme_bw() +
    theme(panel.background = element_rect(fill = 'white', colour = 'black', size = 1.25),
          axis.text = element_text(colour = 'black', size = 11),
          plot.title = element_text(colour = 'black', face = 'bold'),
          rect = element_rect(color = 'black', size = 1.25)) #+
  
}

# Ocean basemap 
esri_ocean <- paste0('https://services.arcgisonline.com/arcgis/rest/services/',
                     'Ocean/World_Ocean_Base/MapServer/tile/${z}/${y}/${x}.jpeg')

# islands shapefile 
coast <- st_read(here("data","shapefiles","FisheriesIslands.shp"))

## create gridded map ## ----------------------------------------------------- ##

# make hexagonal grid with target area of X squared kilometers 
cell_area <- units::as_units(150, "km^2")
grid_spacing <- sqrt(2*cell_area/sqrt(3))

# transform to UTM
dd_sd100_sf2 <- st_transform(dd_sd100_sf, crs = 3857)

hex_grid <- sf::st_make_grid(
  st_as_sfc(st_bbox(dd_sd100_sf2)),
  cellsize = grid_spacing,
  square = F,
  crs = st_crs(dd_sd100_sf2)) %>% st_sf() %>%
  tibble::rowid_to_column(var = "hexid")

# now need to join points with hex grid and count the number of unique dives 
# in the cell
hexbins <- st_join(dd_sd100_sf2, hex_grid, join = st_intersects) %>% 
  sf::st_set_geometry(NULL) %>% 
  group_by(hexid) %>% 
  summarise(n = n(),
            mean_depth = mean(depth_avg50),
            mean_prop = mean(prop_dive_sf),
            mean_dur = mean(dur_mins)) %>% 
  right_join(hex_grid, by = "hexid") %>% 
  sf::st_sf()

# remove those without any dives
hexbins <- filter(hexbins, !is.na(n))
hexbins$prop <- hexbins$n/133

# color palette
pal <- PNWColors::pnw_palette("Shuksan",100)
mush <- PNWColors::pnw_palette("Mushroom",100)

# number of dives oahu to hi
hioahu_n <- ggplot() +
  # ESRI ocean basemap
  annotation_map_tile(type = esri_ocean, zoomin = 1, progress = "none") +
  
  # hexbins
  geom_sf(data = hexbins, aes(fill = n), color = "gray34") +
  
  # coastline
  geom_sf(data = coast, color = "gray50", fill = "gray19") +
  
  # aesthetics
  scale_fill_gradientn(colours = rev(pal)) +
  labs(fill = "# Dives") +
  coord_sf(
    xlim = c(-159,-154.5),
    ylim = c(18.7,22.2),
    crs = 4326,
    expand = F
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black", size = 10),
    #legend.position = c(.1,.3),
    legend.position = "none",
    legend.text = element_text(size = 10, color = "black"),
    legend.background = element_rect(fill = "white", color = "black")
  ) +
  ggspatial::annotation_north_arrow(location = "tr") +
  ggspatial::annotation_scale()

hioahu_n

# number of dives kauai to nihoa
kaunihoa_n <- ggplot() +
  # ESRI ocean basemap
  annotation_map_tile(type = esri_ocean, zoomin = 1, progress = "none") +
  
  # hexbins
  geom_sf(data = hexbins, aes(fill = n), color = "gray34") +
  
  # coastline
  geom_sf(data = coast, color = "gray50", fill = "gray19") +
  
  # aesthetics
  scale_fill_gradientn(colours = rev(pal)) +
  labs(fill = "# Dives") +
  coord_sf(
    xlim = c(-162.6,-159),
    ylim = c(21.5,23.4),
    crs = 4326,
    expand = F
  ) +
  scale_x_continuous(breaks = c(-162,-161,-160,-159)) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black", size = 10),
    #legend.position = c(.1,.3),
    #legend.position = "none",
    legend.text = element_text(size = 10, color = "black"),
    legend.background = element_rect(fill = "white", color = "black")
  ) +
  ggspatial::annotation_north_arrow(location = "tr") +
  ggspatial::annotation_scale()

kaunihoa_n

# mean dive depth HI to Oahu
hioahu_depth <- ggplot() +
  # ESRI ocean basemap
  annotation_map_tile(type = esri_ocean, zoomin = 1, progress = "none") +
  
  # hexbins
  geom_sf(data = hexbins, aes(fill = -mean_depth), color = "gray34") +
  
  # coastline
  geom_sf(data = coast, color = "gray50", fill = "gray19") +
  
  # aesthetics
  scale_fill_viridis_c(option = "G", direction = 1) +
  labs(fill = stringr::str_wrap("Mean dive depth (m)", 10)) +
  coord_sf(
    xlim = c(-159,-154.5),
    ylim = c(18.7,22.2),
    crs = 4326,
    expand = F
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black", size = 10),
    #legend.position = c(.1,.3),
    legend.position = "none",
    legend.text = element_text(size = 10, color = "black"),
    legend.background = element_rect(fill = "white", color = NA)
  ) +
  ggspatial::annotation_north_arrow(location = "tr") +
  ggspatial::annotation_scale()

hioahu_depth

# ggsave(here("outputs","seafloor_depth_plots",
#             "hexgrid_150km_map_mean_depth_n_hioahu_probable_all_wlegend_2025Apr25.png"),
#        width = 6, height = 5, units = "in", bg = "white")

# now kauai to nihoa
kaunihoa_depth <- ggplot() +
  # ESRI ocean basemap
  annotation_map_tile(type = esri_ocean, zoomin = 1, progress = "none") +
  
  # hexbins
  geom_sf(data = hexbins, aes(fill = -mean_depth), color = "gray34") +
  
  # coastline
  geom_sf(data = coast, color = "gray50", fill = "gray19") +
  
  # aesthetics
  scale_fill_viridis_c(option = "G", direction = 1) +
  labs(fill = "Mean dive depth (m)") +
  coord_sf(
    xlim = c(-162.6,-159),
    ylim = c(21.5,23.4),
    crs = 4326,
    expand = F
  ) +
  scale_x_continuous(breaks = c(-162,-161,-160,-159)) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black", size = 10),
    #legend.position = c(.1,.3),
    legend.position = "none",
    legend.text = element_text(size = 10, color = "black"),
    legend.background = element_rect(fill = "white", color = "black")
  ) +
  ggspatial::annotation_north_arrow(location = "tr") +
  ggspatial::annotation_scale()

kaunihoa_depth

# do proportion of water column 
hioahu_prop <- ggplot() +
  # ESRI ocean basemap
  annotation_map_tile(type = esri_ocean, zoomin = 1, progress = "none") +
  
  # hexbins
  geom_sf(data = hexbins, aes(fill = mean_prop), color = "gray34", alpha = 0.7) +
  
  # coastline
  geom_sf(data = coast, color = "gray19", fill = "gray19") +
  
  # aesthetics
  #scale_fill_gradientn(colours = rev(mush)) +
  scale_fill_viridis_c(option = "F", direction = -1) +
  labs(fill = stringr::str_wrap("Mean ratio dive : seafloor depth", 10)) +
  coord_sf(
    xlim = c(-159,-154.5),
    ylim = c(18.7,22.2),
    crs = 4326,
    expand = F
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black", size = 10),
    #legend.position = c(.1,.3),
    #legend.position = "none",
    legend.text = element_text(size = 10, color = "black"),
    legend.background = element_rect(fill = "white", color = NA)
  ) +
  ggspatial::annotation_north_arrow(location = "tr") +
  ggspatial::annotation_scale()

hioahu_prop 

ggsave(here("outputs","seafloor_depth_plots",
            "hexgrid_150km_map_mean_prop_watercol_hioahu_probable_all_wlegend_2025Oct19.png"),
       width = 6, height = 5, units = "in", bg = "white")

# combine them all
ohis <- ggarrange(hioahu_n, hioahu_depth,
                  nrow = 2, ncol = 1)

kanis <- ggarrange(kaunihoa_n, kaunihoa_depth,
                   nrow = 2, ncol = 1)

als2 <- ggarrange(ohis, kanis, nrow = 1, ncol = 2)
als2

ggsave(here("outputs","seafloor_depth_plots",
            "hexgrid_150km_map_mean_depth_n_hioahu_kaunihoa_probable_all_2025Oct19.png"),
       width = 12, height = 8, units = "in", bg = "white")


## make map of imputed tracks and depth distributions ## --------------------- ##

# filter seafloor dives for pctag074
p74_sf_dives <- filter(dd_sd100_sf, DeployID == "PcTag074")

# get the rest of the dives 
pc74_dives <- dives_pseu %>%
  filter(DeployID == "PcTag074") %>%
  st_as_sf(., coords = c("lon","lat"), crs = 4326)

# imputed tracklines
pc74_lines <- imputes %>%
  filter(., DeployID == "PcTag074") %>%
  st_as_sf(., coords = c("lon","lat"), crs = 4326) %>%
  group_by(trk_id) %>%
  summarise(do_union = F) %>%
  st_cast("LINESTRING")

# map: focusing on a specific area and time 
ggplot() +
  annotation_map_tile(type = esri_ocean, zoomin = 1, progress = "none") +
  geom_sf(data = coast, fill = "gray19", color = NA) +
  geom_sf(data = pc74_lines, color = "#5d74a5", linewidth = .5) +
  geom_sf(data = pc74_dives, fill = "#fbdfa2", color = "black", shape = 21, size = 2.75) +
  coord_sf(crs = 4326,
           xlim = c(-158.2,-157.85),
           ylim = c(21.65,21.83)) +
  scale_x_continuous(breaks = c(-157.9, -158,-158.1,-158.2)) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black", size = 10)
  ) +
  ggspatial::annotation_north_arrow(location = "tr") +
  ggspatial::annotation_scale()

ggsave(here("outputs","seafloor_depth_plots","PcTag074_imptrack_with_dives_n_oahu.png"), width = 6, height = 3, units = "in")  

# subset dives 9:13 to plot with distributions
pc74_dives_imp_sub <- filter(dives, DeployID == "PcTag074") %>%
  filter(., dive_id %in% c(9:13))

# get the dive depths for each of the dives
pc74_dives_sub <- pc74_dives_imp_sub %>%
  group_by(dive_id) %>%
  slice(1)

# create ridgeline plot for each dive id
ggplot(pc74_dives_imp_sub, aes(x = abs(sfdepth), y = as.factor(dive_id))) +
  geom_density_ridges(alpha = 0.65, fill = "#5d74a5", color = NA) +
  geom_point(data = pc74_dives_sub, aes(x = depth_avg50, y = as.factor(dive_id)), shape = 21, color = "gray19", fill = "#fbdfa2",
             size = 3) +
  ylab("Dive ID") +
  xlab("Seafloor depth (m)") +
  scale_x_continuous(breaks = seq(100,500, by = 100),
                     limits = c(100,525)) +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black", size = 11),
    axis.title = element_text(color = "black", size = 12),
    panel.grid.major.y = element_line(color = "gray79")
  )

ggsave(here("outputs","seafloor_depth_plots","PcTag074_ridgeline_depthdiff_n_oahu_v4.png"), width = 6, height = 3, units = "in")

