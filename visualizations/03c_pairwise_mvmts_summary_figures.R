## 03c_pairwise_mvmts_summary_figures.R: summary figures of pairwise movements,
## distances, etc. for supplemental material

## Author: Michaela A. Kratofil, Oregon State University, Cascadia Research
## Updated: 29 April 2025

## --------------------------------------------------------------------------- ##

## load packages
library(here)
library(dplyr)
library(ggplot2)

## PcTag030-032 distance summaries and plots ## ------------------------------ ##

# load data from pairwise distance processing script
load(here("pipeline","crawl_pseudotracks",
          "PcTag030_PcTag032_pairwise_dists_locs_2025Apr29.RData"))

# quick summary stats
median(pair_final$dist_km)
mean(pair_final$dist_km)
max(pair_final$dist_km)

# make a subset with locations with more > 1km error removed 
pair_final$avgse <- (pair_final$se.mu.x + pair_final$se.mu.y)/2
pair_sub <- filter(pair_final, avgse < 1000)

# recalculate summary stats
median(pair_sub$dist_km)
mean(pair_sub$dist_km)
max(pair_sub$dist_km)

# add dive record indicator to the distances data
pair$timestamp <- pair$start_utc

# all dives
dive_dists <- filter(pair, What == "Dive") %>%
  left_join(., pair_final[,2:7], by = "timestamp") %>%
  st_drop_geometry()

# dives with error < 1km
dive_dists_sub <- filter(pair, What == "Dive") %>%
  left_join(., pair_sub[,2:7], by = "timestamp") %>%
  st_drop_geometry()

# dive depth vs distance: all dives
ggplot(dive_dists, aes(x = dist_km, y = -depth_avg50, shape = DeployID)) +
  geom_point() +
  theme_classic() +
  ylab("Dive depth (m)") +
  xlab("Distance between pair (km)") +
  scale_shape_manual(values = c(1,16)) +
  scale_x_continuous(breaks = seq(0,25, by = 5)) +
  theme(
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    legend.position = "inside",
    legend.position.inside = c(0.85,0.1)
  )

ggsave(here("outputs","dive_profile_plots_pairs_2025Mar24","PcTag030_PcTag032",
            "PcTag030_PcTag032_dist_pair_vs_divedepth_2025Apr22.png"), width = 5,
       height = 4, units = "in")

# dive depth vs distance: restricted dives
ggplot(dive_dists_sub, aes(x = dist_km, y = -depth_avg50, shape = DeployID)) +
  geom_point() +
  theme_classic() +
  ylab("Dive depth (m)") +
  xlab("Distance between pair (km)") +
  scale_shape_manual(values = c(1,16)) +
  scale_x_continuous(breaks = seq(0,25, by = 5)) +
  theme(
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    legend.position = "inside",
    legend.position.inside = c(0.85,0.1)
  )

ggsave(here("outputs","dive_profile_plots_pairs_2025Mar24","PcTag030_PcTag032",
            "PcTag030_PcTag032_dist_pair_vs_divedepth_under1km_error_2025Apr22.png"), width = 5,
       height = 4, units = "in")

# dive duration vs distance: all dives
ggplot(dive_dists, aes(x = dist_km, y = dur_mins, shape = DeployID)) +
  geom_point() +
  theme_classic() +
  ylab("Dive duration (min)") +
  xlab("Distance between pair (km)") +
  scale_shape_manual(values = c(1,16)) +
  scale_x_continuous(breaks = seq(0,25, by = 5)) +
  scale_y_continuous(breaks = seq(2.5,17.5, by = 2.5)) +
  theme(
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    legend.position = "inside",
    legend.position.inside = c(0.85,0.85)
  )

ggsave(here("outputs","dive_profile_plots_pairs_2025Mar24","PcTag030_PcTag032",
            "PcTag030_PcTag032_dist_pair_vs_divedur_2025Apr22.png"), width = 5,
       height = 4, units = "in")

# dive duration vs distance: restricted dives
ggplot(dive_dists_sub, aes(x = dist_km, y = dur_mins, shape = DeployID)) +
  geom_point() +
  theme_classic() +
  ylab("Dive duration (min)") +
  xlab("Distance between pair (km)") +
  scale_shape_manual(values = c(1,16)) +
  scale_x_continuous(breaks = seq(0,25, by = 5)) +
  scale_y_continuous(breaks = seq(2.5,17.5, by = 2.5)) +
  theme(
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    legend.position = "inside",
    legend.position.inside = c(0.85,0.85)
  )

ggsave(here("outputs","dive_profile_plots_pairs_2025Mar24","PcTag030_PcTag032",
            "PcTag030_PcTag032_dist_pair_vs_divedur_under1km_error_2025Apr22.png"), width = 5,
       height = 4, units = "in")

# histogram of distances for dives and locations: all dives 
pair_final <- pair_final %>%
  mutate(
    dive = ifelse(timestamp %in% dive_dists$timestamp, T, F)
  )

ggplot(pair_final, aes(x = dist_km, fill = dive)) +
  geom_histogram( color = NA, binwidth = 1, boundary = 0, closed = "left") +
  scale_fill_manual(values = c("grey70","grey90"), labels = c("Track","Dives")) +
  theme_classic() +
  xlab("Distance between pair (km)") +
  ylab("Frequency") +
  geom_vline(xintercept = 2.8, linetype = "dotted") +
  geom_vline(xintercept = 4.4, linetype = "dashed") +
  scale_x_continuous(expand = c(0,0), limits = c(0,30),
                     breaks = seq(0,30, by = 5)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,700)) +
  theme(
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    legend.position = "inside",
    legend.position.inside = c(0.85,0.85)
  )

ggsave(here("outputs","dive_profile_plots_pairs_2025Mar24","PcTag030_PcTag032",
            "PcTag030_PcTag032_dist_pair_histogram_2025Mar30.png"), width = 5,
       height = 4, units = "in")

# histogram of distances for dives and locations: restricted dives
pair_sub <- pair_sub %>%
  mutate(
    dive = ifelse(timestamp %in% dive_dists_sub$timestamp, T, F)
  )

ggplot(pair_sub, aes(x = dist_km, fill = dive)) +
  geom_histogram( color = NA, binwidth = 1, boundary = 0, closed = "left") +
  scale_fill_manual(values = c("grey70","grey90"), labels = c("Track","Dives")) +
  theme_classic() +
  xlab("Distance between pair (km)") +
  ylab("Frequency") +
  geom_vline(xintercept = 2.1, linetype = "dotted") +
  geom_vline(xintercept = 3.4, linetype = "dashed") +
  scale_x_continuous(expand = c(0,0), limits = c(0,30),
                     breaks = seq(0,30, by = 5)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,700)) +
  theme(
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    legend.position = "inside",
    legend.position.inside = c(0.85,0.85)
  )

ggsave(here("outputs","dive_profile_plots_pairs_2025Mar24","PcTag030_PcTag032",
            "PcTag030_PcTag032_dist_pair_histogram_under1km_error_2025Apr22.png"), width = 5,
       height = 4, units = "in")

## PcTag090-092 distance summaries and plots ## ------------------------------ ##

# load data from pairwise distance processing script
load(here("pipeline","crawl_pseudotracks",
          "PcTag090_PcTag092_pairwise_dists_locs_2025Apr29.RData"))

# quick summary stats
median(pair_final$dist_km)
mean(pair_final$dist_km)
max(pair_final$dist_km)

# make a subset with locations with more > 1km error removed 
pair_final$avgse <- (pair_final$se.mu.x + pair_final$se.mu.y)/2
pair_sub <- filter(pair_final, avgse < 1000)

# recalculate summary stats
median(pair_sub$dist_km)
mean(pair_sub$dist_km)
max(pair_sub$dist_km)

# add dive record indicator to the distances data
pair$timestamp <- pair$start_utc

# all dives
dive_dists <- filter(pair, What == "Dive") %>%
  left_join(., pair_final[,2:7], by = "timestamp") %>%
  st_drop_geometry()

# dives with error < 1km
dive_dists_sub <- filter(pair, What == "Dive") %>%
  left_join(., pair_sub[,2:7], by = "timestamp") %>%
  st_drop_geometry()

# dive depth vs distance: all dives
ggplot(dive_dists, aes(x = dist_km, y = -depth_avg50, shape = DeployID)) +
  geom_point() +
  theme_classic() +
  ylab("Dive depth (m)") +
  xlab("Distance between pair (km)") +
  scale_shape_manual(values = c(1,16)) +
  scale_x_continuous(limits = c(0,85)) +
  theme(
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    legend.position = "inside",
    legend.position.inside = c(0.85,0.1)
  )

ggsave(here("outputs","dive_profile_plots_pairs_2025Mar24","PcTag090_PcTag092",
            "PcTag090_PcTag092_dist_pair_vs_divedepth_2025Apr22.png"), width = 5,
       height = 4, units = "in")

# dive depth vs distance: restricted dives
ggplot(dive_dists_sub, aes(x = dist_km, y = -depth_avg50, shape = DeployID)) +
  geom_point() +
  theme_classic() +
  ylab("Dive depth (m)") +
  xlab("Distance between pair (km)") +
  scale_shape_manual(values = c(1,16)) +
  scale_x_continuous(limits = c(0,85)) +
  theme(
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    legend.position = "inside",
    legend.position.inside = c(0.85,0.1)
  )

ggsave(here("outputs","dive_profile_plots_pairs_2025Mar24","PcTag090_PcTag092",
            "PcTag090_PcTag092_dist_pair_vs_divedepth_error_under1km_2025Apr22.png"), width = 5,
       height = 4, units = "in")

# dive duration vs distance: all dives
ggplot(dive_dists, aes(x = dist_km, y = dur_mins, shape = DeployID)) +
  geom_point() +
  theme_classic() +
  ylab("Dive duration (min)") +
  xlab("Distance between pair (km)") +
  scale_shape_manual(values = c(1,16)) +
  scale_x_continuous(limits = c(0,85)) +
  theme(
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    legend.position = "inside",
    legend.position.inside = c(0.85,0.1)
  )

ggsave(here("outputs","dive_profile_plots_pairs_2025Mar24","PcTag090_PcTag092",
            "PcTag090_PcTag092_dist_pair_vs_divedur_2025Apr22.png"), width = 5,
       height = 4, units = "in")

# dive duration vs distance: restricted dives
ggplot(dive_dists_sub, aes(x = dist_km, y = dur_mins, shape = DeployID)) +
  geom_point() +
  theme_classic() +
  ylab("Dive duration (min)") +
  xlab("Distance between pair (km)") +
  scale_shape_manual(values = c(1,16)) +
  scale_x_continuous(limits = c(0,85)) +
  theme(
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    legend.position = "inside",
    legend.position.inside = c(0.85,0.1)
  )

ggsave(here("outputs","dive_profile_plots_pairs_2025Mar24","PcTag090_PcTag092",
            "PcTag090_PcTag092_dist_pair_vs_divedur_error_under1km_2025Apr22.png"), width = 5,
       height = 4, units = "in")

# histogram of distances for dives and locations: all dives 
pair_final <- pair_final %>%
  mutate(
    dive = ifelse(timestamp %in% dive_dists$timestamp, T, F)
  )

ggplot(pair_final, aes(x = dist_km, fill = dive)) +
  geom_histogram( color = NA, binwidth = 1, boundary = 0, closed = "left") +
  scale_fill_manual(values = c("grey70","grey90"), labels = c("Track","Dives")) +
  theme_classic() +
  xlab("Distance between pair (km)") +
  ylab("Frequency") +
  geom_vline(xintercept = 6.7, linetype = "dotted") +
  geom_vline(xintercept = 9.7, linetype = "dashed") +
  scale_x_continuous(expand = c(0,0), limits = c(0,85),
                     breaks = seq(0,85, by = 15)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,100)) +
  theme(
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    legend.position = "inside",
    legend.position.inside = c(0.85,0.85)
  )

ggsave(here("outputs","dive_profile_plots_pairs_2025Mar24","PcTag090_PcTag092",
            "PcTag090_PcTag092_dist_pair_histogram_2025Mar24v3.png"), width = 5,
       height = 4, units = "in")

# histogram of distances for dives and locations: restricted dives
pair_sub <- pair_sub %>%
  mutate(
    dive = ifelse(timestamp %in% dive_dists_sub$timestamp, T, F)
  )

ggplot(pair_sub, aes(x = dist_km, fill = dive)) +
  geom_histogram( color = NA, binwidth = 1, boundary = 0, closed = "left") +
  scale_fill_manual(values = c("grey70","grey90"), labels = c("Track","Dives")) +
  theme_classic() +
  xlab("Distance between pair (km)") +
  ylab("Frequency") +
  geom_vline(xintercept = 7.8, linetype = "dotted") +
  geom_vline(xintercept = 10.4, linetype = "dashed") +
  scale_x_continuous(expand = c(0,0), limits = c(0,85),
                     breaks = seq(0,85, by = 15)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,100)) +
  theme(
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    legend.position = "inside",
    legend.position.inside = c(0.85,0.85)
  )

ggsave(here("outputs","dive_profile_plots_pairs_2025Mar24","PcTag090_PcTag092",
            "PcTag090_PcTag092_dist_pair_histogram_under1km_error_2025Apr22.png"), width = 5,
       height = 4, units = "in")
