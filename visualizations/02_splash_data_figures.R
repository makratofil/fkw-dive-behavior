## 02_splash_data_figures.R: satellite-linked depth transmitting tag (SPLASH)
## data figures 

## Author: Michaela A. Kratofil, Oregon State University, Cascadia Research
## Updated: 29 Apr 2025

## ---------------------------------------------------------------------------- ##

## load packages
library(dplyr)
library(lubridate)
library(ggplot2)
library(here)

## read in and prep data ## --------------------------------------------------- ##
load(here("pipeline","splash_data_summary_objects_for_figures_2025Apr29.RData"))

## proportion of dives within depth bins for each ID ## ----------------------- ##
dives <- dives %>%
  mutate(
    depth_bin = case_when(
      depth_avg50 <= 100 ~ "100",
      depth_avg50 > 100 & depth_avg50 <= 200 ~ "200",
      depth_avg50 > 200 & depth_avg50 <= 300 ~ "300",
      depth_avg50 > 300 & depth_avg50 <= 400 ~ "400",
      depth_avg50 > 400 & depth_avg50 <= 500 ~ "500",
      depth_avg50 > 500 & depth_avg50 <= 600 ~ "600",
      depth_avg50 > 600 & depth_avg50 <= 700 ~ "700",
      depth_avg50 > 700 & depth_avg50 <= 800 ~ "800",
      depth_avg50 > 800 & depth_avg50 <= 900 ~ "900",
      depth_avg50 > 900 & depth_avg50 <= 1000 ~ "1000",
      depth_avg50 > 1000 & depth_avg50 <= 1100 ~ "1100",
      depth_avg50 > 1100 & depth_avg50 <= 1200 ~ "1200",
      depth_avg50 > 1200 & depth_avg50 <= 1300 ~ "1300",
      depth_avg50 > 1300 ~ "1400"
    )
  )
unique(dives$depth_bin)

# make a leveled factor for depth bin
dives$depth_bin <- factor(dives$depth_bin, levels = c("1400","1300","1200",
                                                      "1100","1000","900","800",
                                                      "700","600","500","400",
                                                      "300","200","100"))
summary(dives$depth_bin)

# get the proportion of dives per depth bin within ID 
n_dives <- dives %>%
  group_by(DeployID) %>%
  summarise(n_dives = n())

depth_bin_sum <- dives %>%
  group_by(DeployID, depth_bin) %>%
  summarise(
    n = n()
  ) %>%
  left_join(., n_dives, by = "DeployID") %>%
  mutate(
    prop = n/n_dives
  )

# add population back
depth_bin_sum <- depth_bin_sum %>%
  mutate(
    population = case_when(
      DeployID %in% c("PcTag026","PcTag028","PcTag030","PcTag032","PcTag055",
                      "PcTag074") ~ "MHI",
      DeployID %in% c("PcTag035","PcTag037","PcTag049") ~ "NWHI",
      DeployID %in% c("PcTag090","PcTag092","PcTagP09") ~ "Open-ocean"
    )
  )


# check that the proportions were calcualted correctly
depth_bin_sum %>%
  group_by(DeployID) %>%
  summarise(
    sum = sum(prop)
  )

# now plot 
ggplot(depth_bin_sum, aes(y = depth_bin, x = prop, fill = population)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#015b58","#89689d","#e69b99")) +
  facet_wrap(~DeployID, nrow = 3, ncol = 4) +
  xlab("Proportion of dives") +
  ylab("Depth bin (m)") +
  scale_x_continuous(limits = c(0,0.5),
                     expand = c(0,0),
                     breaks = seq(0.05,0.45, by = 0.1)) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    legend.position = "top",
    legend.text = element_text(color = "black", size = 11),
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black", size = 13),
    strip.text = element_text(color = "black", size = 11),
    strip.background = element_rect(fill = NA)
  )

ggsave(here("outputs","misc_plots","dive_bin_prop_byid_bypopulation_v7.png"), width = 9, height = 8, units = "in")


## dive depth vs duration ## ------------------------------------------------- ##
ggplot(dives, aes(x = dur_mins, y = -depth_avg50, color = population)) +
  geom_point(alpha = 0.75) +
  geom_smooth(method = "lm", color = "black") +
  scale_color_manual(values = c("#015b58","#89689d","#e69b99")) +
  facet_wrap(~DeployID) +
  scale_x_continuous(
    breaks = seq(2.5, 19, by =2.5),
    limits = c(1.75, 20),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    breaks = seq(-1500, 0, by = 250),
    limits = c(-1500,0),
    expand = c(0,0)
  ) +
  ylab("Dive depth (m)") +
  xlab("Dive duration (min)") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black", size = 11),
    axis.title = element_text(color = "black", size = 13),
    strip.text = element_text(color = "black", size = 12),
    strip.background = element_rect(fill = NA)
  )


ggsave(here("outputs","misc_plots","dive_depth_dur_point_lmline_byid_pop_2025Apr27.png"), width = 12, height = 7.5, units = "in")


## dive depth and duration by ID and sex ## ----------------------------------- ##

# order ID by sex 
dives$DeployID_s <- factor(dives$DeployID, levels = c("PcTag090","PcTag092","PcTag026",
                                                      "PcTag032","PcTag035","PcTag037",
                                                      "PcTag055","PcTagP09","PcTag028",
                                                      "PcTag030","PcTag049","PcTag074"))

# dive depth
dds <- ggplot(dives, aes(x = DeployID_s, y = -depth_avg50, fill = sex)) +
  geom_boxplot(alpha = 0.8) +
  scale_fill_manual(values = c("#c26a7a","#11c2b5","#efddcf")) +
  scale_y_continuous(
    breaks = seq(-1500, 0, by = 250),
    limits = c(-1500,0),
    expand = c(0,0)
  ) +
  ylab("Dive depth (m)") +
  xlab("") +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black", size = 11),
    axis.text.x = element_blank(),
    axis.title = element_text(color = "black", size = 13),
    legend.text = element_text(color = "black", size = 12),
    legend.title = element_blank()
  )
dds
ggsave(here("outputs","misc_plots","dive_depth_byid_bysex_v2.png"), width = 7, height = 5, units = "in")

# dive duration 
ddus <- ggplot(dives, aes(x = DeployID_s, y = dur_mins, fill = sex)) +
  geom_boxplot(alpha = 0.8) +
  scale_fill_manual(values = c("#c26a7a","#11c2b5","#efddcf")) +
  scale_y_continuous(
    breaks = seq(2.5,17.5, by = 2.5),
    limits = c(1,20),
    expand = c(0,0)
  ) +
  ylab("Dive duration (min)") +
  xlab("") +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black", size = 11),
    axis.text.x = element_text(angle = 45, vjust = .9, hjust = .95),
    axis.title = element_text(color = "black", size = 13),
    legend.text = element_text(color = "black", size = 11),
    legend.title = element_blank()
  )
ddus
ggsave(here("outputs","misc_plots","dive_dur_byid_bysex_v2.png"), width = 7, height = 5, units = "in")

# arrange 
ggarrange(dds, ddus, ncol = 1, nrow = 2, labels = c("(a)","(b)"), common.legend = T)
ggsave(here("outputs","misc_plots","dive_depth_dur_byid_bysex_v4.png"), width = 7, height = 7, units = "in",
       bg = "white")

## dive rate by diel period ## ----------------------------------------------- ## 
ggplot(merge_tod, aes(x = tod, y = rate)) +
  geom_boxplot() +
  ylab("Dive rate (# dives/hr)") +
  xlab("") +
  labs(fill = "") +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black", size = 11),
    legend.position = "none",
    axis.title = element_text(color = "black", size = 13)
  )

ggsave(here("outputs","misc_plots","dive_rate_tod_boxplot_across_indv_bw_v1.png"), width = 5, height = 5, units = "in")

## dive rate by clock hour and depth bin ## ---------------------------------- ##

## read in cleaned and processed behavior log data for HOUR of day ## 
hod <- readRDS(here("pipeline","clean_data_for_analysis",
                    "all_behavlog_pseudotracks_rerouted20mIso_geoprocessed_split_hod_2025Feb18.rds"))

# quick summary/check
str(hod)
summary(hod)

# get duration in hours
hod <- hod %>%
  mutate(
    dur_hrs = dur_mins/60,
    dur_days = dur_hrs/24
  ) 

# get column for hour of the day
hod$hod <- hour(hod$start_hst)

# sum by 1-hour block for all surface and dive records 
hod_sum <- hod %>%
  group_by(DeployID, hod) %>%
  summarise(
    n = n(),
    hrs = sum(dur_hrs, na.rm = T)
  ) 

# filter out dives
dives_hod <- filter(hod, What == "Dive") # 1508 dives

# categorize dives by different depth bins 
dives_hod <- dives_hod %>%
  mutate(
    depth100 = ifelse(depth_avg50 >= 100, T, NA),
    depth250 = ifelse(depth_avg50 >= 250, T, NA),
    depth500 = ifelse(depth_avg50 >= 500, T, NA),
    depth750 = ifelse(depth_avg50 >= 750, T, NA)
  )


# get the number of dives 50m+ per hour per day for each tag 
dives50_hod <- dives_hod %>%
  group_by(DeployID, hod) %>%
  summarise(
    n_dives50 = n(),
    avg_depth = mean(depth_avg50),
    max_depth = max(depth_avg50),
    avg_dur = mean(dur_mins)
  )

# do the same for 100m+ dive depth bin
dives100_hod <- dives_hod %>%
  filter(., depth100 == T) %>%
  group_by(DeployID, hod) %>%
  summarise(
    n_dives100 = n()
  )

# do the same for 250m+ dive depth bin
dives250_hod <- dives_hod %>%
  filter(., depth250 == T) %>%
  group_by(DeployID, hod) %>%
  summarise(
    n_dives250 = n()
  )
# do the same for 500m+ dive depth bin
dives500_hod <- dives_hod %>%
  filter(., depth500 == T) %>%
  group_by(DeployID, hod) %>%
  summarise(
    n_dives500 = n()
  )
# do the same for 750m+ dive depth bin
dives750_hod <- dives_hod %>%
  filter(., depth750 == T) %>%
  group_by(DeployID, hod) %>%
  summarise(
    n_dives750 = n()
  )

# combine the two; any records without a dive can be assigned 0 for no dives
rate_hod <- left_join(hod_sum, dives50_hod, by = c("DeployID","hod")) %>%
  mutate(
    n_dives50 = ifelse(is.na(n_dives50),0, n_dives50),
    avg_depth = ifelse(is.na(avg_depth),0, avg_depth),
    max_depth = ifelse(is.na(max_depth),0, max_depth),
    avg_dur = ifelse(is.na(avg_dur),0,avg_dur),
    dph50 = n_dives50/hrs
  ) %>%
  left_join(., dives100_hod, by = c("DeployID","hod")) %>%
  mutate(
    n_dives100 = ifelse(is.na(n_dives100),0, n_dives100),
    dph100 = n_dives100/hrs
  ) %>%
  left_join(., dives250_hod, by = c("DeployID","hod")) %>%
  mutate(
    n_dives250 = ifelse(is.na(n_dives250),0, n_dives250),
    dph250 = n_dives250/hrs
  ) %>%
  left_join(., dives500_hod, by = c("DeployID","hod")) %>%
  mutate(
    n_dives500 = ifelse(is.na(n_dives500),0, n_dives500),
    dph500 = n_dives500/hrs
  ) %>%
  left_join(., dives750_hod, by = c("DeployID","hod")) %>%
  mutate(
    n_dives750 = ifelse(is.na(n_dives750),0, n_dives750),
    dph750 = n_dives750/hrs
  ) 

## dive rates: overall mean/med ## -------------------------------------------- ##
rate_all <- rate_hod %>%
  group_by(hod) %>%
  summarise(
    n_d50 = sum(n_dives50),
    mean_rate50 = mean(dph50),
    med_rate50 = median(dph50),
    sd_rate50 = sd(dph50),
    n_d100 = sum(n_dives100),
    mean_rate100 = mean(dph100),
    med_rate100 = median(dph100),
    sd_rate100 = sd(dph100),
    n_d250 = sum(n_dives250),
    mean_rate250 = mean(dph250),
    med_rate250 = median(dph250),
    sd_rate250 = sd(dph250),
    n_d500 = sum(n_dives500),
    mean_rate500 = mean(dph500),
    med_rate500 = median(dph500),
    sd_rate500 = sd(dph500),
    n_d750 = sum(n_dives750),
    mean_rate750 = mean(dph750),
    med_rate750 = median(dph750),
    sd_rate750 = sd(dph750)
  )

summary(rate_all)

# 50m+ dive clock plot
rate50 <- ggplot(rate_all, aes(x = as.factor(hod), y = mean_rate50)) +
  geom_hline(yintercept = c(0,.25,0.5,.8), colour = "lightgrey")+
  geom_bar(stat = "identity", fill = "azure2",color = "darkgrey") +
  scale_y_continuous(
    limits= c(-0.25,.8),
    breaks = c(0,.25,0.5,.75)) +
  geom_vline(xintercept = 0:24 -.5, colour = "lightgrey")+
  theme_bw() +
  ggtitle("50m+") +
  labs(fill = "") +
  theme(axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=11),
        axis.title.y = element_text(size = 13),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_blank(),
        title = element_text(size = 12)) +
  xlab("") +
  ylab ("Mean dive rate (#dives/hr)") +
  coord_polar() 

rate50

# 250m+ dive clock plot
rate250 <- ggplot(rate_all, aes(x = as.factor(hod), y = mean_rate250)) +
  geom_hline(yintercept = c(0,0.25,.5), colour = "lightgrey")+
  geom_bar(stat = "identity", fill = "#81a9ad",color = "darkgrey") +
  scale_y_continuous(
    limits= c(-0.15,.5),
    breaks = c(0,0.25,.5)) +
  geom_vline(xintercept = 0:24 -.5, colour = "lightgrey")+
  theme_bw() +
  ggtitle("250m+") +
  labs(fill = "") +
  theme(axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=11),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 13),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_blank(),
        title = element_text(size = 12)) +
  xlab("") +
  ylab ("Mean dive rate (#dives/hr)") +
  coord_polar() 
rate250

# 500m+ dives clock plot
rate500 <- ggplot(rate_all, aes(x = as.factor(hod), y = mean_rate500)) +
  geom_hline(yintercept = c(0,.1,0.2), colour = "lightgrey")+
  geom_bar(stat = "identity", fill = "#537380",color = "darkgrey") +
  scale_y_continuous(
    limits= c(-0.06,.2),
    breaks = c(0,0.1,.2)) +
  geom_vline(xintercept = 0:24 -.5, colour = "lightgrey")+
  theme_bw() +
  ggtitle("500m+") +
  labs(fill = "") +
  theme(axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=11),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 13),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_blank(),
        title = element_text(size = 12)) +
  xlab("") +
  ylab ("Mean dive rate (#dives/hr)") +
  coord_polar() 
rate500

# 750m+ dives clock plot
rate750 <- ggplot(rate_all, aes(x = as.factor(hod), y = mean_rate750)) +
  geom_hline(yintercept = c(0,.05,0.1,.15), colour = "lightgrey")+
  geom_bar(stat = "identity", fill = "#33454e",color = "darkgrey") +
  scale_y_continuous(
    limits= c(-0.04,.15),
    breaks = c(0,.05,0.1,.15)) +
  geom_vline(xintercept = 0:24 -.5, colour = "lightgrey")+
  theme_bw() +
  ggtitle("750m+") +
  labs(fill = "") +
  theme(axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=11),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 13),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_blank(),
        title = element_text(size = 12)) +
  xlab("") +
  ylab ("Mean dive rate (#dives/hr)") +
  coord_polar() 
rate750

# cpmbine all
rate_clocks <- ggarrange(rate50, rate250, rate500, rate750,
                     labels = c("(a)","(b)","(c)","(d)"))
rate_clocks
ggsave(here("outputs","clock_plots","circularplot_dive_rate_mean_by_depth_bin_yaxis_scaled_2025Feb21.png"),
       width = 6, height = 6, units = "in", bg = "white")