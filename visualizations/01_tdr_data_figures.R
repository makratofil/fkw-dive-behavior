## 01_tdr_data_figures.R: time-depth recorder (TDR) data figures 

## Author: Michaela A. Kratofil, Oregon State University, Cascadia Research
## Updated: 29 Apr 2025

## --------------------------------------------------------------------------- ##

## load packages 
library(dplyr)
library(lubridate)
library(ggplot2)
library(here)

## vertical & velocity profile figure ## ------------------------------------- ##

# going to focus on PcTDR02, individual that was observed chasing mahimahi soon
# after tagging
tdr <- read.csv(here("data","PcTDR02_tdr_depth_m_for_r.csv"))
summary(tdr)

## plot both depth and relative velocity ## ---------------------------------- ##
ggplot(tdr, aes(x = Time)) +
  annotate(geom = "rect", xmin = 17.65, xmax = 24.2, ymin = -Inf, ymax = Inf,
            fill = "grey", alpha = 0.5, color = NA) +
  geom_line(aes(y = depth_m/10), color = "#015b58", linewidth = 0.75) +
  geom_line(aes(y = Velocity), color = "#015b58", alpha = 0.45) +
  geom_hline(yintercept = 0, color = "black") +
  scale_y_continuous(
    breaks = c(6,5,4,3,2,1,0,-1,-2,-3,-4,-5),
    labels = c("6","5","4","3","2","1","0","-10","-20","-30","-40","-50")
  ) +
  scale_x_continuous(breaks = seq(10,24, by = 1),
                     limits = c(9.9,24.2),
                     expand = c(0,0)) +
  ylab("") +
  xlab("") +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black", size = 11),
    axis.title = element_text(color = "black", size = 13),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 


ggsave(here("outputs","tdr_plots","PcTDR02_dive_profile_plot_depth_velocity_v4.png"),
       width = 8, height = 5, units = "in")

## zoom in to the mahimahi chase ## ------------------------------------------ ##
ggplot(tdr, aes(x = Time)) +
  geom_rect(aes(xmin = 17.65, xmax = 24.2, ymin = -Inf, ymax = Inf),
            fill = "lightgrey", alpha = 0.5, inherit.aes = F) +
  geom_line(aes(y = depth_m/10), color = "#015b58", linewidth = 0.75) +
  geom_line(aes(y = Velocity), color = "#015b58", alpha = 0.45) +
  geom_hline(yintercept = 0, color = "black") +
  scale_y_continuous(
    breaks = c(6,5,4,3,2,1,0,-1,-2,-3,-4,-5),
    labels = c("6","5","4","3","2","1","0","-10","-20","-30","-40","-50")
  ) +
  scale_x_continuous(breaks = seq(10.3,10.7, by = .1),
                     limits = c(10.3,10.701),
                     expand = c(0,0)) +
  ylab("") +
  xlab("") +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black", size = 11),
    axis.title = element_text(color = "black", size = 13),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

ggsave(here("outputs","tdr_plots","PcTDR02_dive_profile_plot_depth_velocity_mahimahi_v2.png"),
       width = 8, height = 5, units = "in")

## supplemental data figures ## ---------------------------------------------- ##

# read in compiled TDR data
tdr <- read.csv(here("data","PcTagTDRDataCompiled.csv"))

## proportion of dives in depth bin plots ## --------------------------------- ##

# create variable for depth bin
tdr <- tdr %>%
  mutate(
    depth_bin = case_when(
      DepthMeter <= 5 ~ "0-5",
      DepthMeter > 5 & DepthMeter <= 10 ~ "6-10",
      DepthMeter > 10 & DepthMeter <= 15 ~ "11-15",
      DepthMeter > 15 & DepthMeter <= 20 ~ "16-20",
      DepthMeter > 20 & DepthMeter <= 25 ~ "21-25",
      DepthMeter > 25 & DepthMeter <= 30 ~ "26-30",
      DepthMeter > 30 & DepthMeter <= 35 ~ "31-35",
      DepthMeter > 35 & DepthMeter <= 40 ~ "36-40",
      DepthMeter > 40 & DepthMeter <= 45 ~ "41-45",
      DepthMeter > 45 & DepthMeter <= 50 ~ "46-50",
      DepthMeter > 50 ~ "50+"
    )
  )
unique(tdr$depth_bin)

# make a leveled factor for depth bin
tdr$depth_bin <- factor(tdr$depth_bin, levels = c("50+","46-50","41-45",
                                                  "36-40","31-35","26-30","21-25",
                                                  "16-20","11-15","6-10","0-5"))
summary(tdr$depth_bin)

# get the proportion of dives per depth bin
depth_bin_sum <- tdr %>%
  group_by(depth_bin) %>%
  summarise(
    n = n(),
    prop = (n/5215)*100
  ) 

# make the plot 
ggplot(depth_bin_sum, aes(y = depth_bin, x = prop)) +
  geom_col(color = "black", fill = NA) +
  theme_classic() +
  ylab("Depth bin (m)") +
  xlab("Percent time") +
  scale_x_continuous(expand = c(0,0),
                     limits = c(0,75),
                     breaks = seq(0,75, by = 5)) +
  theme(
    axis.text = element_text(color = "black", size = 11),
    axis.title = element_text(color = "black", size = 13)
  )

ggsave(here("outputs","tdr_plots","percent_time_depthbin_v1.png"), width = 5, height = 5, units = "in")


## make boxplot by ID for ids with more data across time of day ## ----------- ##

# filter tags that have data over multiple diel periods
sub <- filter(tdr, TagID %in% c("PcTDR02","PcTDR04","PcTDR05"))

# filter dives >= 10 metrs, and exclude those that max the depth resolution
sub10 <- filter(sub, DepthMeter >= 10) %>%
  filter(., DepthMeter < 230)

# dive depth
dd <- ggplot(sub10, aes(x = TimeOfDay, y = -DepthMeter)) +
  geom_boxplot() +
  theme_bw() +
  facet_wrap(~TagID, scales = "free") +
  xlab("") +
  ylab("Dive depth (m)") +
  theme(
    axis.text = element_text(size = 11, color = "black"),
    strip.text = element_text(size = 12, color = "black"),
    strip.background = element_rect(fill = NA),
    axis.title = element_text(size = 13, color = "black")
  )

# dive duration
ddur <- ggplot(sub10, aes(x = TimeOfDay, y = DifferenceInTimeNumeric/60)) +
  geom_boxplot() +
  facet_wrap(~TagID, scales = "free") +
  xlab("") +
  ylab("Dive duration (min)") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 11, color = "black"),
    strip.text = element_text(size = 12, color = "black"),
    strip.background = element_rect(fill = NA),
    axis.title = element_text(size = 13, color = "black")
  )
ddur

# ascent/descent rates 
rate_df <- sub10 %>%
  tidyr::pivot_longer(cols = c("Descent","Ascent"),names_to = "ascdesc", values_to = "rate")
rate_df$ascdesc <- factor(rate_df$ascdesc, levels = c("Descent","Ascent"))

dr <- ggplot(rate_df, aes(x = TimeOfDay, y = abs(rate), fill = ascdesc)) +
  geom_boxplot() +
  facet_wrap(~TagID, scales = "free") +
  scale_fill_manual(values = c("white","grey")) +
  scale_color_manual(values = c("white","grey")) +
  xlab("") +
  ylab("Rate (m/s)") +
  labs(fill = "") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 11, color = "black"),
    axis.text = element_text(size = 11, color = "black"),
    strip.text = element_text(size = 12, color = "black"),
    strip.background = element_rect(fill = NA),
    axis.title = element_text(size = 13, color = "black")
  )

ggarrange(dd, ddur, dr, labels = c("(a)","(b)","(c)"),
          nrow = 3, ncol = 1)

# save it 
ggsave(here("outputs","tdr_plots","tdr_tod_depth_dur_rate_plot_2025Mar16.png"),
       width = 8, height = 10, units = "in")

