## 06_gamm_prediction_figures.R: make prediction plots from fitted GAMM models 
## for manuscript

## Author: Michaela A. Kratofil, Oregon State University, Cascadia Research
## Updated: 29 Apr 2025

## --------------------------------------------------------------------------- ##

## load packages 
library(dplyr)
library(here)
library(lubridate)
library(sf)
library(ggplot2)
library(mgcv)
library(gratia)
library(marginaleffects)
library(ggpubr)

## dive depth ## ------------------------------------------------------------- ##

# read in the original data that was formatted in the model script
dives_pop <- readRDS(here("pipeline","gamm_models","spatiotemporal_gamms_allpops_model_ff_data_2025Feb21.rds"))

# read in the spatiotemporal model 
depth_pop_gs <- readRDS(here("pipeline","gamm_models","spatiotemporal_gamms_allpops_ff_depth_pop_gs_2025Feb21.rds"))

## time of day (hms) ## 

# extract predictions from model for time of day variable (global smooth)
depth_pop_gs_hms_preds <- plot_predictions(depth_pop_gs, condition = "hms_n", type = "response",
                                           draw = F)

# get the same curve but from the temporal-only predictors model 
depth_t_gs <- readRDS(here("pipeline","gamm_models","temporal_gamms_allpops_depth_gs_2025Feb20.rds"))
depth_gs_hms_preds <- plot_predictions(depth_t_gs, condition = "hms_n", type = "response",
                                       draw = F)


# time of day, zoomed into the effect  
depo_hms <- ggplot(dives_pop, aes(x = hms_n, y = depth_avg50)) +
  
  # observed points
  geom_rug(sides = "b") +
  
  # confidence ribbon for full/spatiotemporal model
  geom_ribbon(aes(ymin=conf.low, ymax = conf.high,x=hms_n), data = depth_pop_gs_hms_preds,
              alpha = 0.2, color = NA) +
  
  # confidence ribbon for temporal predictors only model
  geom_ribbon(aes(ymin=conf.low, ymax = conf.high,x=hms_n), data = depth_gs_hms_preds,
              fill = NA, color = "gray39", linetype = "dashed") +
  
  # estimated line/relationship for spatiotemporal + temporal-predictors only model
  geom_line(aes(y = estimate), data = depth_pop_gs_hms_preds, linewidth = 1) +
  geom_line(aes(y = estimate), data = depth_gs_hms_preds, linetype = "dashed", linewidth = .75) +
  
  # aesthetics
  theme_bw() +
  scale_y_continuous(expand = c(0,0),
                     limits = c(75,450),
                     breaks = seq(100,400, by = 50)) +
  scale_x_continuous(
    limits = c(0,86400),
    breaks = c(01,21600,43200,64800,86399),
    labels = c("00:01","06:00","12:00","18:00","23:59")) +
  ylab("Dive depth (m)") +
  xlab("Time of day (HST)") +
  theme(
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, color = "black"),
    strip.text = element_text(size = 11)
  )
depo_hms

## time of day x tag ID ## 
depth_pop_gs_hms_id_preds <- plot_predictions(depth_pop_gs, condition = c("hms_n","DeployID"), type = "response",
                                              draw = F)

# color by population
id_pop_pal <- c(rep("#015b58",6),rep("#89689d",3),rep("#e69b99",3))

# get the same curves but from the temporal predictors only model
depth_gs_hms_id_preds <- plot_predictions(depth_t_gs, condition = c("hms_n","DeployID"), type = "response",
                                          draw = F)

depo_id_line <- ggplot(dives_pop, aes(x = hms_n, y = depth_avg50)) +
  
  # observed points
  geom_rug(sides = "b") +
  
  # predicted estimates for spatiotemporal + temporal model
  geom_line(aes(y = estimate, color = DeployID), data = depth_pop_gs_hms_id_preds, linewidth = .75,
            alpha = 0.75) +
  geom_line(aes(y = estimate, color = DeployID), data = depth_gs_hms_id_preds, linewidth = .75,
            alpha = 0.3, linetype = "dashed") +
  
  # aesthetics
  scale_color_manual(values = id_pop_pal) +
  theme_bw() +
  scale_y_continuous(expand = c(0,0),
                     limits = c(50,875),
                     breaks = seq(50,800, by = 100)) +
  scale_x_continuous(
    limits = c(0,86400),
    breaks = c(01,21600,43200,64800,86399),
    labels = c("00:01","06:00","12:00","18:00","23:59")) +
  ylab("Dive depth (m)") +
  xlab("Time of day (HST)") +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, color = "black"),
    strip.text = element_text(size = 11)
  )
depo_id_line

## moon phase ##

# extract predictions from model for moon phase (global smooth)
depth_pop_gs_moon_preds <- plot_predictions(depth_pop_gs, condition = "moon_phase_rads", type = "response",
                                            draw = F)

# same curve from temporal predictors only model
depth_gs_moon_preds <- plot_predictions(depth_t_gs, condition = "moon_phase_rads", type = "response",
                                       draw = F)

## moon phase, zoomed in to the effect ## 
depo_moon <- ggplot(dives_pop, aes(x = moon_phase_rads, y = depth_avg50)) +
  
  # observed points
  geom_rug(sides = "b") +
  
  # confidence ribbons for spatiotemporal + temporal only models
  geom_ribbon(aes(ymin=conf.low, ymax = conf.high,x=moon_phase_rads), data = depth_pop_gs_moon_preds,
              alpha = 0.2, color = NA) +
  geom_ribbon(aes(ymin=conf.low, ymax = conf.high,x=moon_phase_rads), data = depth_gs_moon_preds,
              fill = NA, color = "gray39", linetype = "dashed") +
  
  # predicted relationship for spatiotemporal + temporal only models 
  geom_line(aes(y = estimate), data = depth_pop_gs_moon_preds, linewidth = 1) +
  geom_line(aes(y = estimate), data = depth_gs_moon_preds, linetype = "dashed", linewidth = .75) +
  
  # aesthetics
  theme_bw()+
  scale_y_continuous(expand = c(0,0),
                     limits = c(75,450),
                     breaks = seq(100,400, by = 50)) +
  scale_x_continuous(
    limits = c(0,6.3),
    breaks = c(0,pi/2,pi,3*pi/2,6.3),
    labels = c("New","First quarter","Full","Last quarter","New")) +
  ylab("Dive depth (m)") +
  xlab("Moon phase") +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, color = "black"),
    strip.text = element_text(size = 11)
  )
depo_moon

## slope ## 
depth_pop_gs_slope_preds <- plot_predictions(depth_pop_gs, condition = "sf_slope", type = "response",
                                             draw = F)

## slope, zoomed in to the effect ## 
depo_slope <- ggplot(dives_pop, aes(x = sf_slope, y = depth_avg50)) +
  # observed points
  geom_rug(sides = "b") +
  
  # confidence ribbon
  geom_ribbon(aes(ymin=conf.low, ymax = conf.high,x=sf_slope), data = depth_pop_gs_slope_preds,
              alpha = 0.2, color = NA) +
  
  # predicted estimates
  geom_line(aes(y = estimate), data = depth_pop_gs_slope_preds, linewidth = 1) +
  
  # aesthetics
  theme_bw() +
  scale_y_continuous(expand = c(0,0),
                     limits = c(75,450),
                     breaks = seq(100,400, by = 50)) +
  scale_x_continuous(
    breaks = seq(0,40,  by = 10)) +
  ylab("Dive depth (m)") +
  xlab("Slope (degrees)") +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, color = "black"),
    strip.text = element_text(size = 11)
  )
depo_slope

## chla ## 
depth_pop_gs_chla_preds <- plot_predictions(depth_pop_gs, condition = "logchla", type = "response",
                                            draw = F)

## chla, zoomed in to the effect ## 
depo_chla <- ggplot(dives_pop, aes(x = logchla, y = depth_avg50)) +
  # observed points
  geom_rug(sides = "b") +
  
  # confidence ribbon
  geom_ribbon(aes(ymin=conf.low, ymax = conf.high,x=logchla), data = depth_pop_gs_chla_preds,
              alpha = 0.2, color = NA) +
  
  # predicted estimates
  geom_line(aes(y = estimate), data = depth_pop_gs_chla_preds, linewidth = 1) +
  
  # aesthetics
  theme_bw() +
  scale_y_continuous(expand = c(0,0),
                     limits = c(50,600),
                     breaks = seq(50,600, by = 100)) +
  ylab("Dive depth (m)") +
  xlab(bquote('log10(Chlorophyll-a ('*mg^3*'))')) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, color = "black"),
    strip.text = element_text(size = 11)
  )
depo_chla

## mld ## 
depth_pop_gs_mld_preds <- plot_predictions(depth_pop_gs, condition = "mld", type = "response",
                                            draw = F)

## mld but zoomed in to the effect ## 
depo_mld <- ggplot(dives_pop, aes(x = mld, y = depth_avg50)) +
  # observed points
  geom_rug(sides = "b") +
  
  # confidence ribbon
  geom_ribbon(aes(ymin=conf.low, ymax = conf.high,x=mld), data = depth_pop_gs_mld_preds,
              alpha = 0.2, color = NA) +
  
  # predicted estimates
  geom_line(aes(y = estimate), data = depth_pop_gs_mld_preds, linewidth = 1) +
  
  # aesthetics
  theme_bw() +
  scale_y_continuous(expand = c(0,0),
                     limits = c(75,800),
                     breaks = seq(100,750, by = 100)) +
  ylab("Dive depth (m)") +
  xlab("Mixed layer depth (m)") +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, color = "black"),
    strip.text = element_text(size = 11)
  )
depo_mld

# arrange
depo <- ggarrange(depo_hms, depo_id_line, depo_moon, depo_slope, depo_chla, depo_mld)
depo

ggsave(here("outputs","gamm_plots","spatiotemporal_gamm_population_P09",
            "depth_gamm_spatiotemporal_population_wP09_gs_ff_all_predplot_2025Mar15.png"),
       width = 12, height = 6, units = "in")


## dive duration ## ---------------------------------------------------------- ## 

# read in full spatiotemporal model
dur_pop_gs <- readRDS(here("pipeline","gamm_models","spatiotemporal_gamms_allpops_ff_duration_pop_gs_2025Feb21.rds"))

# read in the temporal predictors only model
dur_t_gs <- readRDS(here("pipeline","gamm_models","temporal_gamms_allpops_dur_secs_gs_2025Feb20.rds"))

## hms/time of day ## 

# extract predictions from model for hms (global smooth)
dur_pop_gs_hms_preds <- plot_predictions(dur_pop_gs, condition = "hms_n", type = "response",
                                          draw = F)

# get the same curve from the temporal predictors only model
dur_gs_hms_preds <- plot_predictions(dur_t_gs, condition = "hms_n", type = "response",
                                         draw = F)

dupo_hms <- ggplot(dives_pop, aes(x = hms_n, y = dur_secs)) +
  # observed points
  geom_rug(sides = "b") +
  
  # confidence ribbons 
  geom_ribbon(aes(ymin=conf.low, ymax = conf.high,x=hms_n), data = dur_pop_gs_hms_preds,
              alpha = 0.2, color = NA) +
  geom_ribbon(aes(ymin=conf.low, ymax = conf.high,x=hms_n), data = dur_gs_hms_preds,
              fill = NA, linetype = "dashed", color = "gray39") +
  
  # predicted estimates
  geom_line(aes(y = estimate), data = dur_pop_gs_hms_preds, linewidth = 1) +
  geom_line(aes(y = estimate), data = dur_gs_hms_preds, linewidth = .75,
            linetype = "dashed") +
  
  # aesthetics
  theme_bw() +
  scale_y_continuous(expand = c(0,0),
                     limits = c(175,440),
                     breaks = seq(180,450,by=60),
                     labels = c("3","4","5","6","7")) +
  scale_x_continuous(
    limits = c(0,86400),
    breaks = c(01,21600,43200,64800,86399),
    labels = c("00:01","06:00","12:00","18:00","23:59")) +
  ylab("Dive duration (min)") +
  xlab("Time of day (HST)") +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, color = "black"),
    strip.text = element_text(size = 11)
  )
dupo_hms

## time of day x tag ID ## 
dur_pop_gs_hms_id_preds <- plot_predictions(dur_pop_gs, condition = c("hms_n","DeployID"), type = "response",
                                              draw = F)

# same one from temp model
dur_gs_hms_id_preds <- plot_predictions(dur_t_gs, condition = c("hms_n","DeployID"), type = "response",
                                            draw = F)

# color by population
id_pop_pal <- c(rep("#015b58",6),rep("#89689d",3),rep("#e69b99",3))

dupo_id_line <- ggplot(dives_pop, aes(x = hms_n, y = dur_secs)) +
  
  # observed points
  geom_rug(sides = "b") +
  
  # predicted estimates for spatiotemporal and temporal only models
  geom_line(aes(y = estimate, color = DeployID), data = dur_pop_gs_hms_id_preds, linewidth = .75,
            alpha = 0.75) +
  geom_line(aes(y = estimate, color = DeployID), data = dur_gs_hms_id_preds, linewidth = .75,
            alpha = 0.3, linetype = "dashed") +
  
  # aesthetics
  scale_color_manual(values = id_pop_pal) +
  theme_bw() +
  scale_y_continuous(expand = c(0,0),
                     limits = c(175,720),
                     breaks = seq(180,660,by=120),
                     labels = c("3","5","7","9","11")
                     ) +
  scale_x_continuous(
    limits = c(0,86400),
    breaks = c(01,21600,43200,64800,86399),
    labels = c("00:01","06:00","12:00","18:00","23:59")) +
  ylab("Dive duration (min)") +
  xlab("Time of day (HST)") +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, color = "black"),
    strip.text = element_text(size = 11)
  )

dupo_id_line

## moon phase ##

# extract predictions from model for moon phase (global smooth)
dur_pop_gs_moon_preds <- plot_predictions(dur_pop_gs, condition = "moon_phase_rads", type = "response",
                                            draw = F)

# same from temp model
dur_gs_moon_preds <- plot_predictions(dur_t_gs, condition = "moon_phase_rads", type = "response",
                                          draw = F)

## moon phase, zoomed in to the effect ## 
dupo_moon <- ggplot(dives_pop, aes(x = moon_phase_rads, y = dur_secs)) +
  # observed points
  geom_rug(sides = "b") +
  
  # confidence ribbons
  geom_ribbon(aes(ymin=conf.low, ymax = conf.high,x=moon_phase_rads), data = dur_pop_gs_moon_preds,
              alpha = 0.2, color = NA) +
  geom_ribbon(aes(ymin=conf.low, ymax = conf.high,x=moon_phase_rads), data = dur_gs_moon_preds,
              fill = NA, color = "gray39", linetype = "dashed") +
  
  # predicted estimates
  geom_line(aes(y = estimate), data = dur_pop_gs_moon_preds, linewidth = 1) +
  geom_line(aes(y = estimate), data = dur_gs_moon_preds, linewidth = .75, linetype = "dashed") +
  
  # aesthetics
  theme_bw()+
  scale_y_continuous(expand = c(0,0),
                     limits = c(175,440),
                     breaks = seq(180,450,by=60),
                     labels = c("3","4","5","6","7")) +
  scale_x_continuous(
    limits = c(0,6.3),
    breaks = c(0,pi/2,pi,3*pi/2,6.3),
    labels = c("New","First quarter","Full","Last quarter","New")) +
  ylab("Dive duration (min)") +
  xlab("Moon phase") +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, color = "black"),
    strip.text = element_text(size = 11)
  )
dupo_moon

## slope ## 
dur_pop_gs_slope_preds <- plot_predictions(dur_pop_gs, condition = "sf_slope", type = "response",
                                             draw = F)

## slope, zoomed in to the effect ## 
dupo_slope <- ggplot(dives_pop, aes(x = sf_slope, y = dur_secs)) +
  # observed points
  geom_rug(sides = "b") +
  
  # confidence ribbon
  geom_ribbon(aes(ymin=conf.low, ymax = conf.high,x=sf_slope), data = dur_pop_gs_slope_preds,
              alpha = 0.2, color = NA) +
  
  # predicted estimates
  geom_line(aes(y = estimate), data = dur_pop_gs_slope_preds, linewidth = 1) +
  
  # aesthetics
  theme_bw() +
  scale_y_continuous(expand = c(0,0),
                     limits = c(175,440),
                     breaks = seq(180,450,by=60),
                     labels = c("3","4","5","6","7")) +
  scale_x_continuous(
    breaks = seq(0,40,  by = 10)) +
  ylab("Dive duration (min)") +
  xlab("Slope (degrees)") +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, color = "black"),
    strip.text = element_text(size = 11)
  )
dupo_slope

## chla ## 
dur_pop_gs_chla_preds <- plot_predictions(dur_pop_gs, condition = "logchla", type = "response",
                                          draw = F)
## chla, zoomed in to the effect ## 
dupo_chla <- ggplot(dives_pop, aes(x = logchla, y = dur_secs)) +
  # observed points
  geom_rug(sides = "b") +
  
  # confidence ribbon
  geom_ribbon(aes(ymin=conf.low, ymax = conf.high,x=logchla), data = dur_pop_gs_chla_preds,
              alpha = 0.2, color = NA) +
  
  # predicted estimates
  geom_line(aes(y = estimate), data = dur_pop_gs_chla_preds, linewidth = 1) +
  
  # aesthetics
  theme_bw() +
  scale_y_continuous(expand = c(0,0),
                     limits = c(200,600),
                     breaks = seq(240,540,by=60),
                     labels = c("4","5","6","7","8","9")
                     ) +
  ylab("Dive duration (min)") +
  xlab(bquote('log10(Chlorophyll-a ('*mg^3*'))')) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, color = "black"),
    strip.text = element_text(size = 11)
  )
dupo_chla

## mld ## 
dur_pop_gs_mld_preds <- plot_predictions(dur_pop_gs, condition = "mld", type = "response",
                                          draw = F)
## mld, zoomed in to the effect ## 
dupo_mld <- ggplot(dives_pop, aes(x = mld, y = dur_secs)) +
  # observed points
  geom_rug(sides = "b") +
  
  # confidence ribbon
  geom_ribbon(aes(ymin=conf.low, ymax = conf.high,x=mld), data = dur_pop_gs_mld_preds,
              alpha = 0.2, color = NA) +
  
  # predicted estimates
  geom_line(aes(y = estimate), data = dur_pop_gs_mld_preds, linewidth = 1) +
  
  # aesthetics
  theme_bw() +
  scale_y_continuous(expand = c(0,0),
                     limits = c(175,440),
                     breaks = seq(180,450,by=60),
                     labels = c("3","4","5","6","7")) +
  ylab("Dive duration (mins)") +
  xlab("Mixed layer depth (m)") +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, color = "black"),
    strip.text = element_text(size = 11)
  )
dupo_mld

# arrange
dupo <- ggarrange(dupo_hms, dupo_id_line, dupo_moon, dupo_slope, dupo_chla, dupo_mld)
dupo

ggsave(here("outputs","gamm_plots","spatiotemporal_gamm_population_P09",
            "duration_gamm_spatiotemporal_population_wP09_gs_ff_all_predplot_2025Mar15.png"),
       width = 12, height = 6, units = "in")

## hourly dive rate ## -------------------------------------------------------- ##

# read in the original data that was formatted in the model script
rate_temp <- readRDS(here("pipeline","gamm_models","temporal_gamms_rate_model_data_2025Feb20.rds"))

# read in model object
rate_gs <- readRDS(here("pipeline","gamm_models","temporal_gamms_allpops_rate_gs_2025Feb20.rds"))

## hour of day ## 

# extract predictions from model for time of day variable (global smooth)
rate_gs_hod_preds <- plot_predictions(rate_gs, condition = "hod", type = "response",
                                      draw = F)

## time of day but zoom into the effect ## 
rate_hod <- ggplot(rate_temp, aes(x = hod, y = n_dives/hrs)) +
  # observed points
  geom_rug(sides = "b") +
  
  # confidence ribbon
  geom_ribbon(aes(ymin=conf.low, ymax = conf.high, x=hod), data = rate_gs_hod_preds,
              alpha = 0.2, color = NA) +
  
  # predicted estimates
  geom_line(aes(y = estimate), data = rate_gs_hod_preds, linewidth = 1) +
  
  # aesthetics
  theme_bw() +
  scale_y_continuous(limits = c(0,1),
                     breaks = seq(0,0.75, by = 0.25)) +
  scale_x_continuous(limits = c(0,23), expand = c(0,0)) +
  ylab("Dive rate (# dives/hr)") +
  xlab("Hour of day (HST)") +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, color = "black"),
    strip.text = element_text(size = 11)
  )
rate_hod

## time of day x tag ID ##

# extract predictions 
rate_gs_hod_id_preds <- plot_predictions(rate_gs, condition = c("hod","DeployID"), type = "response",
                                         draw = F)

# population color palette
id_pop_pal <- c(rep("#015b58",6),rep("#89689d",3),rep("#e69b99",3))

rate_hod_id <- ggplot(rate_temp, aes(x = hod, y = n_dives/hrs)) +
  # observed points
  geom_rug(sides = "b") +
  
  # predicted estimates
  geom_line(aes(y = estimate, color = DeployID), data = rate_gs_hod_id_preds, linewidth = 0.75,
            alpha = 0.75) +
  
  # aesthetics
  scale_color_manual(values = id_pop_pal) +
  theme_bw() +
  scale_y_continuous(limits = c(0,1.5),
                     breaks = seq(0,1.25,by = 0.25)) +
  scale_x_continuous(limits = c(0,23), expand = c(0,0)) +
  ylab("Dive rate (# dives/hr)") +
  xlab("Hour of day (HST)") +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, color = "black"),
    strip.text = element_text(size = 11)
  )
rate_hod_id


## moon phase ##

# extract predictions from model for moon phase (global smooth)
rate_gs_moon_preds <- plot_predictions(rate_gs, condition = "moon_phase_rads", type = "response",
                                       draw = F)

rate_moon <- ggplot(rate_temp, aes(x = moon_phase_rads, y = n_dives/hrs)) +
  # observed points
  geom_rug(sides = "b") +
  
  # confidence ribbon
  geom_ribbon(aes(ymin=conf.low, ymax = conf.high,x=moon_phase_rads), data = rate_gs_moon_preds,
              alpha = 0.2, color = NA) +
  
  # predicted estimates
  geom_line(aes(y = estimate), data = rate_gs_moon_preds, linewidth = 1) +
  
  # aesthetics
  theme_bw() +
  scale_y_continuous(
                     limits = c(0,1.5),
                     breaks = seq(0,1.25,by = 0.25)) +
  scale_x_continuous(
                     limits = c(0,6.3),
                     breaks = c(0,pi/2,pi,3*pi/2,6.3),
                     labels = c("New","First quarter","Full","Last quarter","New")) +
  ylab("Dive rate (# dives/hr)") +
  xlab("Moon phase") +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, color = "black"),
    strip.text = element_text(size = 11)
  )
rate_moon

# arrange
rate_all <- ggarrange(rate_hod, rate_hod_id, rate_moon, nrow = 1, ncol = 3)
rate_all


ggsave(here("outputs","gamm_plots","temporal_gamm_P09",
            "rate_nb_gamm_temporal_wP09_gs_all_predplot_2025Mar15.png"),
       width = 12, height = 3, units = "in")
