# False Killer Whale Diving Behavior Analysis
This repository contains the code for processing, analyzing, and visualizing data presented in *Ecological and social contexts of diving behavior reveal foraging flexibility in Hawaiian false killer whales* Kratofil MA, Shaff JF, Hoffbauer HK, Cantor M, Hill MC, Baird RW. (in prep; to be submitted to *Movement Ecology*). Raw location and behavior data are being used for several current analyses and thus are not publicly available at this time. However, location data can be visualized on [Movebank](https://www.movebank.org/cms/movebank-main) and made available from the corresponding author upon reasonable request. 

Code can be cited as: Kratofil MA, Shaff JF, Hoffbauer HK (2025) Code for *Ecological and social contexts of diving behavior reveal foraging flexibility in Hawaiian false killer whales*. Zenodo. https:://doi.org/10.5281/zenodo.15306564 

Corresponding author: Michaela A. Kratofil (michaela.kratofil@oregonstate.edu or mkratofil@cascadiaresearch.org)

## Data processing
This folder contains R scripts and helper functions for reviewing (QA/QC) and processing dive and location data for further analyses. The movement data processing (i.e., raw Argos/GPS location data filtering) is not provided here; see [Kratofil et al. (2023)](https://www.frontiersin.org/journals/marine-science/articles/10.3389/fmars.2023.1053581/full) for details on these methods. Seafloor depth and oceanographic files are not provided due to space limitations, but can be found at the links provided in the manuscript. 

The behavior log QA/QC protocol follows steps reported in [Baird et al. (2019)](https://www.cascadiaresearch.org/files/publications/Bairdetal2019_Hatteras.pdf), [Cioffi et al. (2023)](https://animalbiotelemetry.biomedcentral.com/articles/10.1186/s40317-023-00334-1), and [Shearer et al. (2019)](https://royalsocietypublishing.org/doi/10.1098/rsos.181728), in addition to new steps reported here. The annotated markdown document is provided in this folder, and a brief summary:

* Depth and ZeroDepthOffset values of CRC-checked messages in each tag's -Status.csv file were examined for evidence of pressure transducer failures. Both values should be close to zero, and values greater than +/- 10 were considered as indicators of possible pressure transducer failures.
* Behavior log data occurring after the last CRC-checked status message were removed to ensure behavior log records used in analyses were free of transducer failures
* Message overlaps (> 60 seconds) were examined as these culd indicate corrupt messages that were erroneously retained in the behavior log
* Dive depth, duration, and mean ascent/descent rate (calculated as two times dive depth divided by dive duration) were examined for extreme values.
* Some tags also transmitted time series data, where the depth of the animal was recorded at user-defined intervals (see Additional File 1 in the paper), and these were compared to the behavior log data to additional assess for any indication of tag malfunctioning.

## Analyses
This folder contains R scripts for completing the analyses (e.g., summary statistics, multivariate models) reported in the manuscript. 

## Visualizations
This folder contains R scripts and helper functions for creating the figures that are presented in the manuscript and additional files. Note that for several multi-paneled figures, individual plots were saved in the R scripts and panels were subsequently arranged in Adobe Illustrator. 
