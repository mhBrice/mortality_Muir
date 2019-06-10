# Mortality analyses in an old-growth forest

This repository includes the data and R scripts to reproduce analyses and figures found in the article *Long-term impact of a major ice storm on tree mortality in an old-growth forest* by Deschênes, Brice and Brisson published in Forest Ecology and Management.

All data used for the analyses can be found in the [data](https://github.com/mhBrice/mortality_Muir/tree/master/data) folder. Script #1 prepares and formats the data for the analyses performed in the script #2. The figures produced are in the [results](https://github.com/mhBrice/mortality_Muir/tree/master/results) folder.

## Packages required to run the scripts

library(survival)
library(survminer)

library(dplyr)

library(sf)
library(units)

library(gtools)
library(sjPlot)

library(RColorBrewer)
library(colourlovers)

library(graphicsutils)
