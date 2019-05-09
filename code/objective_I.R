#notes---
#author: Sara E Miller 
#contact: sara.miller@alaska.gov; 907-465-4245
#Last edited: May, 2019

# load libraries----
# devtools::install_github("ben-williams/FNGr")
library(extrafont)
loadfonts(device="win")
library(tidyverse)
library(FNGr)
library(broom)
library(cowplot)

theme_set(theme_sleek())

asintrans <- function(p) {asin(sqrt(p))} # arctransformation function

# load data ---- 
data <- read.csv("data/data_11_26_2018.csv", check.names = FALSE) 
increment <- read.csv("data/increment_measurements.csv", check.names = FALSE) 
