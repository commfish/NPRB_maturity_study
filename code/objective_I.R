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

#data clean ----
data %>% 
  filter(!(maturation_status_histology %in% c("", NA, "no slide", "no score", "4-1", "3-4", "2-3", "1-2")))-> data_clean

# check for duplicate image names in the increment file
increment %>% 
  group_by(image_name) %>% 
  filter(n()>1) %>% 
  summarize(n=n())-> x

# match increments to awl data
merge1<-merge(increment, data_clean, all.y=T)
merge1 %>% 
  filter(!(Increment1 %in% c("", NA)))-> clean_dataset #n=573 useable samples


