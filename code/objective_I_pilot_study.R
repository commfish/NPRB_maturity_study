#notes---
#author: Sara E Miller 
#contact: sara.miller@alaska.gov; 907-465-4245
#Last edited: February, 2019

# load ----
devtools::install_github("ben-williams/FNGr")
#font_import() #only do this one time - it takes a while
loadfonts(device="win")

is_installed <- function(mypkg){ is.element(mypkg, installed.packages()[,1])}

load_or_install <- function(package_names){
  for(package_name in package_names){
    if(!is_installed(package_name)){install.packages(package_name,repos="http://lib.stat.cmu.edu/R/CRAN")}
    library(package_name,character.only=TRUE,quietly=TRUE,verbose=FALSE)
  }}

load_or_install(c("extrafont",
                  "tidyverse",
                  "devtools",
                  "FNGr",
                  "cowplot"))
library(extrafont)
library(tidyverse)
library(devtools)
library(FNGr)
library(cowplot)

loadfonts(device="win")
windowsFonts(Times=windowsFont("TT Times New Roman"))
theme_set(theme_bw(base_size=12, base_family='Times New Roman') +
            theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  plot.title = element_text(family='Times New Roman', hjust = 0.5, size=12)))

# load data ---- 
data <- read.csv("data/prop_pilot_study.csv", check.names = FALSE) 

#clean data ----
data %>%
  filter(total >0)-> data  #delete samples without measurements 

#analysis----
#summarize by region and individual fish
data %>%
  mutate(prop1=(Focus_An1)/total,
         prop2=(Focus_An1+An1_An2)/total,
         prop3=(Focus_An1+An1_An2+An2_An3)/total,
         prop4=(Focus_An1+An1_An2+An2_An3+An3_An4)/total,
         prop5=(Focus_An1+An1_An2+An2_An3+An3_An4+An4_An5)/total,
         prop6=(Focus_An1+An1_An2+An2_An3+An3_An4+An4_An5+An5_An6)/total,
         prop7=(Focus_An1+An1_An2+An2_An3+An3_An4+An4_An5+An5_An6+An6_An7)/total,
         prop8=(Focus_An1+An1_An2+An2_An3+An3_An4+An4_An5+An5_An6+An6_An7+An7_An8)/total,
         prop9=(Focus_An1+An1_An2+An2_An3+An3_An4+An4_An5+An5_An6+An6_An7+An7_An8+An8_An9)/total) %>% 
    group_by(slide, zone) %>% 
    summarise(mean1=mean(prop1),
              mean2=mean(prop2),
              mean3=mean(prop3),
              mean4=mean(prop4),
              mean5=mean(prop5),
              mean6=mean(prop6),
              mean7=mean(prop7),
              mean8=mean(prop8),
              mean9=mean(prop9),
              stdev1=sd(prop1),
              stdev2=sd(prop2),
              stdev3=sd(prop3),
              stdev4=sd(prop4),
              stdev5=sd(prop5),
              stdev6=sd(prop6),
              stdev7=sd(prop7),
              stdev8=sd(prop8),
              stdev9=sd(prop9),
              n = n(),
              se1 = stdev1/sqrt(n()),
              se2 = stdev2/sqrt(n()),
              se3 = stdev3/sqrt(n()),
              se4 = stdev4/sqrt(n()),
              se5 = stdev5/sqrt(n()),
              se6 = stdev6/sqrt(n()),
              se7 = stdev7/sqrt(n()),
              se8 = stdev8/sqrt(n()),
              se9 = stdev9/sqrt(n())) -> data_fish

#data across fish (by region summaries)
data %>%
  mutate(prop1=(Focus_An1)/total,
         prop2=(Focus_An1+An1_An2)/total,
         prop3=(Focus_An1+An1_An2+An2_An3)/total,
         prop4=(Focus_An1+An1_An2+An2_An3+An3_An4)/total,
         prop5=(Focus_An1+An1_An2+An2_An3+An3_An4+An4_An5)/total,
         prop6=(Focus_An1+An1_An2+An2_An3+An3_An4+An4_An5+An5_An6)/total,
         prop7=(Focus_An1+An1_An2+An2_An3+An3_An4+An4_An5+An5_An6+An6_An7)/total,
         prop8=(Focus_An1+An1_An2+An2_An3+An3_An4+An4_An5+An5_An6+An6_An7+An7_An8)/total,
         prop9=(Focus_An1+An1_An2+An2_An3+An3_An4+An4_An5+An5_An6+An6_An7+An7_An8+An8_An9)/total) %>% 
  
  group_by(zone) %>% 
  summarise(mean1=mean(prop1,na.rm=T),
            mean2=mean(prop2,na.rm=T),
            mean3=mean(prop3,na.rm=T),
            mean4=mean(prop4,na.rm=T),
            mean5=mean(prop5,na.rm=T),
            mean6=mean(prop6,na.rm=T),
            mean7=mean(prop7,na.rm=T),
            mean8=mean(prop8,na.rm=T),
            mean9=mean(prop9,na.rm=T),
            stdev1=sd(prop1,na.rm=T),
            stdev2=sd(prop2,na.rm=T),
            stdev3=sd(prop3,na.rm=T),
            stdev4=sd(prop4,na.rm=T),
            stdev5=sd(prop5,na.rm=T),
            stdev6=sd(prop6,na.rm=T),
            stdev7=sd(prop7,na.rm=T),
            stdev8=sd(prop8,na.rm=T),
            stdev9=sd(prop9,na.rm=T),
            n = n(),
            se1 = stdev1/sqrt(n()),
            se2 = stdev2/sqrt(n()),
            se3 = stdev3/sqrt(n()),
            se4 = stdev4/sqrt(n()),
            se5 = stdev5/sqrt(n()),
            se6 = stdev6/sqrt(n()),
            se7 = stdev7/sqrt(n()),
            se8 = stdev8/sqrt(n()),
            se9 = stdev9/sqrt(n())) -> data_region