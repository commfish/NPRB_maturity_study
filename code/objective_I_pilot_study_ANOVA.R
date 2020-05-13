#notes---
#author: Sara E Miller 
#contact: sara.miller@alaska.gov; 907-465-4245
#Last edited: February, 2019

# load libraries----
# devtools::install_github("commfish/fngr")
# devtools::install_github('droglenc/RFishBC')
# devtools::install_github('droglenc/FSA')

library(extrafont)
windowsFonts(Times=windowsFont("Times New Roman"))
library(tidyverse)
library(fngr)
library(broom)
library(gridExtra)
library(cowplot)
library(FSA)
library(RFishBC)
library(magrittr)
library(ggpmisc)

theme_sleek <- function(base_size = 12, base_family = "Arial") {
  half_line <- base_size/2
  theme_light(base_size = 12, base_family = "Arial") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks.length = unit(half_line / 2.2, "pt"),
      strip.background = element_rect(fill = NA, colour = NA),
      strip.text.x = element_text(colour = "black"),
      strip.text.y = element_text(colour = "black"),
      panel.border = element_rect(fill = NA),
      legend.key.size = unit(0.9, "lines"),
      legend.key = element_rect(colour = NA, fill = NA),
      legend.background = element_rect(colour = NA, fill = NA)
    )
}


theme_set(theme_sleek())

asintrans <- function(p) {asin(sqrt(p))} # arctransformation function

# load data ---- 
data <- read.csv("data/prop_pilot_study.csv", check.names = FALSE) %>%
  mutate(rad1 = anu1,#calculate radial units (ie. focus to each annulus)
         rad2 = anu1 + anu2,
         rad3 = anu1+anu2+anu3,
         rad4 = anu1+anu2+anu3+anu4,
         rad5 = anu1+anu2+anu3+anu4+anu5,
         rad6 = anu1+anu2+anu3+anu4+anu5+anu6,
         rad7 = anu1+anu2+anu3+anu4+anu5+anu6+anu7,
         rad8 = anu1+anu2+anu3+anu4+anu5+anu6+anu7+anu8,
         rad9 = anu1+anu2+anu3+anu4+anu5+anu6+anu7+anu8+anu9,
         prop1 = rad1 / radcap,
         prop2 = rad2 / radcap,
         prop3 = rad3 / radcap,
         prop4 = rad4 / radcap,
         prop5 = rad5 / radcap,
         prop6 = rad6 / radcap,
         prop7 = rad7 / radcap,
         prop8 = rad8 / radcap,
         prop9 = rad9 / radcap) -> data 

#analysis----
#summarize by region and individual fish
data %>%
  group_by(zone) %>% 
  mutate(n = n()) %>% 
  group_by(zone, n) %>% 
  summarise_at(vars(contains('prop')), funs(mean, sd), na.rm = T) %>% 
  mutate( se1 = prop1_sd / sqrt(n()),
          se2 = prop2_sd / sqrt(n()),
          se3 = prop3_sd / sqrt(n()),
          se4 = prop4_sd / sqrt(n()),
          se5 = prop5_sd / sqrt(n()),
          se6 = prop6_sd / sqrt(n()),
          se7 = prop7_sd / sqrt(n()),
          se8 = prop8_sd / sqrt(n()),
          se9 = prop9_sd / sqrt(n())) -> data_region

# Randomly sample two scales from each region
set.seed(167) 

data %>% 
  dplyr::select(id, fish, agecap, lencap, contains('anu'), 
                contains('rad'), zone) %>% 
  group_by(fish, zone) %>% 
  sample_n(2) -> sample1   # sample two scales from each fish zone (5 scales sampled from each fish)

set.seed(167) 
data %>% 
  dplyr::select(id, fish, agecap, lencap, contains('anu'), 
                contains('rad'), zone) %>% 
  group_by(fish, zone) %>% 
  sample_n(2) %>%  # sample two scales from each fish zone (5 scales sampled from each fish)
  gather(agei, radi, rad1:rad9) %>%
  arrange(id, agei) %>% 
  dplyr::select(id, fish, agecap, lencap, radcap, agei, radi, zone) %>% 
  mutate(agei = as.numeric(stringr::str_sub(agei, 4, 4))) %>% 
  filter(!is.na(radi), agei<=agecap) -> sample2

# ANOVA
sample1 %>% 
  filter(agecap == 4) %>%
  group_by(fish)  %>%
  filter(radcap > 0) %>%
  filter(!fish %in% c(28, 47,77,205)) %>%
 summarize (n=n()) -> table1 # sample size for age-4
 
sample1 %>% 
  filter(agecap == 4) %>%
  group_by(fish)  %>%
  filter(radcap > 0) %>%
  filter(!fish %in% c(28, 47,77,205)) %>%
  mutate(prop_age2 =rad2/radcap) -> sample_data

# Model 1
delivery.mod1 = aov(radcap ~ zone * fish, data = sample_data)
summary(delivery.mod1)
delivery.res = sample_data
delivery.res$M1.fit = fitted(delivery.mod1)
delivery.res$M1.Resid = resid(delivery.mod1)

ggplot(delivery.res, aes(M1.fit, M1.resid, colour = fish)) + geom_point() +
  xlab("Fitted values") + ylab("Residuals")

ggplot(delivery.res, aes(M1.fit, M1.resid, colour = fish)) + geom_point() +
  xlab("Fitted values") + ylab("Residuals") + facet_wrap (~zone)

ggplot(delivery.res, aes(sample = M1.resid)) + stat_qq()

TukeyHSD(delivery.mod1, which = "zone")

# Model 2
delivery.mod2 = aov(prop_age2 ~ zone * fish, data = sample_data)
summary(delivery.mod2)
delivery.res = sample_data
delivery.res$M2.fit = fitted(delivery.mod2)
delivery.res$M2.Resid = resid(delivery.mod2)

ggplot(delivery.res, aes(M2.fit, M2.resid, colour = fish)) + geom_point() +
  xlab("Fitted values") + ylab("Residuals")

ggplot(delivery.res, aes(M2.fit, M2.resid, colour = fish)) + geom_point() +
  xlab("Fitted values") + ylab("Residuals") + facet_wrap (~zone)

ggplot(delivery.res, aes(sample = M2.resid)) + stat_qq()

TukeyHSD(delivery.mod2, which = "zone")
