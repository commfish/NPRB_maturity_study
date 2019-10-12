#notes---
#author: Sara E Miller 
#contact: sara.miller@alaska.gov; 907-465-4245
#Last edited: July, 2019
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4885900/
#https://stats.stackexchange.com/questions/82105/mcfaddens-pseudo-r2-interpretation

# load libraries----
# devtools::install_github("ben-williams/FNGr")
source("code/functions.R")
library(extrafont)
loadfonts(device="win")
library(nlme)
library(tidyverse)
library(FNGr)
library(broom)
library(cowplot)
library(mixtools)
library(car)
library(modEvA) #pseudo-R squared
library(MASS)
library(digest)
#library(lmtest)
#library(mgcv)
#library(visreg)
#library(boot)
#library(AICcmodavg) #AICc
#library(rms)
#library(mixdist)
#library(fitdistrplus)
#library(gam)
theme_set(theme_sleek())

#Measurements and maturity II and above are mature
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
  filter(!(Increment1 %in% c("", NA))) %>% 
  filter(!(scale_region%in% c("OOA (OUT OF AREA)", "A", "C", "D", "E", "G", "H")))-> clean_dataset #n=296 useable samples

#calculate proportion of outer scale ring
clean_dataset %>% 
  mutate(mature = ifelse(maturity_state_histology == 'I', 'immature', 'mature'),
         anu1 = ifelse(is.na(Increment1), 0 , Increment1),
         anu2 = ifelse(is.na(Increment2), 0 , Increment2),
         anu3 = ifelse(is.na(Increment3), 0 , Increment3),
         anu4 = ifelse(is.na(Increment4), 0 , Increment4),
         anu5 = ifelse(is.na(Increment5), 0 , Increment5),
         anu6 = ifelse(is.na(Increment6), 0 , Increment6),
         anu7 = ifelse(is.na(Increment7), 0 , Increment7)) %>%
mutate(annulus = ifelse(anu7 >0, anu7,
                        ifelse(anu6 >0, anu6,
                               ifelse(anu5 >0, anu5,
                                      ifelse(anu4 >0, anu4,
                                             ifelse(anu3 >0, anu3,
                                                    ifelse(anu2 >0, anu2, anu1))))))) %>%
  mutate(age=as.factor(age),
         maturity = ifelse(maturation_status_histology == 1, 0, 1),
         mature = ifelse(maturation_status_histology == 1, "immature", "mature")) %>%
  dplyr::select(annulus, mature, maturity, maturation_status_histology, age, sex_histology) -> merge
  write.csv(merge, file= "data/ouput.csv") 

# Boxplot Figures---- 
merge %>% 
  ggplot(data=., aes(x = age, y = annulus, color = mature)) +
  geom_jitter(size = 1) +
  labs(x = "Age",y = "outer increment proportion") + 
  theme(legend.title=element_blank(), legend.position=c(.8,.9)) +
  scale_color_manual(values=c("#999999", "black")) +
  scale_fill_manual(values=c("#999999", "black" )) +
  scale_y_continuous(breaks = c(0, 1,2,3), limits = c(0, 3)) +
  geom_text(aes(x = 1, y = 0.6, label="a)", hjust = 1),family="Times New Roman", colour="black", size=5)-> plot1


merge %>%
  filter(age == 3) %>% 
  ggplot(data=.,aes(x = mature, y = annulus, color=mature)) + labs(x = "",
  y = "Outer increment proportion (age-3)")  +
  geom_boxplot() + theme(legend.position= "none") +
  scale_color_manual(values=c("#999999", "black")) +
  scale_fill_manual(values=c("#999999", "black" )) +
  scale_y_continuous(breaks = c(0, 1,2,3), limits = c(0, 3)) +
  geom_text(aes(x = 1.2, y = 0.6, label="b)", hjust = 4),family="Times New Roman", colour="black", size=5)-> plot2

merge %>%
  filter(age == 3) %>%
  filter(sex_histology == 'Female') %>%
  ggplot(data=.,aes(x = mature, y = annulus, color=mature)) + labs(x = "",
                                                                      y = "Outer increment proportion (females age-3)")  +
  geom_boxplot() +
  scale_color_manual(values=c("#999999", "black")) + theme(legend.position= "none") +
  
  scale_fill_manual(values=c("#999999", "black" )) + 
  scale_y_continuous(breaks = c(0, 1,2,3), limits = c(0, 3)) +
  geom_text(aes(x = 1.2, y = 0.6, label="c)", hjust = 4),family="Times New Roman", colour="black", size=5)-> plot3

cowplot::plot_grid(plot1, plot2, plot3, align = "vh", nrow = 1, ncol=3)
ggsave("figs/boxplot_measurements_II_above_mature.png", dpi = 500, height = 4, width = 8, units = "in")


#Generalized Linear models (age-3 only)----
merge %>%
  filter(age == 3) %>%
  filter(sex_histology == 'Female') -> age3

merge %>%
  filter(age == 3) -> age3

fit <- glm(maturity ~ (annulus) , family = binomial, data = age3) 
Anova(fit)
RsqGLM(fit)#peudo R2 
summary(fit)

#Measurements and maturity III and above are mature
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
  filter(!(Increment1 %in% c("", NA))) %>% 
  filter(!(scale_region%in% c("OOA (OUT OF AREA)", "A", "C", "D", "E", "G", "H")))-> clean_dataset #n=296 useable samples

#calculate proportion of outer scale ring
clean_dataset %>% 
  mutate(mature = ifelse(maturity_state_histology == 'I', 'immature', 'mature'),
         anu1 = ifelse(is.na(Increment1), 0 , Increment1),
         anu2 = ifelse(is.na(Increment2), 0 , Increment2),
         anu3 = ifelse(is.na(Increment3), 0 , Increment3),
         anu4 = ifelse(is.na(Increment4), 0 , Increment4),
         anu5 = ifelse(is.na(Increment5), 0 , Increment5),
         anu6 = ifelse(is.na(Increment6), 0 , Increment6),
         anu7 = ifelse(is.na(Increment7), 0 , Increment7)) %>%
  mutate(annulus = ifelse(anu7 >0, anu7,
                          ifelse(anu6 >0, anu6,
                                 ifelse(anu5 >0, anu5,
                                        ifelse(anu4 >0, anu4,
                                               ifelse(anu3 >0, anu3,
                                                      ifelse(anu2 >0, anu2, anu1))))))) %>%
  mutate(age=as.factor(age),
         maturity = ifelse(maturation_status_histology == 1, 0,
                           ifelse(maturation_status_histology == 2, 0, 1))) %>%
  mutate(mature = ifelse(maturity == 1, "mature", "immature")) %>%
  dplyr::select(annulus, mature, maturity, maturation_status_histology, age, sex_histology) -> merge
write.csv(merge, file= "data/ouput.csv") 

# Boxplot Figures---- 
merge %>% 
  ggplot(data=., aes(x = age, y = annulus, color = mature)) +
  geom_jitter(size = 1) +
  labs(x = "Age",y = "outer increment proportion") + 
  theme(legend.title=element_blank(), legend.position=c(.8,.9)) +
  scale_color_manual(values=c("#999999", "black")) +
  scale_fill_manual(values=c("#999999", "black" )) +
  scale_y_continuous(breaks = c(0, 1,2,3), limits = c(0, 3)) +
  geom_text(aes(x = 1, y = 0.6, label="a)", hjust = 1),family="Times New Roman", colour="black", size=5)-> plot1


merge %>%
  filter(age == 3) %>% 
  ggplot(data=.,aes(x = mature, y = annulus, color=mature)) + labs(x = "",
                                                                   y = "Outer increment proportion (age-3)")  +
  geom_boxplot() + theme(legend.position= "none") +
  scale_color_manual(values=c("#999999", "black")) +
  scale_fill_manual(values=c("#999999", "black" )) +
  scale_y_continuous(breaks = c(0, 1,2,3), limits = c(0, 3)) +
  geom_text(aes(x = 1.2, y = 0.6, label="b)", hjust = 4),family="Times New Roman", colour="black", size=5)-> plot2

merge %>%
  filter(age == 3) %>%
  filter(sex_histology == 'Female') %>%
  ggplot(data=.,aes(x = mature, y = annulus, color=mature)) + labs(x = "",
                                                                   y = "Outer increment proportion (females age-3)")  +
  geom_boxplot() +
  scale_color_manual(values=c("#999999", "black")) + theme(legend.position= "none") +
  
  scale_fill_manual(values=c("#999999", "black" )) + 
  scale_y_continuous(breaks = c(0, 1,2,3), limits = c(0, 3)) +
  geom_text(aes(x = 1.2, y = 0.6, label="c)", hjust = 4),family="Times New Roman", colour="black", size=5)-> plot3

cowplot::plot_grid(plot1, plot2, plot3, align = "vh", nrow = 1, ncol=3)
ggsave("figs/boxplot_measurements_III_above_mature.png", dpi = 500, height = 4, width = 8, units = "in")


#Generalized Linear models (age-3 only)----
merge %>%
  filter(age == 3) %>%
  filter(sex_histology == 'Female') -> age3

merge %>%
  filter(age == 3) -> age3

fit <- glm(maturity ~ (annulus) , family = binomial, data = age3) 
Anova(fit)
RsqGLM(fit)#peudo R2 
summary(fit)
