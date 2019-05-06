#notes---
#author: Sara E Miller 
#contact: sara.miller@alaska.gov; 907-465-4245
#Last edited: December, 2018

# load ----
#font_import() #only do this one time - it takes a while
loadfonts(device="win")

fun_length <- function(x){
  return(data.frame(y=median(x),label= paste0("n=", length(x))))
}

library(extrafont)
library(tidyverse)
library(dgof)
library(Matching)
library(devtools)
library(FNGr)
library(cowplot)
library(psych)
devtools::install_github("ben-williams/FNGr")


# load data ---- 
data <- read.csv("data/data_11_26_2018.csv", check.names = FALSE) 

#clean data ----

data %>% filter(!(maturation_status_histology %in% c("", NA, "no slide", "no score", "4-1", "3-4", "2-3", "1-2")))-> data_clean

data_clean %>%
  group_by(maturity_state_field, maturation_status_histology) %>%
  summarise(count = n()) -> table3 

#Macroscopic versus Histology
#cohen kappa evaluation 
data %>% filter(!(maturation_status_histology %in% c("", NA, "no slide", "no score", "4-1", "3-4", "2-3", "1-2", "7", "6"))) %>%
  dplyr::select(maturation_status_histology = maturation_status_histology,
                maturity_state_field = maturity_state_field) -> data_clean
data_clean %>%
  group_by(maturity_state_field, maturation_status_histology) %>%
  summarise(count = n()) -> table3a 
cohen.kappa(data_clean) # all stages

#Two proportions z-test with continuity correction
#The function returns:
  #the value of Pearson’s chi-squared test statistic.
  #a p-value
  #a 95% confidence intervals
  #an estimated probability of success 

#Pearson's X^2
data %>% filter(!(maturation_status_histology %in% c("", NA, "no slide", "no score", "4-1", "3-4", "2-3", "1-2", "7", "6"))) %>%
  mutate(maturation_status_histology = as.numeric(maturation_status_histology),
                      maturity_state_field = as.numeric(maturity_state_field)) %>%
  mutate(status_h = ifelse(maturation_status_histology>2, 'mature','immature'),
         status_m = ifelse(maturity_state_field>2, 'mature','immature')) -> data_clean

data_clean %>%
  dplyr::select(status_h) %>%
  group_by(status_h) %>%
  summarise(count = n()) -> tablex  

data_clean %>%
  dplyr::select(status_m) %>%
  group_by(status_m) %>%
  summarise(count = n()) -> tabley

res <- prop.test(x = c(337, 335), n = c(726, 726),  correct = FALSE) #immature

data %>% filter(!(maturation_status_histology %in% c("", NA, "no slide", "no score", "4-1", "3-4", "2-3", "1-2", "7", "6"))) %>%
  mutate(maturation_status_histology = as.numeric(maturation_status_histology),
         maturity_state_field = as.numeric(maturity_state_field)) %>%
  mutate(status_h = ifelse(maturation_status_histology>2, 'mature','immature'),
         status_m = ifelse(maturity_state_field>2, 'mature','immature')) %>%
dplyr::select(status_h, status_m) -> data_clean

cohen.kappa(data_clean) #this is higher since only by immature/mature

#GSI versus histology
data %>% filter(!(maturation_status_histology %in% c("", NA, "no slide", "no score", "4-1", "3-4", "2-3", "1-2", "7", "6"))) %>%
  mutate(maturation_status_histology = as.numeric(maturation_status_histology),
                      GSI = as.numeric(GSI)) %>%
  mutate(status_h = ifelse(maturation_status_histology>2, 'mature','immature'),
         status_gsi = ifelse(GSI>0.0499, 'mature','immature')) -> data_clean

data_clean %>%
  dplyr::select(status_h) %>%
  group_by(status_h) %>%
  summarise(count = n()) -> tablex  

data_clean %>%
  dplyr::select(status_gsi) %>%
  group_by(status_gsi) %>%
  summarise(count = n()) -> tabley

res <- prop.test(x = c(337, 709), n = c(726, 726),  correct = FALSE)

data %>% filter(!(maturation_status_histology %in% c("", NA, "no slide", "no score", "4-1", "3-4", "2-3", "1-2", "7", "6"))) %>%
  mutate(maturation_status_histology = as.numeric(maturation_status_histology),
         GSI = as.numeric(GSI)) %>%
  mutate(status_h = ifelse(maturation_status_histology>2, 'mature','immature'),
         status_gsi = ifelse(GSI>0.0499, 'mature','immature'))%>%
  dplyr::select(status_h, status_gsi) -> data_clean

cohen.kappa(data_clean)

#GSI table (Table 5)
data %>% filter(!(maturation_status_histology %in% c("", NA, "no slide", "no score", "4-1", "3-4", "2-3", "1-2", "7", "6"))) %>%
 filter(!(GSI %in% c(0))) %>%
  mutate(maturation_status_histology = as.numeric(maturation_status_histology),
         GSI = as.numeric(GSI)) %>%
  mutate(status_h = ifelse(maturation_status_histology>2, 'mature','immature')) -> data_clean

data_clean %>%
  dplyr::select(status_h, GSI, sex_histology, catch_date) %>%
  group_by(status_h,sex_histology, catch_date ) %>%
  summarise_all(funs(n(),median, mean,sd,se=sd(.)/sqrt(n())))-> tablex  

data_clean %>%
  dplyr::select(status_h, GSI, sex_histology) %>%
  group_by(status_h,sex_histology) %>%
  summarise_all(funs(n(),median, mean,sd,se=sd(.)/sqrt(n())))-> tabley  

#Figures 
loadfonts(device="win")
windowsFonts(Times=windowsFont("TT Times New Roman"))
theme_set(theme_bw(base_size=12, base_family='Times New Roman') +
            theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  plot.title = element_text(family='Times New Roman', hjust = 0.5, size=12),
                  axis.text.y = element_text(size=12,colour="black",family="Times New Roman"),
                  axis.title.y = element_text(size=12,colour="black",family="Times New Roman")))
tickr_length <- data.frame(length_millimeters = 200:1000)
axisf <- tickr(tickr_length, length_millimeters, 100)

#bubble plot of age versus histology
bb <- c(201,100,50,20,10,5) # define breaks.
ll <- c("200+","100","50","20","10", "5") # labels.
p<-ggplot(table1)+ geom_point(aes(x = age, y = maturation_status_histology, size = count),shape=16, alpha=0.80) +
  labs(x="Age", y = "Maturation Stage (histology)") +  scale_size_continuous(name = "",
                                                                             breaks = bb,
                                                                             limits = c(0, 250),
                                                                             labels = ll,
                                                                             range = c(0, 10))
ggsave(filename = 'figs/bubbleplot.png', dpi =200, width=6, height=8, units = "in")

#Box plot of GSI cutoff
data %>% filter(!(maturation_status_histology %in% c("", NA, "no slide", "no score", "4-1", "3-4", "2-3", "1-2", "7", "6"))) %>%
  filter(!(GSI %in% c(0))) %>%
  mutate(maturation_status_histology = as.numeric(maturation_status_histology),
         GSI = as.numeric(GSI)) %>%
  mutate(status_h = ifelse(maturation_status_histology>1, 'Mature','Immature')) -> data_clean

data_clean %>% filter(sex_histology == "Female") -> female 
female<- as.data.frame(female)
data_clean %>% filter(sex_histology == "Male") -> male 

female<-ggplot(female)+ geom_boxplot(aes(x = status_h, y = GSI))+
  labs(x="", y = "GSI") + theme(legend.position="none")+
  geom_text(aes(x =10/25/2017, y= 0.075, label="Females"),hjust = -0.6,  family="Times New Roman", colour="black")

male<-ggplot(male)+ geom_boxplot(aes(x = status_h, y = GSI))+
  labs(x="", y = "GSI") + theme(legend.position="none")+
  geom_text(aes(x =10/25/2017, y= 0.025, label="Males"),hjust = -0.6,  family="Times New Roman", colour="black")
cowplot::plot_grid(female, male, nrow = 1, ncol=2) 
ggsave(filename = 'figs/gsi.png', dpi =200, width=8, height=6, units = "in")

#Boxplot of GSI cutoff by catch_date and sex
data %>% filter(!(maturation_status_histology %in% c("", NA, "no slide", "no score", "4-1", "3-4", "2-3", "1-2", "7", "6"))) %>%
  filter(!(GSI %in% c(0))) %>%
  mutate(maturation_status_histology = as.numeric(maturation_status_histology),
         GSI = as.numeric(GSI)) %>%
  mutate(status_h = ifelse(maturation_status_histology>1, 'Mature','Immature')) -> data_clean

data_clean %>% filter(sex_histology == "Female" )->female
female<- as.data.frame(female)

females<-ggplot(female, aes(status_h, GSI, colour = status_h)) +
  geom_boxplot(position=position_dodge(width=0.9)) +
  stat_summary( fun.data = fun_length, geom = "text", 
                position=position_dodge(width=0.9), vjust=2, family="Times New Roman", size=3)+
  labs(x="", y = "GSI") + theme(legend.position="none")+scale_colour_grey(start = .6, end = .1)+
  geom_text(aes(x =10/25/2017, y= 0.075, label="Females"),hjust = -0.6,  family="Times New Roman", colour="black")

data_clean %>% filter(sex_histology == "Male") -> male 
library(plyr)
ddply(male, .(catch_date), summarise, label = length(GSI))
      
males<-ggplot(male, aes(status_h, GSI, colour = status_h)) +
  geom_boxplot(position=position_dodge(width=0.9)) +
  stat_summary( fun.data = fun_length, geom = "text", 
                position=position_dodge(width=0.9), vjust=2, family="Times New Roman", size=3)+
  labs(x="", y = "GSI") + theme(legend.position="none")+scale_colour_grey(start = .6, end = .1)+
  geom_text(aes(x =10/25/2017, y= 0.075, label="Males"),hjust = -0.6,  family="Times New Roman", colour="black")
cowplot::plot_grid(females, males, nrow = 1, ncol=2)
ggsave(filename = 'figs/gsi_boxplots.png', dpi =200, width=8, height=6, units = "in")
ggplot_build(females)$data #boxplot data
ggplot_build(males)$data #boxplot data

#GSI sample size by catch date, sex, and histology stage
data %>% filter(!(maturation_status_histology %in% c("", NA, "no slide", "no score", "4-1", "3-4", "2-3", "1-2", "7", "6"))) %>%
  filter(!(GSI %in% c(0))) %>%
  mutate(maturation_status_histology = as.numeric(maturation_status_histology),
         GSI = as.numeric(GSI)) %>%
  mutate(status_h = ifelse(maturation_status_histology>1, 'Mature','Immature')) -> data_clean

data_clean %>%
  dplyr::select(maturation_status_histology , GSI, sex_histology, catch_date) %>%
  mutate(maturation_status_histology = ifelse(maturation_status_histology == 1, 'I',
                                              ifelse(maturation_status_histology == 2, 'II',
                                                     ifelse(maturation_status_histology == 3, 'III',
                                                            ifelse(maturation_status_histology == 4, 'IV', 'V'))))) %>%
  group_by(maturation_status_histology ,sex_histology, catch_date ) %>%
  summarise_all(funs(n()))-> tablez  
tablez  %>% filter(sex_histology == "Female") -> female 

females<-ggplot(female, aes(maturation_status_histology, GSI)) +
  geom_point(size=4, aes(shape = catch_date))+scale_shape_manual(values=c(16,6,3))+
  labs(x="Female maturity stage (based on histology)", y = "Sample size") + theme(legend.title=element_blank(), legend.position=c(.9,.9))

ggsave(filename = 'figs/gsi_sample_size_females.png', dpi =200, width=8, height=6, units = "in")
