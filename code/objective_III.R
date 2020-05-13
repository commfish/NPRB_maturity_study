#author: Sara E Miller 
#contact: sara.miller@alaska.gov; 907-465-4245
#Last edited: May, 2020

# load ----
#font_import() #only do this one time - it takes a while
loadfonts(device="win")

source("code/functions.R")
devtools::install_github("commfish/fngr")
library(extrafont)
library(tidyverse)
library(dgof)
library(Matching)
library(devtools)
library(fngr)
library(cowplot)
library(psych)
library(broom)
library(lubridate)
library(plyr)
theme_set(theme_sleek())

# load data ---- 
data <- read.csv("data/data_11_26_2018.csv", check.names = FALSE) 
spawn <- read.csv("data/spawn_timing_Sitka.csv", check.names = FALSE) 

#clean data ----
data %>% filter(!(maturation_status_histology %in% c("", NA, "no slide", "no score", "4-1", "3-4", "2-3", "1-2"))) %>%
         filter(!(sex_histology %in% c("Male"))) %>%
         filter(!(sample_no %in% c("743", "629", "334", "289"))) -> data_clean
write.csv(data_clean, "data/cpue_new_test.csv") #dataset for GSI/histology by age

# Macroscopic versus Histology (all stages)
# cohen kappa evaluation 
data %>% 
  filter(!(maturation_status_histology %in% c("", NA, "no slide", "no score", "4-1", "3-4", "2-3", "1-2", "7", "6"))) %>%
  filter(!(sex_histology %in% c("Male"))) %>%
  filter(!(sample_no %in% c("743", "629", "334", "289"))) %>%
  dplyr::mutate(maturity_histology = ifelse(maturation_status_histology == 1, "immature" , "mature")) -> data_clean
data_clean %>%
  dplyr::select(maturation_status_histology,maturity_state_field) %>%
  as.data.frame() %>%
  dplyr::group_by(maturity_state_field, maturation_status_histology) %>%
  dplyr::summarise(count = n()) -> table1
data_clean %>%
  dplyr::select(maturation_status_histology,maturity_state_field) %>%
  cohen.kappa(.) # all stages (psych package)

# Histology by age
data %>% filter(!(maturation_status_histology %in% c("", NA, "no slide", "no score", "4-1", "3-4", "2-3", "1-2", "7", "6"))) %>%
  filter(!(sex_histology %in% c("Male"))) %>%
  filter(!(sample_no %in% c("743", "629", "334", "289"))) %>%
  dplyr::mutate(maturity_histology = ifelse(maturation_status_histology == 1, "immature" , "mature")) -> data_clean
data_clean %>%
  dplyr::select(maturity_histology,age) %>%
  group_by(maturity_histology, age) %>%
  summarise(count = n()) -> table2

# Macroscopic versus Histology (mature/immature)
# cohen kappa evaluation 
data %>% 
  filter(!(maturation_status_histology %in% c("", NA, "no slide", "no score", "4-1", "3-4", "2-3", "1-2", "7", "6"))) %>%
  filter(!(sex_histology %in% c("Male"))) %>%
  filter(!(sample_no %in% c("743", "629", "334", "289"))) %>%
  dplyr::mutate(maturity_histology = ifelse(maturation_status_histology == 1, "immature" , "mature")) %>%
  dplyr::mutate(maturity_field = ifelse(maturity_state_field == 1, "immature" , "mature")) %>%
  dplyr::select(maturity_histology, maturity_field) -> data_clean

data_clean %>%
  group_by(maturity_histology, maturity_field) %>%
  summarise(count = n()) -> table3 
cohen.kappa(data_clean) # all stages

# GSI versus histology
data %>% filter(!(maturation_status_histology %in% c("", NA, "no slide", "no score", "4-1", "3-4", "2-3", "1-2", "7", "6"))) %>%
  filter(!(sex_histology %in% c("Male"))) %>%
  filter(!(sample_no %in% c("743", "629", "334", "289"))) %>%
  dplyr::mutate(maturity_histology = ifelse(maturation_status_histology == 1, "immature" , "mature")) %>%
  mutate(GSI=as.numeric(GSI)) %>%
  dplyr::mutate(maturity_GSI = ifelse(GSI > 0.05, "mature" , "immature")) %>%
  dplyr::select(maturity_histology, maturity_GSI)  -> data_clean
cohen.kappa(data_clean) # all stages
data_clean %>%
  dplyr::select(maturity_histology, maturity_GSI) %>%
  dplyr::group_by(maturity_histology, maturity_GSI) %>%
  dplyr::summarise(count = n()) -> table4 

## filter GSI = 0
data %>% filter(!(maturation_status_histology %in% c("", NA, "no slide", "no score", "4-1", "3-4", "2-3", "1-2", "7", "6"))) %>%
  filter(!(sex_histology %in% c("Male"))) %>%
  filter(!(sample_no %in% c("743", "629", "334", "289"))) %>%
  dplyr::mutate(maturity_histology = ifelse(maturation_status_histology == 1, "immature" , "mature")) %>%
  mutate(GSI=as.numeric(GSI)) %>%
  dplyr::mutate(maturity_GSI = ifelse(GSI > 0.05, "mature" , "immature")) %>%
  filter(!(GSI %in% c(0))) -> data_clean

data_clean %>%
  dplyr::select(maturity_histology, GSI, sex_histology) %>%
  dplyr::group_by(maturity_histology, sex_histology) %>%
  dplyr::summarise(n = n(),
            median = median(GSI, na.rm=T), mean =mean(GSI),sd=sd(GSI),se=sd(GSI)/sqrt(n()), min = min(GSI), max =max(GSI))%>%
  as.data.frame() -> table5  #average, se, min, max of GSI

# boxplot of GSI cutoff by sex
data %>% filter(!(maturation_status_histology %in% c("", NA, "no slide", "no score", "4-1", "3-4", "2-3", "1-2", "7", "6"))) %>%
  filter(!(sex_histology %in% c("Male"))) %>%
  filter(!(sample_no %in% c("743", "629", "334", "289"))) %>%
  dplyr::mutate(maturity_histology = ifelse(maturation_status_histology == 1, "immature" , "mature")) %>%
  mutate(GSI=as.numeric(GSI)) %>%
  dplyr::mutate(maturity_GSI = ifelse(GSI > 0.05, "mature" , "immature")) %>%
ggplot(., aes(maturity_histology, GSI, colour = maturity_histology)) +
  geom_boxplot(position=position_dodge(width=0.9)) +
  stat_summary( fun.data = fun_length, geom = "text", 
                position=position_dodge(width=0.9), vjust=3, family="Times New Roman", size=3) +
  labs(x="", y = "GSI") + theme(legend.position="none")+scale_colour_grey(start = .1, end = .1)
ggsave(filename = 'figs/gsi_boxplots.png', dpi =200, width=8, height=6, units = "in")

# GSI regression (females)
data %>% filter(!(maturation_status_histology %in% c("", NA, "no slide", "no score", "4-1", "3-4", "2-3", "1-2", "7", "6"))) %>%
  filter(!(sex_histology %in% c("Male"))) %>%
  filter(!(sample_no %in% c("743", "629", "334", "289"))) %>%
  dplyr::mutate(maturity_histology = ifelse(maturation_status_histology == 1, "0" , "1")) %>%
  mutate(GSI=as.numeric(GSI)) %>%
  dplyr::mutate(maturity_GSI = ifelse(GSI > 0.05, "mature" , "immature")) %>%
  filter(!(GSI %in% c(0))) %>%
  arrange(maturity_histology, GSI)-> data_clean

## cut-off of 1.4%
data_clean %>%
  filter(!(maturity_histology %in% c(0))) %>%
  mutate(value_I = ifelse(GSI>0.014,1,0)) %>%
  dplyr::group_by(maturity_histology) %>%
  dplyr::summarise(n = n(),
            value1 = sum(value_I)/n(),
            sumv=sum(value_I),
            min=min(GSI))-> cutoff_1.4_mature


data_clean %>%
  filter(!(maturity_histology %in% c(1))) %>%
  mutate(value_I = ifelse(GSI>0.014,1,0)) %>%
  dplyr::group_by(maturity_histology) %>%
  dplyr::summarise(n = n(),
            value1 = sum(value_I)/n(),
            sumv=sum(value_I),
            min=min(GSI))-> cutoff_1.4_immature

## cut-off of 0.6%
data_clean %>%
  filter(!(maturity_histology %in% c(0))) %>%
  mutate(value_I = ifelse(GSI>0.006,1,0)) %>%
  dplyr::group_by(maturity_histology) %>%
  dplyr::summarise(n = n(),
            value1 = sum(value_I)/n(),
            sumv=sum(value_I),
            min=min(GSI))-> cutoff_0.6_mature


data_clean %>%
  filter(!(maturity_histology %in% c(1))) %>%
  mutate(value_I = ifelse(GSI>0.006,1,0)) %>%
  dplyr::group_by(maturity_histology) %>%
  dplyr::summarise(n = n(),
            value1 = sum(value_I)/n(),
            sumv=sum(value_I),
            min=min(GSI))-> cutoff_0.6_immature

data_clean %>% 
  mutate(log_gsi = log(GSI),
         log_length = log(length_mm)) -> log_data

log_data %>% 
  do(A1 = lm(log_gsi ~ log_length, data = log_data[log_data$maturity_histology==1,]),
     A2 = lm(log_gsi ~ log_length, data = log_data[log_data$maturity_histology==0,])) -> lm_out

lm_out %>% 
  augment(A1) %>% 
  mutate(gsi = exp(log_gsi), 
  length = exp(log_length), 
  fit = exp(.fitted)) -> data1 
lm_out %>% 
  augment(A2) %>% 
  mutate(gsi = exp(log_gsi), 
         length = exp(log_length), 
         fit = exp(.fitted)) -> data2

# GSI frequency and cut-off plot
ggplot() +
  geom_point(data=data1, aes(x = length, y = gsi), color ="grey50") + 
  geom_point(data=data2, aes(length, gsi), color = "black")+ 
  labs(y = "GSI", x =  "Length (mm)") + geom_smooth(method= 'lm', data=data1,aes(length, gsi), col = "black") +
  geom_smooth(method= 'lm', data=data2, aes(length, gsi), col = "black") +
  geom_hline(yintercept = 0.014, linetype="dotted", color = "grey50", size=1) +
  geom_hline(yintercept = 0.006, linetype="dotted", color = "grey50", size=1) +
  geom_text(aes(x =135, y= 0.10, label="A)"),family="Times New Roman", colour="black", size=5) +
  geom_text(aes(x =225, y= 0.055, label="mature"),family="Times New Roman", colour="black") +
  geom_text(aes(x =200, y= 0.005, label="immature"),family="Times New Roman", colour="black") +
  geom_text(aes(x =136, y= 0.017, label="cut-off = 1.4%"),family="Times New Roman", colour="black", size=3) +
  geom_text(aes(x =136, y= 0.009, label="cut-off = 0.6%"),family="Times New Roman", colour="black", size=3) ->plot1

data %>% filter(!(maturation_status_histology %in% c("", NA, "no slide", "no score", "4-1", "3-4", "2-3", "1-2", "7", "6"))) %>%
  filter(!(sex_histology %in% c("Male"))) %>%
  filter(!(sample_no %in% c("743", "629", "334", "289"))) %>%
  dplyr::mutate(maturity_histology = ifelse(maturation_status_histology == 1, "immature" , "mature")) %>%
  mutate(GSI=as.numeric(GSI))-> data_clean

cdat<- ddply(data_clean, "maturity_histology", summarise, rating.mean=mean(GSI))
ggplot(data_clean, aes(x=GSI, color=maturity_histology, fill=maturity_histology)) + geom_histogram(alpha=0.5, position = 'identity') +
  ylab("Frequency")+ xlab("GSI") +
  geom_text(aes(x =0, y= 100, label="B)"),family="Times New Roman", colour="black", size=5) +
  #scale_x_continuous(breaks = c(0, 0.02, 0.04, 0.06, 0.08, 0.10), limits = c(0,0.10))+
  #ggtitle("Females; n = 567") + 
  scale_color_manual(values=c("#999999", "black")) +theme(legend.title=element_blank(), legend.position=c(.85,.85)) +
  scale_fill_manual(values=c("#999999", "black" )) +geom_vline(data=cdat, aes(xintercept=rating.mean,  colour=maturity_histology),
                                                               linetype="dashed", size=1) -> plot2

cowplot::plot_grid(plot1, plot2, nrow = 1, ncol=2)
ggsave("figs/freq.png", dpi = 500, height = 6, width =10, units = "in")

#bubble plot of age versus histology
spawn  %>%
  mutate(date = mdy(date),
  month = month(date),
         day=day(date)) -> data1
data1$date2 <- paste(data1$year, data1$month, data1$day, sep="-") %>% ymd() %>% as.Date()
data1 %>% 
  mutate(julian = yday(date2),
         spawn1 = replace(spawn, spawn == 0, NA)) %>%
  dplyr::select(year, julian, spawn, spawn1) -> data1

axisx <- tickr(data1, year, 2)
axisy <- tickr(data1, julian, 5)
ggplot(data=data1, aes(x = year, y = julian, size = spawn1)) +
  geom_point(shape = 21, colour = "black", fill = "black") +
  scale_size(range = c(1, 5)) +
  xlab('') +
  ylab('Julian Day') +
  guides(size = FALSE) +
  scale_x_continuous(breaks = axisx$breaks, labels = axisx$labels) +
  scale_y_continuous(breaks = axisy$breaks, labels = axisy$labels) 
ggsave(filename = 'figs/bubbleplot.png', dpi =200, width=8, height=6, units = "in")


merge %>%
  filter(age == 2) %>%
  ggplot(., aes(x=anu_adj, color=maturity, fill=maturity)) + geom_histogram(alpha=0.5, position = 'identity') +
  ylab("Frequency")+ xlab("Outer increment measurement (mm)") +
  scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2), limits = c(0,2))+
  ggtitle("Age 2; n = 64") + theme(legend.position="none") +
  scale_color_manual(values=c("#999999", "black")) +theme(legend.title=element_blank(), legend.position=c(.15,.85)) +
  scale_fill_manual(values=c("#999999", "black" )) -> plot1

