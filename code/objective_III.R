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
library(broom)
library(lubridate)
devtools::install_github("ben-williams/FNGr")
theme_set(theme_sleek())

# load data ---- 
data <- read.csv("data/data_11_26_2018.csv", check.names = FALSE) 
spawn <- read.csv("data/spawn_timing_Sitka.csv", check.names = FALSE) 

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
  filter(!(GSI %in% c(0))) %>%
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

res <- prop.test(x = c(113, 484), n = c(501, 501),  correct = FALSE)
res #Pearson's and p-value output

data %>% filter(!(maturation_status_histology %in% c("", NA, "no slide", "no score", "4-1", "3-4", "2-3", "1-2", "7", "6"))) %>%
  mutate(maturation_status_histology = as.numeric(maturation_status_histology),
         GSI = as.numeric(GSI)) %>%
  filter(!(GSI %in% c(0))) %>%
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
  group_by(status_h, sex_histology, catch_date) %>%
  summarise(n = n(),
            median = median(GSI), mean =mean(GSI),sd=sd(GSI),se=sd(GSI)/sqrt(n()), min = min(GSI))%>%
  as.data.frame() -> tablex  

data_clean %>%
  dplyr::select(status_h, GSI, sex_histology) %>%
  group_by(status_h,sex_histology) %>%
  summarise(n = n(),
            median = median(GSI, na.rm=T), mean =mean(GSI, na.rm=T),sd=sd(GSI, na.rm=T),se=sd(GSI, na.rm=T)/sqrt(n()),min = min(GSI,na.rm=T))%>%
  as.data.frame() -> tabley 

#Figures 
#Boxplot of GSI cutoff by sex
data %>% filter(!(maturation_status_histology %in% c("", NA, "no slide", "no score", "4-1", "3-4", "2-3", "1-2", "7", "6"))) %>%
  filter(!(GSI %in% c(0))) %>%
  mutate(maturation_status_histology = as.numeric(maturation_status_histology),
         GSI = as.numeric(GSI)) %>%
  mutate(status_h = ifelse(maturation_status_histology>2, 'Mature','Immature')) -> data_clean

data_clean %>% filter(sex_histology == "Female" )->female
female<- as.data.frame(female)

females<-ggplot(female, aes(status_h, GSI, colour = status_h)) +
  geom_boxplot(position=position_dodge(width=0.9)) +
  stat_summary( fun.data = fun_length, geom = "text", 
                position=position_dodge(width=0.9), vjust=3, family="Times New Roman", size=3)+
  labs(x="", y = "GSI") + theme(legend.position="none")+scale_colour_grey(start = .1, end = .1)+
  geom_text(aes(x =10/25/2017, y= 0.075, label="Females"),hjust = -0.6,  family="Times New Roman", colour="black")

data_clean %>% filter(sex_histology == "Male") -> male 
      
males<-ggplot(male, aes(status_h, GSI, colour = status_h)) +
  geom_boxplot(position=position_dodge(width=0.9)) +
  stat_summary(fun.data = fun_length, geom = "text", 
                position=position_dodge(width=0.9), vjust=2, family="Times New Roman", size=3)+
  labs(x="", y = "GSI") + theme(legend.position="none")+scale_colour_grey(start = .1, end = .1)+
  geom_text(aes(x =10/25/2017, y= 0.075, label="Males"),hjust = -0.6,  family="Times New Roman", colour="black")
cowplot::plot_grid(females, males, nrow = 1, ncol=2)
ggsave(filename = 'figs/gsi_boxplots.png', dpi =200, width=8, height=6, units = "in")

#GSI sample size by females only
data %>% filter(!(maturation_status_histology %in% c("", NA, "no slide", "no score", "4-1", "3-4", "2-3", "1-2", "7", "6"))) %>%
  filter(!(GSI %in% c(0))) %>%
  mutate(maturation_status_histology = as.numeric(maturation_status_histology),
         GSI = as.numeric(GSI)) %>%
  mutate(status_h = ifelse(maturation_status_histology>2, 'Mature','Immature')) -> data_clean

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

#regression (females)
data %>% 
  filter(!(maturation_status_histology %in% c("", NA, "no slide", "no score", "4-1", "3-4", "2-3", "1-2", "7", "6"))) %>%
  filter(!(GSI %in% c(0))) %>%
  mutate(maturation_status_histology = as.numeric(maturation_status_histology),
         GSI = as.numeric(GSI)) %>%
  mutate(status_h = ifelse(maturation_status_histology>2, 1,0)) %>% 
  filter(sex_histology == "Female") %>%
  arrange(status_h, GSI)-> data_clean

data_clean %>%
  mutate(value_I = ifelse(GSI>0.016,1,0)) %>%
  group_by(status_h)%>%
  summarise(n = n(),
            value1 = sum(value_I)/n(),
            value2 =(n()-sum(value_I))/n(),
            value_I = sum(value_I),
            value_M = n()-value_I,
            min=min(GSI))-> cutoff_1.6

data_clean %>%
  mutate(value_I = ifelse(GSI>0.012,1,0)) %>%
  group_by(status_h)%>%
  summarise(n = n(),
            value1 = sum(value_I)/n(),
            value2 =(n()-sum(value_I))/n(),
            value_I = sum(value_I),
            value_M = n()-value_I,
            min=min(GSI))-> cutoff_1.2
  

data_clean %>% 
  mutate(log_gsi = log(GSI),
         log_length = log(length_mm)) -> log_data

log_data %>% 
  do(A1 = lm(log_gsi ~ log_length, data = log_data[log_data$status_h==1,]),
     A2 = lm(log_gsi ~ log_length, data = log_data[log_data$status_h==0,])) -> lm_out

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
ggplot() +
  geom_point(data=data1, aes(x = length, y = gsi), color ="grey50") + 
  geom_point(data=data2, aes(length, gsi), color = "black")+ 
  labs(y = "GSI", x =  "Length (mm)") + geom_smooth(method= 'lm', data=data1,aes(length, gsi), col = "black") +
  geom_smooth(method= 'lm', data=data2, aes(length, gsi), col = "black") +
  geom_hline(yintercept = 0.016, linetype="dotted", color = "grey50", size=1) +
  geom_hline(yintercept = 0.012, linetype="dotted", color = "grey50", size=1) +
  geom_text(aes(x =135, y= 0.07, label="Females"),family="Times New Roman", colour="black", size=5) +
  geom_text(aes(x =225, y= 0.05, label="mature"),family="Times New Roman", colour="black") +
  geom_text(aes(x =135, y= 0.007, label="immature"),family="Times New Roman", colour="black") +
  geom_text(aes(x =136, y= 0.0175, label="cut-off = 1.6%"),family="Times New Roman", colour="black", size=3) +
  geom_text(aes(x =136, y= 0.0135, label="cut-off = 1.2%"),family="Times New Roman", colour="black", size=3)

ggsave(file="figs/cut_off.png", dpi=500, width=6, height=4)


tickr_length <- data.frame(length_millimeters = 200:1000)
axisf <- tickr(tickr_length, length_millimeters, 100)

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
