#notes---
#author: Sara E Miller 
#contact: sara.miller@alaska.gov; 907-465-4245
#Last edited: December, 2018

# load ----
is_installed <- function(mypkg){ is.element(mypkg, installed.packages()[,1])}

load_or_install <- function(package_names){
  for(package_name in package_names){
    if(!is_installed(package_name)){install.packages(package_name,repos="http://lib.stat.cmu.edu/R/CRAN")}
    library(package_name,character.only=TRUE,quietly=TRUE,verbose=FALSE)
  }}

load_or_install(c("extrafont",
                  "tidyverse",
                  "dgof",
                  "Matching",
                  "devtools",
                  "FNGr",
                  "cowplot",
                  "psych"))
library(extrafont)
library(tidyverse)
library(dgof)
library(Matching)
library(devtools)
library(FNGr)
library(cowplot)
library(psych)

tickr <- function(
  data, # dataframe
  var, # column of interest
  to # break point definition 
){
  
  VAR <- enquo(var) # makes VAR a dynamic variable
  
  data %>% 
    distinct(!!VAR) %>%
    #    ungroup(!!VAR) %>% 
    mutate(labels = ifelse(!!VAR %in% seq(to * round(min(!!VAR) / to), max(!!VAR), to),
                           !!VAR, "")) %>%
    dplyr::select(breaks = UQ(VAR), labels)
}


# load data ---- 
data <- read.csv("data/data_11_26_2018.csv", check.names = FALSE) 

#clean data ----
data %>% filter(!(maturation_status_histology %in% c("", NA, "no slide", "no score", "4-1", "3-4", "2-3", "1-2"))) %>%
         filter(!(age %in% c(NA)))  -> data_clean

data_clean %>%
  mutate (age = as.numeric(age)) %>%
  group_by(age, maturation_status_histology) %>%
  summarise(count = n()) -> table1 

data %>% filter(!(maturation_status_histology %in% c("", NA, "no slide", "no score", "4-1", "3-4", "2-3", "1-2")))-> data_clean

data_clean %>%
  group_by(maturity_state_GSI, maturation_status_histology) %>%
  summarise(count = n()) -> table2 

data_clean %>%
  group_by(maturity_state_field, maturation_status_histology) %>%
  summarise(count = n()) -> table3 

data %>% filter(!(sex_histology %in% c("", NA, "no slide")))  -> data_clean

maturity_state_field

#figures 
loadfonts(device="win")
windowsFonts(Times=windowsFont("TT Times New Roman"))
theme_set(theme_bw(base_size=12, base_family='Times New Roman') +
            theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  plot.title = element_text(family='Times New Roman', hjust = 0.5, size=12)))
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

#cohen kappa evaluation
data %>% filter(!(maturation_status_histology %in% c("", NA, "no slide", "no score", "4-1", "3-4", "2-3", "1-2", "7"))) %>%
dplyr::select(maturation_status_histology = maturation_status_histology,
              maturity_state_field = maturity_state_field) -> data_clean

cohen.kappa(data_clean)

#Two proportions z-test with continuity correction
#The function returns:
  #the value of Pearson’s chi-squared test statistic.
  #a p-value
  #a 95% confidence intervals
  #an estimated probability of success 

data_clean %>% mutate(maturation_status_histology = as.numeric(maturation_status_histology),
                      maturity_state_field = as.numeric(maturity_state_field)) %>%
  mutate(status_h = ifelse(maturation_status_histology>1, 'mature','immature'),
         status_m = ifelse(maturity_state_field>1, 'mature','immature')) -> data_clean

data_clean %>%
  dplyr::select(status_h) %>%
  group_by(status_h) %>%
  summarise(count = n()) -> tablex  

data_clean %>%
  dplyr::select(status_m) %>%
  group_by(status_m) %>%
  summarise(count = n()) -> tabley

res <- prop.test(x = c(264, 141), n = c(727, 727),  correct = FALSE)
