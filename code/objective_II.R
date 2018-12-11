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
                  "readxl",
                  "Matching",
                  "devtools",
                  "FNGr",
                  "cowplot"))

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
data %>% filter(maturation_status_histology != "") %>% 
  filter(maturation_status_histology != "NA") %>% 
  filter(maturation_status_histology != "no slide") %>% 
  filter(maturation_status_histology != "no score") %>% 
  filter(maturation_status_histology != "4-1") %>% 
  filter(maturation_status_histology != "3-4") %>%
  filter(maturation_status_histology != "2-3") %>% 
  filter(maturation_status_histology != "1-2")  -> data_clean

data.clean %>%
  mutate (age = as.numeric(age)) %>%
  group_by(age, maturation_status_histology) %>%
  summarise(count = n()) -> table1 

data.clean %>%
  group_by(maturity_state_GSI, maturation_status_histology) %>%
  summarise(count = n()) -> table2 


data.clean %>%
  group_by(maturity_state_field, maturation_status_histology) %>%
  summarise(count = n()) -> table3 

data %>% filter(sex_histology != "") %>% 
  filter(sex_histology != "NA") %>% 
  filter(sex_histology != "no slide") -> data_clean

data_clean %>%
  mutate (date = as.Date(catch_date, format='%m/%d/%Y')) %>%
  group_by(maturity_state_GSI, date, sex_histology) %>%
  summarise(count = n(), 
            mean = mean (GSI)) -> table4 

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
  labs(x="Age", y = "Maturation Status (histology)") +  scale_size_continuous(name = "",
                                                                                breaks = bb,
                                                                                limits = c(0, 250),
                                                                                labels = ll,
                                                                                range = c(0, 10))
ggsave(filename = 'figs/bubbleplot.png', dpi =200, width=6, height=8, units = "in")



#confidence intervals for a difference in proportions (package PropCIs) (see package for refs)
#diffscoreci(x1, n1, x2, n2, conf.level)
#wald2ci(x1, n1, x2, n2, conf.level, adjust)
#prop.test (MASS package)
prop.test(x=c(122,99), n=c(600,400), correct=F)

diffscoreci(7, 21, 13, 17, conf.level=0.95)

wald2ci(7, 21, 13, 17, conf.level=0.95, adjust = "Wald")