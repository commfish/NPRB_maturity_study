#notes---
#author: Sara E Miller 
#contact: sara.miller@alaska.gov; 907-465-4245
#Last edited: February, 2019

# load ----
#devtools::install_github("ben-williams/FNGr")
#devtools::install_github('droglenc/RFishBC')
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
                  "cowplot",
                  "lme4",
                  "FSA",
                  "RFishBC",
                  "stringr",
                  "magrittr"))
library(extrafont)
library(tidyverse)
library(devtools)
library(FNGr)
library(cowplot)
library(lme4)
library(FSA)
library(RFishBC)
library(magrittr)
library(stringr)

windowsFonts(Times=windowsFont("TT Times New Roman"))
theme_set(theme_bw(base_size=12, base_family='Times New Roman') +
            theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  plot.title = element_text(family='Times New Roman', hjust = 0.5, size=12)))

asintrans <- function(p) {asin(sqrt(p))} # arctransformation function
# load data ---- 
data <- read.csv("data/prop_pilot_study.csv", check.names = FALSE) 

#clean data ----
data %>%
  filter(radcap >0)-> data  #delete samples without measurements 
data %>% #calculate proportions
  mutate(prop1=(anu1)/radcap,
         prop2=(anu1+anu2)/radcap,
         prop3=(anu1+anu2+anu3)/radcap,
         prop4=(anu1+anu2+anu3+anu4)/radcap,
         prop5=(anu1+anu2+anu3+anu4+anu5)/radcap,
         prop6=(anu1+anu2+anu3+anu4+anu5+anu6)/radcap,
         prop7=(anu1+anu2+anu3+anu4+anu5+anu6+anu7)/radcap,
         prop8=(anu1+anu2+anu3+anu4+anu5+anu6+anu7+anu8)/radcap,
         prop9=(anu1+anu2+anu3+anu4+anu5+anu6+anu7+anu8+anu9)/radcap) %>%
  mutate(aprop1=asintrans(prop1), #transform data by arctransform
         aprop2=asintrans(prop2),
         aprop3=asintrans(prop3),
         aprop4=asintrans(prop4),
         aprop5=asintrans(prop5),
         aprop6=asintrans(prop6),
         aprop7=asintrans(prop7),
         aprop8=asintrans(prop8),
         aprop9=asintrans(prop9)) %>%
mutate(rad1=(anu1),#calculate radial units (ie. focus to each annulus)
       rad2=(anu1+anu2),
       rad3=(anu1+anu2+anu3),
       rad4=(anu1+anu2+anu3+anu4),
       rad5=(anu1+anu2+anu3+anu4+anu5),
       rad6=(anu1+anu2+anu3+anu4+anu5+anu6),
       rad7=(anu1+anu2+anu3+anu4+anu5+anu6+anu7),
       rad8=(anu1+anu2+anu3+anu4+anu5+anu6+anu7+anu8),
       rad9=(anu1+anu2+anu3+anu4+anu5+anu6+anu7+anu8+anu9),
       radcap=radcap)-> data 
         
#analysis----
#summarize by region and individual fish
data %>% 
    group_by(fish, zone) %>% 
    summarise(lencap=mean(lencap),
              mean1=mean(prop1),
              mean2=mean(prop2),
              mean3=mean(prop3),
              mean4=mean(prop4),
              mean5=mean(prop5),
              mean6=mean(prop6),
              mean7=mean(prop7),
              mean8=mean(prop8),
              mean9=mean(prop9),
              amean1=mean(aprop1),
              amean2=mean(aprop2),
              amean3=mean(aprop3),
              amean4=mean(aprop4),
              amean5=mean(aprop5),
              amean6=mean(aprop6),
              amean7=mean(aprop7),
              amean8=mean(aprop8),
              amean9=mean(aprop9),
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

ggplot(data = data, aes(x = lencap, y = aprop1, colour = slide)) +
  geom_point() +
  # geom_smooth(method = "lm") +
  facet_grid(~zone) +
  labs(x = "lencap", y = 
         "Transformed prop (focus to first annulus)")-> x

ggplot(data = data, aes(x = lencap, y = aprop1, colour = zone)) +
  geom_point() +
  # geom_smooth(method = "lm") +
  facet_grid(~ zone) +
  labs(x = "lencap", y = 
         "Transformed prop (focus to first annulus)")-> y

#Test #1: mixed effects models
lm_simple_i <- lmer(aprop1 ~ lencap, data = data)
lm_simple_s <- lmer(aprop1 ~ lencap, data = data)

lm_fish_i <- lmer(aprop1 ~ lencap + (1 | fish), data = data) #Model with a random intercept and slope for each fish 
lm_fish_s <- lmer(aprop1 ~ lencap + (1 + lencap | fish), data = data)

lm_zone_i <- lmer(aprop1 ~ lencap + (1 | zone), data = data) #Model with a random intercept and slope for each zone
lm_zone_s <- lmer(aprop1 ~ lencap + (1 + lencap | zone), data = data)

# ranef(lm_diver_s)
AIC(lm_fish_i);AIC(lm_fish_s);AIC(lm_zone_i);AIC(lm_zone_s)
#The model that explained the most variation and had the lowest AIC was the one with a random intercept and slope for fish.

#Test #2: Back Calculation of Length by BSH and SPH Methods
#http://derekogle.com/IFAR/supplements/backcalculation/
data %>% 
  dplyr::select(fish, agecap, lencap, anu1, anu2, anu3, anu4, anu5, anu6, anu7, anu8, anu9,rad1, rad2, rad3, rad4, rad5, rad6, rad7, rad8, rad9, radcap, zone)->dataR 

dataR <- gather(dataR,agei,radi,anu1:anu9) %>%
  arrange(fish,agei)%>%
  dplyr::select(fish, agecap, lencap, radcap, agei, radi, zone)-> dataR 

stringr::str_sub(dataR$agei,1,3)<-""

dataR %<>% mutate(agei=as.numeric(agei)) %>%
  filterD(!is.na(radi)) %>%
  filterD(agei<=agecap) 

#calculate starting values
lm.sl <- lm(radcap~lencap,data=dataR)
a <- coef(lm.sl)[[1]] 
b <- coef(lm.sl)[[2]] 
lm.ls <- lm(lencap~radcap,data=dataR)
c <- coef(lm.sl)[[1]] 
d <- coef(lm.sl)[[2]] 

#Francis 1990 pg. 897 recommends the SPH and BPH methods; the difference btw the back-calculated lengths be taken as a 
#minimum meaure of imprecision of back-calculation
dataR %<>% mutate(#DL.len=(radi/radcap)*lencap, #Dahl Lee
                  #FL.len=(radi/radcap)*(lencap-c)+c, #Fraser-Lee
                    SPH.len=(-a/b)+(lencap+a/b)*(radi/radcap), #Scale Proportional Hypothesis (Hile 1941:212)
                    BPH.len=lencap*(c+d*radi)/(c+d*radcap)) #Body Proportional Hypothesis (Whitney and Carlander 1956) 

#summary of data by fish, zone, and age (three samples/zone and age for each fish)
dataR %>%
  group_by(fish, zone, agei) %>%
  summarize(n.SPH=validn(SPH.len),
            n.BPH=validn(BPH.len),
            mean.SPH=round(mean(SPH.len),0),
            mean.BPH=round(mean(BPH.len),0),
            sd.SPH=round(sd(SPH.len),1),
            sd.BPH=round(sd(BPH.len),1)) %>%
  as.data.frame() -> data_long

#reshape data so can calculate the difference in zones
data_long %>% dplyr::select(fish, agei, zone, mean.SPH) %>%
              spread(key = zone, value = mean.SPH) -> data_wide_SPH

data_long %>% dplyr::select(fish, agei, zone, n.SPH) %>%
  spread(key = zone, value = n.SPH) -> data_wide_nSPH

data_long %>% dplyr::select(fish, agei, zone, mean.BPH) %>%
  spread(key = zone, value = mean.BPH) -> data_wide_BPH

data_long %>% dplyr::select(fish, agei, zone, n.BPH) %>%
  spread(key = zone, value = n.BPH) -> data_wide_nBPH

#calcualte difference in back-calcualted zones
#http://www2.phy.ilstu.edu/~wenning/slh/Percent%20Difference%20Error.pdf
data_wide_SPH %>% 
  mutate(A1_A2= ((A1-A2)/(0.5*(A1+A2))*100),
         A2_A3= ((A2-A3)/(0.5*(A2+A3))*100),
         A3_A4= ((A3-A4)/(0.5*(A3+A4))*100),
         A4_A6= ((A4-A6)/(0.5*(A4+A6))*100),
         A1_A3= ((A1-A3)/(0.5*(A1+A3))*100),
         A1_A4= ((A1-A4)/(0.5*(A1+A4))*100),
         A1_A6= ((A1-A6)/(0.5*(A1+A6))*100),
         A2_A4= ((A2-A4)/(0.5*(A2+A4))*100),
         A2_A6= ((A2-A6)/(0.5*(A2+A6))*100),
         A3_A6= ((A3-A6)/(0.5*(A3+A6))*100))-> data_wide_SPH

data_wide_BPH %>% 
  mutate(A1_A2= ((A1-A2)/(0.5*(A1+A2))*100),
         A2_A3= ((A2-A3)/(0.5*(A2+A3))*100),
         A3_A4= ((A3-A4)/(0.5*(A3+A4))*100),
         A4_A6= ((A4-A6)/(0.5*(A4+A6))*100),
         A1_A3= ((A1-A3)/(0.5*(A1+A3))*100),
         A1_A4= ((A1-A4)/(0.5*(A1+A4))*100),
         A1_A6= ((A1-A6)/(0.5*(A1+A6))*100),
         A2_A4= ((A2-A4)/(0.5*(A2+A4))*100),
         A2_A6= ((A2-A6)/(0.5*(A2+A6))*100),
         A3_A6= ((A3-A6)/(0.5*(A3+A6))*100))-> data_wide_BPH

data_wide_SPH %>% dplyr::select(fish, agei,A1_A2, A2_A3, A3_A4,
                                A4_A6, A1_A3, A1_A4, A1_A6, A2_A4, A2_A6,
                                A3_A6) -> data_wide_SPH
data_wide_SPH %>% gather(variable, value, -agei, -fish) -> data_wide_SPH
  
data_wide_BPH %>% dplyr::select(fish, agei,A1_A2, A2_A3, A3_A4,
                                A4_A6, A1_A3, A1_A4, A1_A6, A2_A4, A2_A6,
                                A3_A6) -> data_wide_BPH
data_wide_BPH %>% gather(variable, value, -agei, -fish) -> data_wide_BPH

ggplot(data = data_wide_SPH, aes(x = agei, y = aprop1, colour = zone)) +
  geom_bar() +
  # geom_smooth(method = "lm") +
  facet_grid(~ zone) +
  labs(x = "lencap", y = 
         "Transformed prop (focus to first annulus)")-> y