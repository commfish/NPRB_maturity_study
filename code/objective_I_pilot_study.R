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
                  "magrittr", 
                  "arm",
                  "nlme", 
                  "multcomp"))
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
library(arm)
library(nlme)
library(multcomp)

windowsFonts(Times=windowsFont("TT Times New Roman"))
theme_set(theme_bw(base_size=12, base_family='Times New Roman') +
            theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  #legend.title=element_blank(),
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


#Test #1: Back Calculation of Length by BSH and SPH Methods Compared to Preferred Area
#http://derekogle.com/IFAR/supplements/backcalculation/
data %>% 
  dplyr::select(id, fish, agecap, lencap, anu1, anu2, anu3, anu4, anu5, anu6, anu7, anu8, anu9,rad1, rad2, rad3, rad4, rad5, rad6, rad7, rad8, rad9, radcap, zone) -> dataR 

# sample one from each fish and zone combo since have three samples from each fish/zone combo
set.seed(167) #set seed keeps random sample
dataR %>% 
  group_by(fish, zone) %>% 
  sample_n(1) -> sample1

#create long dataframe from wide
sample1 <-gather(sample1,agei,radi, rad1:rad9) %>%
  arrange(id,agei)%>% 
  dplyr::select(id, fish, agecap, lencap, radcap, agei, radi, zone)-> sample1
stringr::str_sub(sample1$agei,1,3)<-"" #annulus

#delete rows with NA and where agei>agecap (radial measurements with no information)
sample1 %<>% mutate(agei=as.numeric(agei)) %>%
  filterD(!is.na(radi)) %>%
  filterD(agei<=agecap) 

#calculate starting values for back-calculation methods based on just zone A1 "preferred area"
dataR %>%
  filter(zone == "A1")-> dataR
lm.sl <- lm(radcap~lencap,data=dataR)
a <- coef(lm.sl)[[1]] 
b <- coef(lm.sl)[[2]] 
lm.ls <- lm(lencap~radcap,data=dataR)
c <- coef(lm.sl)[[1]] 
d <- coef(lm.sl)[[2]] 

#back-calculation methods with ratios
#Francis 1990 pg. 897 recommends the SPH and BPH methods; the difference btw the back-calculated lengths be taken as a 
#minimum meaure of imprecision of back-calculation
sample1 %<>% mutate(SPH.len=(-a/b)+(lencap+a/b)*(radi/radcap), #Scale Proportional Hypothesis (Hile 1941:212)
                    BPH.len=lencap*(c+d*radi)/(c+d*radcap)) #Body Proportional Hypothesis (Whitney and Carlander 1956) 

sample1 %>% dplyr::select(fish, agei, zone, SPH.len) %>%
             spread(key = zone, value = SPH.len) -> data_wide_SPH

sample1 %>% dplyr::select(fish, agei, zone, BPH.len) %>%
  spread(key = zone, value = BPH.len) -> data_wide_BPH

#calculate difference in back-calculated zones in percents from zone A1
#http://www2.phy.ilstu.edu/~wenning/slh/Percent%20Difference%20Error.pdf
#SCALE PROPOTIONAL HYPOTHESIS
data_wide_SPH %<>% 
  mutate(A2= abs((A1-A2)/A1)*100,
         A3= abs((A1-A3)/A1)*100,
         A4= abs((A1-A4)/A1)*100,
         A6= abs((A1-A6)/A1)*100) %>% 
  dplyr::select(fish, agei,A2, A3, A4, A6) 

data_wide_SPH %<>% gather(variable, value, -agei, -fish) %>% 
  group_by(agei, variable) %>% 
  summarise(mean.SPH=mean(value, na.rm=T),
            sd.SPH=sd(value, na.rm=T),
            n.SPH=n(),
            se.SPH=sd(value, na.rm=T)/sqrt(n()))%>%
  mutate(age = as.factor(agei),
         zone=variable)

#BODY PROPOTIONAL HYPOTHESIS
data_wide_BPH %<>%
  mutate(A2= abs((A1-A2)/A1)*100,
         A3= abs((A1-A3)/A1)*100,
         A4= abs((A1-A4)/A1)*100,
         A6= abs((A1-A6)/A1)*100) %>% 
  dplyr::select(fish, agei,A2, A3, A4, A6)

data_wide_BPH %<>% gather(variable, value, -agei, -fish) %>% 
  group_by(agei, variable) %>% 
  summarise(mean.BPH=mean(value, na.rm=T),
            sd.BPH=sd(value, na.rm=T),
            n.BPH=n(),
            se.BPH=sd(value, na.rm=T)/sqrt(n())) %>% 
  mutate(age = as.factor(agei),
         zone=variable) -> data_wide_BPH 

#plots by BPH and SPH
tickr_length <- data.frame(mean.SPH = 0:8)
axisb <- tickr(tickr_length, mean.SPH, 0.5)
ggplot(data = data_wide_SPH, aes(x = age, y = mean.SPH)) +
  geom_bar(aes(fill=zone), stat="identity", position="dodge",alpha=0.9) +
  scale_fill_grey(start = 0, end = .8)+theme(legend.position=c(.9,.75))+
  geom_text(aes(x = 3, y=8, label="A) Scale Proportional Hypothesis"), family="Times New Roman")+
  scale_y_continuous(breaks = axisb$breaks, labels = axisb$labels) +
  labs(x = "Age (Annulus)", y =  "Mean Percent Difference from Length Calculated at Zone A1")-> SPH

tickr_length <- data.frame(mean = 0:3)
axisb <- tickr(tickr_length, mean, 0.5)
ggplot(data = data_wide_BPH, aes(x = age, y = mean.BPH)) +
  geom_bar(aes(fill=zone), stat="identity", position="dodge",alpha=0.9) +
  scale_fill_grey(start = 0, end = .8)+theme(legend.position=c(.9,.75))+
  geom_text(aes(x = 3, y=3, label="B) Body Proportional Hypothesis"), family="Times New Roman")+
  scale_y_continuous(breaks = axisb$breaks, labels = axisb$labels) +
  labs(x = "Age (Annulus)", y =  "Mean Percent Difference from Length Calculated at Zone A1")-> BPH

cowplot::plot_grid(SPH, BPH,   align = "h", nrow = 1, ncol=2) 
ggsave("figs/SPH_BPH2.png", dpi = 300, height = 6, width = 10, units = "in")

#Test #2: ANOVA with repeated treatments (between and within group variability) with multilevels 
          #(a regression that allows for the errors to be dependent on eachother (as our conditions of Valence were repeated within each participant). 
#https://sapa-project.org/blog/2013/06/28/repeated-measures-anova-in-r/
model.fixed = gls(aprop1 ~ zone, data=data, method="ML") #fixed effects only
lm_simple <- lm(aprop1 ~ zone, data = data) #fixed effects only

baseline = lme(aprop1 ~ 1, random = ~1|fish/zone, data=data, method="ML")
zone_model = lme(aprop1 ~ zone, random = ~1|fish/zone, data=data, method="ML")

anova(baseline, zone_model) #zone had a significant impact on the measured aprop1 of the fish
posthoc <- glht(zone_model, linfct = mcp(zone = "Tukey"))

summary(lm_simple)
summary(baseline)
summary(zone_model)
summary(posthoc)

library(psych)
describeBy(data$aprop1, group = data$zone)
# y= proportions, x = prop1, propr 2, and bars are A1 through A6

#figure for proportions by zone for all increments
qplot(zone, value, data = data, fill=animals)+ 
  geom_boxplot() + facet_grid(~region) + scale_fill_brewer()