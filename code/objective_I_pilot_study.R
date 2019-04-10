#notes---
#author: Sara E Miller 
#contact: sara.miller@alaska.gov; 907-465-4245
#Last edited: February, 2019

# load libraries----
devtools::install_github("ben-williams/FNGr")
devtools::install_github('droglenc/RFishBC')
devtools::install_github('droglenc/FSA')

#font_import() #only do this one time - it takes a while
loadfonts(device="win")

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
library(FinCal)
library(msmtools)
library(broom)
library(psych)
library(WRS2)
library(cowplot)
library(ggbiplot)
library(vqv)
library(nlme)

windowsFonts(Times=windowsFont("Times New Roman"))
theme_sleek <- function(base_size = 12, base_family = "Times") {
  half_line <- base_size/2
  theme_light(base_size = 12, base_family = "Times") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks.length = unit(half_line / 2.2, "pt"),
      strip.background = element_rect(fill = NA, colour = NA),
      strip.text.x = element_text(colour = "black"),
      strip.text.y = element_text(colour = "black"),
      panel.border = element_rect(fill = NA),
      legend.text=element_text(size=12),
      legend.key.size = unit(0.9, "lines"),
      legend.key = element_rect(colour = NA, fill = NA),
      legend.background = element_rect(colour = NA, fill = NA)
    )
}

theme_set(theme_sleek())

asintrans <- function(p) {asin(sqrt(p))} # arctransformation function

# load data ---- 
data <- read.csv("data/prop_pilot_study.csv", check.names = FALSE) 

#clean data ----
data %>%
  filter(radcap >0) %>% #delete samples without measurements 
#calculate proportions
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

ggplot(data = data, aes(x = lencap, y = aprop1, colour = specimen)) +
  geom_point() +
  facet_grid(~zone) +
  labs(x = "lencap", y = 
         "Transformed prop (focus to first annulus)")-> x

ggplot(data = data, aes(x = lencap, y = aprop1, colour = zone)) +
  geom_point() +
  facet_grid(~ zone) +
  labs(x = "lencap", y = 
         "Transformed prop (focus to first annulus)")-> y


#Test #1: Back Calculation of Length by BPH and SPH Methods Compared to Preferred Area
#http://derekogle.com/IFAR/supplements/backcalculation/
data %>% 
  dplyr::select(id, fish, agecap, lencap, anu1, anu2, anu3, anu4, anu5, anu6, anu7, anu8, anu9,rad1, rad2, rad3, rad4, rad5, rad6, rad7, rad8, rad9, radcap, zone) -> dataR 

# sample one from each fish and zone combo since have three samples from each fish/zone combo
set.seed(167) #set seed keeps random sample the same
dataR %>% 
  group_by(fish, zone) %>% 
  sample_n(1) -> sample1 #sample one scale from each fish zone (5 scales sampled from each fish)

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
#back-calculation methods with ratios
#Francis 1990 pg. 897 recommends the SPH and BPH methods; the difference btw the back-calculated lengths be taken as a 
#minimum meaure of imprecision of back-calculation
#calculate a, b, c, d for each zone
sample1 %>%
  filter(zone == "A1" & agei == 1)-> lm_data
sample1 %>%
  filter(agei == 1) %>% 
  as.data.frame() -> lm_data_zone

lm.sl <- lm(radcap~lencap,data=lm_data)
a <- coef(lm.sl)[[1]] 
b <- coef(lm.sl)[[2]] 
lm.ls <- lm(lencap~radcap,data=lm_data)
c <- coef(lm.ls)[[1]] 
d <- coef(lm.ls)[[2]] 

sample1 %>%
    filter(zone == "A1") %>%
    mutate(SPH.len=(-a/b)+(lencap+a/b)*(radi/radcap), #Scale Proportional Hypothesis (Hile 1941:212)
                    BPH.len=lencap*((c+d*radi)/(c+d*radcap))) -> sample_A1 #Body Proportional Hypothesis (Whitney and Carlander 1956) 

sample1 %>%
  filter(zone == "A2" & agei == 1)-> lm_data
lm.sl <- lm(radcap~lencap,data=lm_data)
a <- coef(lm.sl)[[1]] 
b <- coef(lm.sl)[[2]] 
lm.ls <- lm(lencap~radcap,data=lm_data)
c <- coef(lm.ls)[[1]] 
d <- coef(lm.ls)[[2]] 

sample1 %>%
  filter(zone == "A2") %>%
  mutate(SPH.len=(-a/b)+(lencap+a/b)*(radi/radcap), 
         BPH.len=lencap*((c+d*radi)/(c+d*radcap))) -> sample_A2

sample1 %>%
  filter(zone == "A3" & agei == 1)-> lm_data #filter by agei so only one sample/region/fish
lm.sl <- lm(radcap~lencap,data=lm_data)
a <- coef(lm.sl)[[1]] 
b <- coef(lm.sl)[[2]] 
lm.ls <- lm(lencap~radcap,data=lm_data)
c <- coef(lm.ls)[[1]] 
d <- coef(lm.ls)[[2]] 

sample1 %>%
  filter(zone == "A3") %>%
  mutate(SPH.len=(-a/b)+(lencap+a/b)*(radi/radcap), 
         BPH.len=lencap*((c+d*radi)/(c+d*radcap))) -> sample_A3

sample1 %>%
  filter(zone == "A4" & agei == 1)-> lm_data
lm.sl <- lm(radcap~lencap,data=lm_data)
a <- coef(lm.sl)[[1]] 
b <- coef(lm.sl)[[2]] 
lm.ls <- lm(lencap~radcap,data=lm_data)
c <- coef(lm.ls)[[1]] 
d <- coef(lm.ls)[[2]] 

sample1 %>%
  filter(zone == "A4") %>%
  mutate(SPH.len=(-a/b)+(lencap+a/b)*(radi/radcap), 
         BPH.len=lencap*((c+d*radi)/(c+d*radcap))) -> sample_A4

sample1 %>%
  filter(zone == "A6" & agei == 1)-> lm_data
lm.sl <- lm(radcap~lencap,data=lm_data)
a <- coef(lm.sl)[[1]] 
b <- coef(lm.sl)[[2]] 
lm.ls <- lm(lencap~radcap,data=lm_data)
c <- coef(lm.ls)[[1]] 
d <- coef(lm.ls)[[2]] 

sample1 %>%
  filter(zone == "A6") %>%
  mutate(SPH.len=(-a/b)+(lencap+a/b)*(radi/radcap), 
         BPH.len=lencap*((c+d*radi)/(c+d*radcap))) -> sample_A6

x<- rbind(sample_A1, sample_A2) #combine data for all zones
x<- rbind(x, sample_A3)
x<- rbind(x, sample_A4)
x<- rbind(x, sample_A6)

x %>%
  group_by(agei, zone) %>%
  summarize(n.radcap=validn(radcap), mn.radcap=mean(radcap),sd.radcap=sd(radcap),
            ss.radcap = sum(radcap^2),cv.radcap = coefficient.variation(sd.radcap, mn.radcap),
            n.SPH=validn(SPH.len),mn.SPH.len=mean(SPH.len),sd.SPH.len=sd(SPH.len),ss.SPH.len = sum(SPH.len^2),
            cv.SPH = coefficient.variation(sd.SPH.len, mn.SPH.len),
            n.BPH=validn(BPH.len),mn.BPH.len=mean(BPH.len),sd.BPH.len=sd(BPH.len), ss.BPH.len = sum(BPH.len^2),
            cv.BPH = coefficient.variation(sd.BPH.len, mn.BPH.len)) %>%
  as.data.frame() -> sample2
write.csv(sample2, "output/table_summary_obj1.csv") 

x %>% dplyr::select(fish, agei, zone, SPH.len) %>%
             spread(key = zone, value = SPH.len) -> data_wide_SPH

  
x %>% dplyr::select(fish, agei, zone, BPH.len) %>%
  spread(key = zone, value = BPH.len) -> data_wide_BPH

#calculate difference in back-calculated zones in percents from zone A1
#http://www2.phy.ilstu.edu/~wenning/slh/Percent%20Difference%20Error.pdf

#SCALE PROPOTIONAL HYPOTHESIS
data_wide_SPH %<>% 
  mutate(A2= abs((A1-A2)/A1)*100,
         A3= abs((A1-A3)/A1)*100,
         A4= abs((A1-A4)/A1)*100,
         A6= abs((A1-A6)/A1)*100)%>% 
  dplyr::select(fish, agei,A2, A3, A4, A6) 

data_wide_SPH %>% 
  gather(variable, value, -agei, -fish) %>% 
  group_by(agei, variable) %>% 
  summarise(mean.SPH=mean(value, na.rm=T),
            sd.SPH=sd(value, na.rm=T),
            n.SPH=n(),
            se.SPH=sd(value, na.rm=T)/sqrt(n()))%>%
  mutate(age = as.factor(agei),
         zone=as.factor(variable)) -> data_wide_SPH 

#BODY PROPOTIONAL HYPOTHESIS
data_wide_BPH %<>%
  mutate(A2= abs((A1-A2)/A1)*100,
         A3= abs((A1-A3)/A1)*100,
         A4= abs((A1-A4)/A1)*100,
         A6= abs((A1-A6)/A1)*100) %>% 
  dplyr::select(fish, agei,A2, A3, A4, A6)

data_wide_BPH  %>%  
  gather(variable, value, -agei, -fish) %>% 
  group_by(agei, variable) %>% 
  summarise(mean.BPH=mean(value, na.rm=T),
            sd.BPH=sd(value, na.rm=T),
            n.BPH=n(),
            se.BPH=sd(value, na.rm=T)/sqrt(n())) %>% 
  mutate(age = as.factor(agei),
         zone=as.factor(variable)) -> data_wide_BPH 

#plots by SPH
tickr_length <- data.frame(mean.SPH = 0:18)
axisb <- tickr(tickr_length, mean.SPH, 2)
ggplot(data = data_wide_SPH, aes(x = age, y = mean.SPH)) +
  geom_bar(aes(fill=zone), stat="identity", position="dodge",alpha=0.9) +
  scale_fill_grey(start = 0, end = .8)+theme(legend.position=c(.9,.75), legend.title=element_blank())+
  annotate("text", x = 1.9, y=18, label="A) Scale Proportional Hypothesis", family="Times New Roman")+
  scale_y_continuous(breaks = axisb$breaks, labels = axisb$labels) +
  labs(x = "Age", y =  "Mean % Difference")-> SPH

tickr_length <- data.frame(mean.SPH = 0:300)
axisb <- tickr(tickr_length, mean.SPH, 50)
ggplot(data = sample2, aes(x = as.factor(agei), y = mn.SPH.len)) +
  geom_bar(aes(fill=zone), stat="identity", position="dodge",alpha=0.9) +
  scale_fill_grey(start = 0, end = .8)+theme(legend.position=c(.05,.7), legend.title=element_blank())+
  annotate("text", x = 2, y=250, label="A) Scale Proportional Hypothesis", family="Times New Roman")+
  scale_y_continuous(breaks = c(0, 50, 100, 150, 200, 250), limits = c(0,250))+
  labs(x = "Age", y =  "Back-Calculated Length (mm)")-> SPH2

#ANCOVA Models (SPH)
lm_data_zone %>% 
  do(taggingr = lm(radcap~lencap + zone, data = sample1)) -> lm_out

lm_out %>% 
  tidy(taggingr) %>% 
  write_csv("output/SPH_ancova.csv")

lm_out %>% 
  glance(taggingr) %>% 
  write_csv("output/SPH_ancova_R2.csv")

lm_out %>% 
  augment(taggingr) %>% 
  mutate(fit = (.fitted)) %>% 
  ggplot(aes(x = lencap, y = radcap)) +
  geom_point(color="grey50")+facet_wrap(~zone)+geom_line(aes(lencap, fit), color = "black") + 
  annotate("text", x = 175, y=7, label="Scale Proportional Hypothesis", family="Times New Roman")+
  scale_x_continuous(breaks = c(100, 150, 200, 250), limits = c(100,250))+
  labs(y = "Scale Radius (mm)", x =  "Capture Length (mm)")-> SPH3
ggsave("figs/SPH_ANCOVA.png", dpi = 500, height = 6, width = 8, units = "in")

#Linear models by zone (SPH)
lm_data_zone %>% 
   do(A1 = lm(radcap ~ lencap, data = sample1[sample1$zone=="A1",]),
      A2 = lm(radcap ~ lencap, data = sample1[sample1$zone=="A2",]),
      A3 = lm(radcap ~ lencap, data = sample1[sample1$zone=="A3",]),
      A4 = lm(radcap ~ lencap, data = sample1[sample1$zone=="A4",]),
      A6 = lm(radcap ~ lencap, data = sample1[sample1$zone=="A6",])) -> lm_out

lm_out %>% 
  tidy(A1) %>% 
  mutate(zone = "A1") -> A1
lm_out %>% 
  tidy(A2) %>% 
mutate(zone = "A2") -> A2
lm_out %>% 
  tidy(A3) %>% 
mutate(zone = "A3") -> A3
lm_out %>% 
  tidy(A4) %>% 
mutate(zone = "A4") -> A4
lm_out %>% 
  tidy(A6) %>% 
  mutate(zone = "A6") -> A6

x<- rbind(A1, A2) #combine data for all zones
x<- rbind(x, A3)
x<- rbind(x, A4)
x<- rbind(x, A6)
write_csv(x, "output/SPH_lm.csv")

lm_out %>% 
  glance(A1) %>% 
  mutate(zone = "A1") -> A1
lm_out %>% 
  glance(A2) %>% 
  mutate(zone = "A2") -> A2
lm_out %>% 
  glance(A3) %>% 
  mutate(zone = "A3") -> A3
lm_out %>% 
  glance(A4) %>% 
  mutate(zone = "A4") -> A4
lm_out %>% 
  glance(A6) %>% 
  mutate(zone = "A6") -> A6

x<- rbind(A1, A2) #combine data for all zones
x<- rbind(x, A3)
x<- rbind(x, A4)
x<- rbind(x, A6)
write_csv(x, "output/SPH_lm_R2.csv")

lm_out %>% 
  augment(A1) %>% 
  mutate(fit = (.fitted)) %>% 
  ggplot(aes(x = lencap, y = radcap)) +
  geom_point(color ="grey50")+geom_line(aes(lencap, fit), color = "black") + 
  annotate("text", x = 200, y=7, label="SPH (A1)", family="Times New Roman")+
  scale_x_continuous(breaks = c(150, 200, 250), limits = c(150,250))+
  labs(y = "Scale Radius (mm)", x =  "Capture Length (mm)")-> A1

lm_out %>% 
  augment(A2) %>% 
  mutate(fit = (.fitted)) %>% 
  ggplot(aes(x = lencap, y = radcap)) +
  geom_point(color ="grey50")+geom_line(aes(lencap, fit), color = "black") + 
  annotate("text", x = 200, y=7, label="SPH (A2)", family="Times New Roman")+
  scale_x_continuous(breaks = c(150, 200, 250), limits = c(150,250))+
  labs(y = "Scale Radius (mm)", x =  "Capture Length (mm)")-> A2

lm_out %>% 
  augment(A3) %>% 
  mutate(fit = (.fitted)) %>% 
  ggplot(aes(x = lencap, y = radcap)) +
  geom_point(color ="grey50")+geom_line(aes(lencap, fit), color = "black") + 
  annotate("text", x = 200, y=7, label="SPH (A3)", family="Times New Roman")+
  scale_x_continuous(breaks = c(150, 200, 250), limits = c(150,250))+
  labs(y = "Scale Radius (mm)", x =  "Capture Length (mm)")-> A3

lm_out %>% 
  augment(A4) %>% 
  mutate(fit = (.fitted)) %>% 
  ggplot(aes(x = lencap, y = radcap)) +
  geom_point(color ="grey50")+geom_line(aes(lencap, fit), color = "black") + 
  annotate("text", x = 200, y=7, label="SPH (A4)", family="Times New Roman")+
  scale_x_continuous(breaks = c(150, 200, 250), limits = c(150,250))+
  labs(y = "Scale Radius (mm)", x =  "Capture Length (mm)")-> A4

lm_out %>% 
  augment(A6) %>% 
  mutate(fit = (.fitted)) %>% 
  ggplot(aes(x = lencap, y = radcap)) +
  geom_point(color ="grey50")+geom_line(aes(lencap, fit), color = "black") + 
  annotate("text", x = 200, y=7, label="SPH (A6)", family="Times New Roman")+
  scale_x_continuous(breaks = c(150, 200, 250), limits = c(150,250))+
  labs(y = "Scale Radius (mm)", x =  "Capture Length (mm)")-> A6
cowplot::plot_grid(A1, A2, A3, A4, A6,  align = "vh", nrow = 2, ncol=3)
ggsave("figs/SPH_regression.png", dpi = 500, height = 6, width = 8, units = "in")

#plots by BPH
tickr_length <- data.frame(mean.SPH = 0:35)
axisb <- tickr(tickr_length, mean.SPH, 5)
ggplot(data = data_wide_BPH, aes(x = age, y = mean.BPH)) +
  geom_bar(aes(fill=zone), stat="identity", position="dodge",alpha=0.9) +
  scale_fill_grey(start = 0, end = .8)+theme(legend.position="none")+
  annotate("text", x = 1.9, y=35, label="B) Body Proportional Hypothesis", family="Times New Roman")+
  scale_y_continuous(breaks = axisb$breaks, labels = axisb$labels) +
  labs(x = "Age", y =  "Mean % Difference")-> BPH
cowplot::plot_grid(SPH, BPH,   align = "hv", nrow = 2, ncol=1) 
ggsave("figs/length_diff.png", dpi = 500, height = 6, width = 8, units = "in")

tickr_length <- data.frame(mean.BPH = 0:350)
axisb <- tickr(tickr_length, mean.BPH, 50)
ggplot(data = sample2, aes(x = as.factor(agei), y = mn.BPH.len)) +
  geom_bar(aes(fill=zone), stat="identity", position="dodge",alpha=0.9) +
  scale_fill_grey(start = 0, end = .8)+theme(legend.position="none")+
  annotate("text", x = 2, y=350, label="B) Body Proportional Hypothesis", family="Times New Roman")+
  scale_y_continuous(breaks = c(0, 50, 100, 150, 200, 250, 300, 350), limits = c(0,350))+
  labs(x = "Age", y =  "Back-Calculated Length (mm)")-> BPH2
cowplot::plot_grid(SPH2, BPH2,   align = "hv", nrow = 2, ncol=1) 
ggsave("figs/length.png", dpi = 500, height = 6, width = 8, units = "in")


#ANCOVA Models (BPH)
lm_data_zone %>% 
  do(taggingr = lm(lencap ~ radcap + zone, data = .)) -> lm_out

lm_out %>% 
  tidy(taggingr) %>% 
  write_csv("output/BPH_ancova.csv")

lm_out %>% 
  tidy(taggingr) %>% 
  write_csv("output/BPH_ancova_R2.csv")

lm_out %>% 
  augment(taggingr) %>% 
  mutate(fit = (.fitted)) %>% 
  ggplot(aes(x = radcap, y = lencap)) +
  geom_point()+facet_wrap(~zone)+geom_line(aes(radcap, fit), color ="grey50") + 
  annotate("text", x = 4, y=250, label="Body Proportional Hypothesis", family="Times New Roman")+
  scale_y_continuous(breaks = c(100, 150, 200, 250), limits = c(100,250))+
  labs(x = "Scale Radius (mm)", y =  "Capture Length (mm)")-> BPH3
ggsave("figs/BPH_ANCOVA.png", dpi = 500, height = 6, width = 8, units = "in")

#Linear models by zone (BPH)
lm_data_zone %>% 
  do(A1 = lm(lencap ~ radcap, data = sample1[sample1$zone=="A1",]),
     A2 = lm(lencap ~ radcap, data = sample1[sample1$zone=="A2",]),
     A3 = lm(lencap ~ radcap, data = sample1[sample1$zone=="A3",]),
     A4 = lm(lencap ~ radcap, data = sample1[sample1$zone=="A4",]),
     A6 = lm(lencap ~ radcap, data = sample1[sample1$zone=="A6",])) -> lm_out

lm_out %>% 
  tidy(A1) %>% 
  mutate(zone = "A1") -> A1
lm_out %>% 
  tidy(A2) %>% 
  mutate(zone = "A2") -> A2
lm_out %>% 
  tidy(A3) %>% 
  mutate(zone = "A3") -> A3
lm_out %>% 
  tidy(A4) %>% 
  mutate(zone = "A4") -> A4
lm_out %>% 
  tidy(A6) %>% 
  mutate(zone = "A6") -> A6

x<- rbind(A1, A2) #combine data for all zones
x<- rbind(x, A3)
x<- rbind(x, A4)
x<- rbind(x, A6)
write_csv(x, "output/BPH_lm.csv")

lm_out %>% 
  glance(A1) %>% 
  mutate(zone = "A1") -> A1
lm_out %>% 
  glance(A2) %>% 
  mutate(zone = "A2") -> A2
lm_out %>% 
  glance(A3) %>% 
  mutate(zone = "A3") -> A3
lm_out %>% 
  glance(A4) %>% 
  mutate(zone = "A4") -> A4
lm_out %>% 
  glance(A6) %>% 
  mutate(zone = "A6") -> A6

x<- rbind(A1, A2) #combine data for all zones
x<- rbind(x, A3)
x<- rbind(x, A4)
x<- rbind(x, A6)
write_csv(x, "output/BPH_lm_R2.csv")

lm_out %>% 
  augment(A1) %>% 
  mutate(fit = (.fitted)) %>% 
  ggplot(aes(x = radcap, y = lencap)) +
  geom_point(color ="grey50")+geom_line(aes(radcap, fit), color = "black") + 
  annotate("text", x = 4.5, y=250, label="BPH (A1)", family="Times New Roman")+
  scale_y_continuous(breaks = c(150, 200, 250), limits = c(150,250))+
  labs(x = "Scale Radius (mm)", y =  "Capture Length (mm)")-> A1

lm_out %>% 
  augment(A2) %>% 
  mutate(fit = (.fitted)) %>% 
  ggplot(aes(x = radcap, y = lencap)) +
  geom_point(color ="grey50")+geom_line(aes(radcap, fit), color = "black") + 
  annotate("text", x = 4.5, y=250, label="BPH (A2)", family="Times New Roman")+
  scale_y_continuous(breaks = c(150, 200, 250), limits = c(150,250))+
  labs(x = "Scale Radius (mm)", y =  "Capture Length (mm)")-> A2

lm_out %>% 
  augment(A3) %>% 
  mutate(fit = (.fitted)) %>% 
  ggplot(aes(x = radcap, y = lencap)) +
  geom_point(color ="grey50")+geom_line(aes(radcap, fit), color = "black") + 
  annotate("text", x = 4.75, y=250, label="BPH (A3)", family="Times New Roman")+
  scale_y_continuous(breaks = c(150, 200, 250), limits = c(150,250))+
  labs(x = "Scale Radius (mm)", y =  "Capture Length (mm)")-> A3

lm_out %>% 
  augment(A4) %>% 
  mutate(fit = (.fitted)) %>% 
  ggplot(aes(x = radcap, y = lencap)) +
  geom_point(color ="grey50")+geom_line(aes(radcap, fit), color = "black") + 
  annotate("text", x = 3.5, y=250, label="BPH (A4)", family="Times New Roman")+
  scale_y_continuous(breaks = c(150, 200, 250), limits = c(150,250))+
  labs(x = "Scale Radius (mm)", y =  "Capture Length (mm)")-> A4

lm_out %>% 
  augment(A6) %>% 
  mutate(fit = (.fitted)) %>% 
  ggplot(aes(x = radcap, y = lencap)) +
  geom_point(color ="grey50")+geom_line(aes(radcap, fit), color = "black") + 
  annotate("text", x = 3.5, y=250, label="BPH (A6)", family="Times New Roman")+
  scale_y_continuous(breaks = c(150, 200, 250), limits = c(150,250))+
  labs(x = "Scale Radius (mm)", y =  "Capture Length (mm)")-> A6
cowplot::plot_grid(A1, A2, A3, A4, A6,  align = "vh", nrow = 2, ncol=3)
ggsave("figs/BPH_regression.png", dpi = 500, height = 6, width =8, units = "in")

#Test #2: PCA 
#test for normality
eda.norm <- function(x, ...)
{
  par(mfrow=c(2,2))
  if(sum(is.na(x)) > 0)
    warning("NA's were removed before plotting")
  x <- x[!is.na(x)]
  hist(x, main = "Histogram and non-\nparametric density estimate", prob = T)
  iqd <- summary(x)[5] - summary(x)[2]
  lines(density(x, width = 2 * iqd))
  boxplot(x, main = "Boxplot", ...)
  qqnorm(x)
  qqline(x)
  plot.ecdf(x, main="Empirical and normal cdf")
  LIM <- par("usr")
  y <- seq(LIM[1],LIM[2],length=100)
  lines(y, pnorm(y, mean(x), sqrt(var(x))))
  shapiro.test(x)
}
data %>% 
  dplyr::select(id, fish, agecap, lencap, radcap, aprop1, aprop2, aprop3, aprop4, aprop5, aprop6, aprop7, aprop8, aprop9, zone) -> dataR 

# sample one from each fish and zone combo since have three samples from each fish/zone combo
set.seed(167) #set seed keeps random sample the same
dataR %>% 
  group_by(fish, zone) %>% 
  sample_n(1) %>%
  filterD(!is.na(aprop1)) %>% 
  as.data.frame() -> sample1

#correlation matrix of data
sample1 %>% dplyr::select(fish, aprop1, zone) %>%
  spread(key = zone, value = aprop1) -> data_wide
data <- data_wide[,2:length(data_wide)]
round((cor(data, use = "complete.obs")),2)


#Test #3: ANOVA with repeated treatments (between and within group variability) with multilevels 
#(a regression that allows for the errors to be dependent on eachother (as our conditions of Valence were repeated within each participant). 
#https://sapa-project.org/blog/2013/06/28/repeated-measures-anova-in-r/
data %>% 
  dplyr::select(id, fish, agecap, lencap, radcap, aprop1, aprop2, aprop3, aprop4, aprop5, aprop6, aprop7, aprop8, aprop9, zone) -> dataR 

# sample one from each fish and zone combo since have three samples from each fish/zone combo
set.seed(167) #set seed keeps random sample the same
dataR %>% 
  group_by(fish, zone) %>% 
  sample_n(1) %>%
  as.data.frame() -> sample1

#analysis of variance
a1 <- aov(sample1$aprop1 ~ sample1$zone)
summary(a1)
posthoc <- TukeyHSD(x=a1, 'sample1$zone', conf.level=0.95)
posthoc

#analysis of variance
a2 <- aov(sample1$aprop2 ~ sample1$zone)
summary(a2)
posthoc <- TukeyHSD(x=a2, 'sample1$zone', conf.level=0.95)
posthoc

#analysis of variance
a3 <- aov(sample1$aprop3 ~ sample1$zone)
summary(a3)
posthoc <- TukeyHSD(x=a3, 'sample1$zone', conf.level=0.95)
posthoc

#analysis of variance
a4 <- aov(sample1$aprop4 ~ sample1$zone)
summary(a4)
posthoc <- TukeyHSD(x=a4, 'sample1$zone', conf.level=0.95)
posthoc

#analysis of variance
a5 <- aov(sample1$aprop5 ~ sample1$zone)
summary(a5)
posthoc <- TukeyHSD(x=a5, 'sample1$zone', conf.level=0.95)
posthoc

#analysis of variance
a6 <- aov(sample1$aprop2 ~ sample1$zone)
summary(a8)
posthoc <- TukeyHSD(x=a8, 'sample1$zone', conf.level=0.95)
posthoc


#random effects model in r
fish <-as.factor(fish)
model <- lme(aprop1 ~ zone, random=~1|fish, data=sample1, method="REML") #model with random effects of fish
model.fixed <-  aov(aprop1 ~ zone, data=sample1)

anova(model.fixed) #compare models
posthoc <- glht(model.fixed, linfct = mcp(zone = "Tukey"))
summary(posthoc)
ref <- lsmeans(model,specs = c("zone"))






#Test #4: One-way repeated measures ANOVA for trimmed means
#https://cran.r-project.org/web/packages/WRS2/WRS2.pdf
eda.norm(sample1$aprop1)
sample1 %>% dplyr::select(fish, aprop1, zone) -> sample1
  spread(key = zone, value = aprop1) -> data_wide

x <- sample1[complete.cases(sample1), ]
rmanova(y=x$aprop1, groups=x$zone, x$fish)
rmmcp(sample1$aprop1, sample1$zone, sample1$fish)

rmanova(sample1$aprop2, sample1$zone, sample1$fish)
rmmcp(sample1$aprop2, sample1$zone, sample1$fish)

rmanova(sample1$aprop1, sample1$zone, sample1$fish)
rmmcp(sample1$aprop1, sample1$zone, sample1$fish)

