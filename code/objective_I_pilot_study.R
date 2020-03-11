#notes---
#author: Sara E Miller 
#contact: sara.miller@alaska.gov; 907-465-4245
#Last edited: February, 2019

# load libraries----
 devtools::install_github("commfish/fngr")
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
  filter(radcap >0) %>% #delete samples without measurements 
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
         prop9 = rad9 / radcap,
         aprop1 = asintrans(prop1), #transform data by arctransform
         aprop2 = asintrans(prop2),
         aprop3 = asintrans(prop3),
         aprop4 = asintrans(prop4),
         aprop5 = asintrans(prop5),
         aprop6 = asintrans(prop6),
         aprop7 = asintrans(prop7),
         aprop8 = asintrans(prop8),
         aprop9 = asintrans(prop9)) -> data 

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

ggplot(data = data, aes(x = lencap, y = aprop1, colour = specimen)) +
  geom_point() +
  facet_grid(~zone) +
  labs(x = "capture length", y = 
         "Transformed prop (focus to first annulus)")-> x

ggplot(data = data, aes(x = lencap, y = aprop1, colour = zone)) +
  geom_point() +
  facet_grid(~ zone) +
  labs(x = "capture", y = 
         "Transformed prop (focus to first annulus)")-> y


# Test #1: Back Calculation of Length by BPH and SPH Methods Compared to Preferred Area
# source: http://derekogle.com/IFAR/supplements/backcalculation/
set.seed(167) 

data %>% 
  dplyr::select(id, fish, agecap, lencap, contains('anu'), 
                contains('rad'), zone) %>% 
  group_by(fish, zone) %>% 
  sample_n(1)  %>%  # sample one scale from each fish zone (5 scales sampled from each fish)
  gather(agei, radi, rad1:rad9) %>%
  arrange(id, agei) %>% 
  dplyr::select(id, fish, agecap, lencap, radcap, agei, radi, zone) %>% 
  mutate(agei = as.numeric(stringr::str_sub(agei, 4, 4))) %>% 
  filter(!is.na(radi), agei<=agecap) -> sample1

# calculate starting values for back-calculation methods based on just zone A1 "preferred area"
# back-calculation methods with ratios
# Francis 1990 pg. 897 recommends the SPH and BPH methods; the difference btw the back-calculated lengths be taken as a 
# minimum meaure of imprecision of back-calculation
# calculate a, b, c, d for each zone

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
lm.age <- lm(radcap~lencap+agecap,data=lm_data)
e <- coef(lm.age)[[1]] 
f <- coef(lm.age)[[2]]
g <- coef(lm.age)[[3]]
h<-e/f
i<-g/f

sample1 %>%
    filter(zone == "A1") %>%
    mutate(SPH_len=(-a/b)+(lencap+a/b)*(radi/radcap), #Scale Proportional Hypothesis (Hile 1941:212)
                    BPH_len=lencap*((c+d*radi)/(c+d*radcap)),#Body Proportional Hypothesis (Whitney and Carlander 1956) 
                          SPH_age=((-h)+(lencap+h+i*agecap)*(radi/radcap)-(i*agei)))-> sample_A1 #age effect (Morita & Matsuishi 2001)


sample1 %>%
  filter(zone == "A2" & agei == 1)-> lm_data
lm.sl <- lm(radcap~lencap,data=lm_data)
a <- coef(lm.sl)[[1]] 
b <- coef(lm.sl)[[2]] 
lm.ls <- lm(lencap~radcap,data=lm_data)
c <- coef(lm.ls)[[1]] 
d <- coef(lm.ls)[[2]] 
lm.age <- lm(radcap~lencap+agecap,data=lm_data)
e <- coef(lm.age)[[1]] 
f <- coef(lm.age)[[2]]
g <- coef(lm.age)[[3]]

sample1 %>%
  filter(zone == "A2") %>%
  mutate(SPH_len=(-a/b)+(lencap+a/b)*(radi/radcap), #Scale Proportional Hypothesis (Hile 1941:212)
         BPH_len=lencap*((c+d*radi)/(c+d*radcap)),#Body Proportional Hypothesis (Whitney and Carlander 1956) 
         SPH_age=((-h)+(lencap+h+i*agecap)*(radi/radcap)-(i*agei)))-> sample_A2 #age effect (Morita & Matsuishi 2001)


sample1 %>%
  filter(zone == "A3" & agei == 1)-> lm_data #filter by agei so only one sample/region/fish
lm.sl <- lm(radcap~lencap,data=lm_data)
a <- coef(lm.sl)[[1]] 
b <- coef(lm.sl)[[2]] 
lm.ls <- lm(lencap~radcap,data=lm_data)
c <- coef(lm.ls)[[1]] 
d <- coef(lm.ls)[[2]] 
lm.age <- lm(radcap~lencap+agecap,data=lm_data)
e <- coef(lm.age)[[1]] 
f <- coef(lm.age)[[2]]
g <- coef(lm.age)[[3]]

sample1 %>%
  filter(zone == "A3") %>%
  mutate(SPH_len=(-a/b)+(lencap+a/b)*(radi/radcap), #Scale Proportional Hypothesis (Hile 1941:212)
         BPH_len=lencap*((c+d*radi)/(c+d*radcap)),#Body Proportional Hypothesis (Whitney and Carlander 1956) 
         SPH_age=((-h)+(lencap+h+i*agecap)*(radi/radcap)-(i*agei)))-> sample_A3 #age effect (Morita & Matsuishi 2001)


sample1 %>%
  filter(zone == "A4" & agei == 1)-> lm_data
lm.sl <- lm(radcap~lencap,data=lm_data)
a <- coef(lm.sl)[[1]] 
b <- coef(lm.sl)[[2]] 
lm.ls <- lm(lencap~radcap,data=lm_data)
c <- coef(lm.ls)[[1]] 
d <- coef(lm.ls)[[2]] 
lm.age <- lm(radcap~lencap+agecap,data=lm_data)
e <- coef(lm.age)[[1]] 
f <- coef(lm.age)[[2]]
g <- coef(lm.age)[[3]]

sample1 %>%
  filter(zone == "A4") %>%
  mutate(SPH_len=(-a/b)+(lencap+a/b)*(radi/radcap), #Scale Proportional Hypothesis (Hile 1941:212)
         BPH_len=lencap*((c+d*radi)/(c+d*radcap)),#Body Proportional Hypothesis (Whitney and Carlander 1956) 
         SPH_age=((-h)+(lencap+h+i*agecap)*(radi/radcap)-(i*agei)))-> sample_A4 #age effect (Morita & Matsuishi 2001)



sample1 %>%
  filter(zone == "A6" & agei == 1)-> lm_data
lm.sl <- lm(radcap~lencap,data=lm_data)
a <- coef(lm.sl)[[1]] 
b <- coef(lm.sl)[[2]] 
lm.ls <- lm(lencap~radcap,data=lm_data)
c <- coef(lm.ls)[[1]] 
d <- coef(lm.ls)[[2]] 
lm.age <- lm(radcap~lencap+agecap,data=lm_data)
e <- coef(lm.age)[[1]] 
f <- coef(lm.age)[[2]]
g <- coef(lm.age)[[3]]

sample1 %>%
  filter(zone == "A6") %>%
  mutate(SPH_len=(-a/b)+(lencap+a/b)*(radi/radcap), #Scale Proportional Hypothesis (Hile 1941:212)
         BPH_len=lencap*((c+d*radi)/(c+d*radcap)),#Body Proportional Hypothesis (Whitney and Carlander 1956) 
         SPH_age=((-h)+(lencap+h+i*agecap)*(radi/radcap)-(i*agei)))-> sample_A6 #age effect (Morita & Matsuishi 2001)

x<- rbind(sample_A1, sample_A2) #combine data for all zones
x<- rbind(x, sample_A3)
x<- rbind(x, sample_A4)
x<- rbind(x, sample_A6)
write.csv(x, "output/test.csv")
x %>%
  group_by(agei, zone) %>%
  summarize(n = n(),
             mn_radcap = mean(radcap),
             sd_radcap = sd(radcap),
             ss_radcap = sum(radcap^2),
             cv_radcap = sd_radcap / mn_radcap,
             mn_SPH_len = mean(SPH_len),
             sd_SPH_len = sd(SPH_len),
             ss_SPH_len = sum(SPH_len^2),
             cv_SPH = sd_SPH_len / mn_SPH_len,
             error_SPH = qt(0.975, df= n - 1) * sd_SPH_len / sqrt(n),
             upper_SPH = mn_SPH_len + error_SPH,
             lower_SPH = mn_SPH_len - error_SPH,
             mn_BPH_len = mean(BPH_len),
             sd_BPH_len = sd(BPH_len), 
             ss_BPH_len = sum(BPH_len^2),
             cv_BPH = sd_BPH_len / mn_BPH_len,
             error_BPH = qt(0.975,df = n - 1) * sd_BPH_len / sqrt(n),
             upper_BPH = mn_BPH_len + error_BPH,
             lower_BPH = mn_BPH_len - error_BPH, 
             mn_SPH_age = mean(SPH_age),
            sd_SPH_age = sd(SPH_age), 
            ss_SPH_age = sum(SPH_age^2),
            cv_SPH_age = sd_SPH_age / mn_SPH_age,
            error_SPH_age = qt(0.975,df = n - 1) * sd_SPH_age / sqrt(n),
            upper_SPH_age = mn_SPH_age + error_SPH_age,
            lower_SPH_age = mn_SPH_age - error_SPH_age) %>%
  as.data.frame() -> sample2
write.csv(sample2, "output/table_summary_obj1.csv") # summary output for table 2

x %>% dplyr::select(fish, agei, zone, SPH_len) %>%
             spread(key = zone, value = SPH_len)  %>%
  as.data.frame() -> data_wide_sph
  
x %>% dplyr::select(fish, agei, zone, BPH_len) %>%
  spread(key = zone, value = BPH_len) %>%
  as.data.frame() -> data_wide_bph

x %>% dplyr::select(fish, agei, zone, SPH_age) %>%
  spread(key = zone, value = SPH_age) %>%
  as.data.frame() -> data_wide_sph_age
# figure 4 output found in data_wide_sph, data_wide_bph, and data_wide_sph_age

# calculate difference in back-calculated zones in percents from zone A1
# source: http://www2.phy.ilstu.edu/~wenning/slh/Percent%20Difference%20Error.pdf

# SCALE PROPORTIONAL HYPOTHESIS
data_wide_sph %>% 
  mutate(A2= abs((A1-A2)/A1)*100,
         A3= abs((A1-A3)/A1)*100,
         A4= abs((A1-A4)/A1)*100,
         A6= abs((A1-A6)/A1)*100)%>% 
  dplyr::select(fish, agei,A2, A3, A4, A6) %>% 
  gather(variable, value, -agei, -fish) %>% 
  group_by(agei, variable) %>% 
  summarise(mean_SPH=mean(value, na.rm=T),
            sd_SPH=sd(value, na.rm=T),
            n_SPH=n(),
            se_SPH=sd(value, na.rm=T)/sqrt(n()))%>%
  mutate(age = as.factor(agei),
         zone=as.factor(variable)) -> data_wide_sph

# BODY PROPORTIONAL HYPOTHESIS
data_wide_bph %>% 
  mutate(A2= abs((A1-A2)/A1)*100,
         A3= abs((A1-A3)/A1)*100,
         A4= abs((A1-A4)/A1)*100,
         A6= abs((A1-A6)/A1)*100) %>% 
  dplyr::select(fish, agei,A2, A3, A4, A6) %>%  
  gather(variable, value, -agei, -fish) %>% 
  group_by(agei, variable) %>% 
  summarise(mean_BPH=mean(value, na.rm=T),
            sd_BPH=sd(value, na.rm=T),
            n_BPH=n(),
            se_BPH=sd(value, na.rm=T)/sqrt(n())) %>% 
  mutate(age = as.factor(agei),
         zone=as.factor(variable)) -> data_wide_bph 

# AGE EFFECT SPH
data_wide_sph_age %>% 
  mutate(A2= abs((A1-A2)/A1)*100,
         A3= abs((A1-A3)/A1)*100,
         A4= abs((A1-A4)/A1)*100,
         A6= abs((A1-A6)/A1)*100) %>% 
  dplyr::select(fish, agei,A2, A3, A4, A6) %>%  
  gather(variable, value, -agei, -fish) %>% 
  group_by(agei, variable) %>% 
  summarise(mean_SPH_age=mean(value, na.rm=T),
            sd_SPH_age=sd(value, na.rm=T),
            n_SPH_age=n(),
            se_SPH_age=sd(value, na.rm=T)/sqrt(n())) %>% 
  mutate(age = as.factor(agei),
         zone=as.factor(variable)) -> data_wide_sph_age 

# plots by SPH
data_wide_sph %>% 
  group_by(age) %>% 
  summarise(labels = mean(n_SPH)) -> labels

data_wide_sph %>% 
  ggplot(aes(age, mean_SPH)) +
  geom_bar(aes(fill=zone), stat="identity", position="dodge", alpha=0.9) +
  scale_fill_grey(start = 0, end = .8, name = "") + 
  theme(legend.position=c(.9,.75)) +
  annotate("text", x = 2, y=60, label="A) Scale Proportional Hypothesis", 
           family="Times") +
  geom_text(data = labels, aes(age, y = -1.2, label=paste ("n = " ,labels, sep =""), group=age), 
            size = 2.5) +
  xlab("Age") +
  ylab("Mean % Difference") -> SPH

tickr_length <- data.frame(mean.SPH = 0:300)
axisb <- tickr(tickr_length, mean.SPH, 50)
ggplot(data = sample2, aes(x = as.factor(agei), y = mn_SPH_len)) +
  geom_bar(aes(fill=zone), stat="identity", position="dodge",alpha=0.9) +
  scale_fill_grey(start = 0, end = .8)+theme(legend.position=c(.05,.7), legend.title=element_blank())+
  annotate("text", x = 2, y=250, label="A) Scale Proportional Hypothesis", family="Times New Roman")+
  scale_y_continuous(breaks = c(0, 50, 100, 150, 200, 250), limits = c(0,250))+
  geom_errorbar(aes(ymin = lower_SPH, ymax = upper_SPH, group=zone),
                width = 0.2,
                linetype = "solid",
                position = position_dodge(width = 1),
                color="black", size=1)+
  labs(x = "Age", y =  "Back-Calculated Length (mm)")-> SPH2

# Linear models by zone (SPH)
lm_data_zone %>% 
   do(A1 = lm(radcap ~ lencap, data = lm_data_zone[lm_data_zone$zone=="A1",]),
      A2 = lm(radcap ~ lencap, data = lm_data_zone[lm_data_zone$zone=="A2",]),
      A3 = lm(radcap ~ lencap, data = lm_data_zone[lm_data_zone$zone=="A3",]),
      A4 = lm(radcap ~ lencap, data = lm_data_zone[lm_data_zone$zone=="A4",]),
      A6 = lm(radcap ~ lencap, data = lm_data_zone[lm_data_zone$zone=="A6",])) -> lm_out_SPH

lm_out_SPH %>% 
  tidy(A1) %>% 
  mutate(zone = "A1") -> A1
lm_out_SPH %>% 
  tidy(A2) %>% 
mutate(zone = "A2") -> A2
lm_out_SPH %>% 
  tidy(A3) %>% 
mutate(zone = "A3") -> A3
lm_out_SPH %>% 
  tidy(A4) %>% 
mutate(zone = "A4") -> A4
lm_out_SPH %>% 
  tidy(A6) %>% 
  mutate(zone = "A6") -> A6

x<- rbind(A1, A2) #combine data for all zones
x<- rbind(x, A3)
x<- rbind(x, A4)
x<- rbind(x, A6)
write_csv(x, "output/SPH_lm.csv")

lm_out_SPH %>% 
  glance(A1) %>% 
  mutate(zone = "A1") -> A1
lm_out_SPH %>% 
  glance(A2) %>% 
  mutate(zone = "A2") -> A2
lm_out_SPH %>% 
  glance(A3) %>% 
  mutate(zone = "A3") -> A3
lm_out_SPH %>% 
  glance(A4) %>% 
  mutate(zone = "A4") -> A4
lm_out_SPH %>% 
  glance(A6) %>% 
  mutate(zone = "A6") -> A6

x<- rbind(A1, A2) #combine data for all zones
x<- rbind(x, A3)
x<- rbind(x, A4)
x<- rbind(x, A6)
write_csv(x, "output/SPH_lm_R2.csv")

lm_out_SPH %>% 
  augment(A1) %>% 
  mutate(fit = (.fitted)) %>% 
  ggplot(aes(x = lencap, y = radcap)) +
  geom_point(color ="grey50")+geom_line(aes(lencap, fit), color = "black") + 
  annotate("text", x = 200, y=7, label="SPH (A1)", family="Times New Roman")+
  scale_x_continuous(breaks = c(150, 200, 250), limits = c(150,250))+
  labs(y = "Scale Radius (mm)", x =  "Capture Length (mm)")-> A1

lm_out_SPH %>% 
  augment(A2) %>% 
  mutate(fit = (.fitted)) %>% 
  ggplot(aes(x = lencap, y = radcap)) +
  geom_point(color ="grey50")+geom_line(aes(lencap, fit), color = "black") + 
  annotate("text", x = 200, y=7, label="SPH (A2)", family="Times New Roman")+
  scale_x_continuous(breaks = c(150, 200, 250), limits = c(150,250))+
  labs(y = "Scale Radius (mm)", x =  "Capture Length (mm)")-> A2

lm_out_SPH %>% 
  augment(A3) %>% 
  mutate(fit = (.fitted)) %>% 
  ggplot(aes(x = lencap, y = radcap)) +
  geom_point(color ="grey50")+geom_line(aes(lencap, fit), color = "black") + 
  annotate("text", x = 200, y=7, label="SPH (A3)", family="Times New Roman")+
  scale_x_continuous(breaks = c(150, 200, 250), limits = c(150,250))+
  labs(y = "Scale Radius (mm)", x =  "Capture Length (mm)")-> A3

lm_out_SPH %>% 
  augment(A4) %>% 
  mutate(fit = (.fitted)) %>% 
  ggplot(aes(x = lencap, y = radcap)) +
  geom_point(color ="grey50")+geom_line(aes(lencap, fit), color = "black") + 
  annotate("text", x = 200, y=7, label="SPH (A4)", family="Times New Roman")+
  scale_x_continuous(breaks = c(150, 200, 250), limits = c(150,250))+
  labs(y = "Scale Radius (mm)", x =  "Capture Length (mm)")-> A4

lm_out_SPH %>% 
  augment(A6) %>% 
  mutate(fit = (.fitted)) %>% 
  ggplot(aes(x = lencap, y = radcap)) +
  geom_point(color ="grey50")+geom_line(aes(lencap, fit), color = "black") + 
  annotate("text", x = 200, y=7, label="SPH (A6)", family="Times New Roman")+
  scale_x_continuous(breaks = c(150, 200, 250), limits = c(150,250))+
  labs(y = "Scale Radius (mm)", x =  "Capture Length (mm)")-> A6
cowplot::plot_grid(A1, A2, A3, A4, A6,  align = "vh", nrow = 2, ncol=3)
ggsave("figs/SPH_regression.png", dpi = 500, height = 6, width = 8, units = "in")

# plots by BPH
data_wide_bph %>% 
  group_by(age) %>% 
  summarise(labels = mean(n_BPH)) -> labels

data_wide_bph %>% 
  ggplot(aes(age, mean_BPH)) +
  geom_bar(aes(fill=zone), stat="identity", position="dodge", alpha=0.9) +
  scale_fill_grey(start = 0, end = .8, guide = F) + 
  theme(legend.position=c(.9,.75)) +
  annotate("text", x = 2, y=40, 
           label="B) Body Proportional Hypothesis", family="Times") +
  geom_text(data = labels, aes(age, y = -1.2, label=paste ("n = " ,labels, sep =""), group=age),
            size = 2.5) +
  xlab("Age") +
  ylab("Mean % Difference") -> BPH

tickr_length <- data.frame(mean.BPH = 0:350)
axisb <- tickr(tickr_length, mean.BPH, 50)
ggplot(data = sample2, aes(x = as.factor(agei), y = mn_BPH_len)) +
  geom_bar(aes(fill=zone), stat="identity", position="dodge",alpha=0.9) +
  scale_fill_grey(start = 0, end = .8)+theme(legend.position="none")+
  annotate("text", x = 2, y=250, label="B) Body Proportional Hypothesis", family="Times New Roman")+
  scale_y_continuous(breaks = c(0, 50, 100, 150, 200, 250), limits = c(0,250))+
  geom_errorbar(aes(ymin = lower_BPH, ymax = upper_BPH, group=zone),
                width = 0.2,
                linetype = "solid",
                position = position_dodge(width = 1),
                color="black", size=1)+
  labs(x = "Age", y =  "Back-Calculated Length (mm)")-> BPH2

# Linear models by zone (BPH)
lm_data_zone %>% 
  do(A1 = lm(lencap ~ radcap, data = lm_data_zone[lm_data_zone$zone=="A1",]),
     A2 = lm(lencap ~ radcap, data = lm_data_zone[lm_data_zone$zone=="A2",]),
     A3 = lm(lencap ~ radcap, data = lm_data_zone[lm_data_zone$zone=="A3",]),
     A4 = lm(lencap ~ radcap, data = lm_data_zone[lm_data_zone$zone=="A4",]),
     A6 = lm(lencap ~ radcap, data = lm_data_zone[lm_data_zone$zone=="A6",])) -> lm_out_BPH

lm_out_BPH %>% 
  tidy(A1) %>% 
  mutate(zone = "A1") -> A1
lm_out_BPH %>% 
  tidy(A2) %>% 
  mutate(zone = "A2") -> A2
lm_out_BPH %>% 
  tidy(A3) %>% 
  mutate(zone = "A3") -> A3
lm_out_BPH %>% 
  tidy(A4) %>% 
  mutate(zone = "A4") -> A4
lm_out_BPH %>% 
  tidy(A6) %>% 
  mutate(zone = "A6") -> A6

x<- rbind(A1, A2) #combine data for all zones
x<- rbind(x, A3)
x<- rbind(x, A4)
x<- rbind(x, A6)
write_csv(x, "output/BPH_lm.csv")

lm_out_BPH %>% 
  glance(A1) %>% 
  mutate(zone = "A1") -> A1
lm_out_BPH %>% 
  glance(A2) %>% 
  mutate(zone = "A2") -> A2
lm_out_BPH %>% 
  glance(A3) %>% 
  mutate(zone = "A3") -> A3
lm_out_BPH %>% 
  glance(A4) %>% 
  mutate(zone = "A4") -> A4
lm_out_BPH %>% 
  glance(A6) %>% 
  mutate(zone = "A6") -> A6

x<- rbind(A1, A2) #combine data for all zones
x<- rbind(x, A3)
x<- rbind(x, A4)
x<- rbind(x, A6)
write_csv(x, "output/BPH_lm_R2.csv")

lm_out_BPH %>% 
  augment(A1) %>% 
  mutate(fit = (.fitted)) %>% 
  ggplot(aes(x = radcap, y = lencap)) +
  geom_point(color ="grey50")+geom_line(aes(radcap, fit), color = "black") + 
  annotate("text", x = 4.5, y=250, label="BPH (A1)", family="Times New Roman")+
  scale_y_continuous(breaks = c(150, 200, 250), limits = c(150,250))+
  labs(x = "Scale Radius (mm)", y =  "Capture Length (mm)")-> A1

lm_out_BPH %>% 
  augment(A2) %>% 
  mutate(fit = (.fitted)) %>% 
  ggplot(aes(x = radcap, y = lencap)) +
  geom_point(color ="grey50")+geom_line(aes(radcap, fit), color = "black") + 
  annotate("text", x = 4.5, y=250, label="BPH (A2)", family="Times New Roman")+
  scale_y_continuous(breaks = c(150, 200, 250), limits = c(150,250))+
  labs(x = "Scale Radius (mm)", y =  "Capture Length (mm)")-> A2

lm_out_BPH %>% 
  augment(A3) %>% 
  mutate(fit = (.fitted)) %>% 
  ggplot(aes(x = radcap, y = lencap)) +
  geom_point(color ="grey50")+geom_line(aes(radcap, fit), color = "black") + 
  annotate("text", x = 4.75, y=250, label="BPH (A3)", family="Times New Roman")+
  scale_y_continuous(breaks = c(150, 200, 250), limits = c(150,250))+
  labs(x = "Scale Radius (mm)", y =  "Capture Length (mm)")-> A3

lm_out_BPH %>% 
  augment(A4) %>% 
  mutate(fit = (.fitted)) %>% 
  ggplot(aes(x = radcap, y = lencap)) +
  geom_point(color ="grey50")+geom_line(aes(radcap, fit), color = "black") + 
  annotate("text", x = 3.5, y=250, label="BPH (A4)", family="Times New Roman")+
  scale_y_continuous(breaks = c(150, 200, 250), limits = c(150,250))+
  labs(x = "Scale Radius (mm)", y =  "Capture Length (mm)")-> A4

lm_out_BPH %>% 
  augment(A6) %>% 
  mutate(fit = (.fitted)) %>% 
  ggplot(aes(x = radcap, y = lencap)) +
  geom_point(color ="grey50")+geom_line(aes(radcap, fit), color = "black") + 
  annotate("text", x = 3.5, y=250, label="BPH (A6)", family="Times New Roman")+
  scale_y_continuous(breaks = c(150, 200, 250), limits = c(150,250))+
  labs(x = "Scale Radius (mm)", y =  "Capture Length (mm)")-> A6
cowplot::plot_grid(A1, A2, A3, A4, A6,  align = "vh", nrow = 2, ncol=3)
ggsave("figs/BPH_regression.png", dpi = 500, height = 6, width =8, units = "in")

# plots by age effect SPH
data_wide_sph_age %>% 
  group_by(age) %>% 
  summarise(labels = mean(n_SPH_age)) -> labels

data_wide_sph_age %>% 
  ggplot(aes(age, mean_SPH_age)) +
  geom_bar(aes(fill=zone), stat="identity", position="dodge", alpha=0.9) +
  scale_fill_grey(start = 0, end = .8, guide = F) + 
  theme(legend.position=c(.9,.75)) +
  annotate("text", x = 2.5, y=10, 
           label="C) Age Effect Scale Proportional Hypothesis", family="Times") +
  geom_text(data = labels, aes(age, y = -1.2, label=paste ("n = " ,labels, sep =""), group=age),
            size = 2.5) +
  xlab("Age") +
  ylab("Mean % Difference") -> SPH_age

cowplot::plot_grid(SPH, BPH, SPH_age,   align = "hv", nrow = 3, ncol=1) 
ggsave("figs/length_diff.png", dpi = 500, height = 8, width = 8, units = "in")

tickr_length <- data.frame(mean.SPH.age = 0:350)
axisb <- tickr(tickr_length, mean.SPH.age, 50)
ggplot(data = sample2, aes(x = as.factor(agei), y = mn_SPH_age)) +
  geom_bar(aes(fill=zone), stat="identity", position="dodge",alpha=0.9) +
  scale_fill_grey(start = 0, end = .8)+theme(legend.position="none")+
  annotate("text", x = 2.5, y=250, label="C) Age Effect Scale Proportional Hypothesis", family="Times New Roman")+
  scale_y_continuous(breaks = c(0, 50, 100, 150, 200, 250), limits = c(0,250))+
  geom_errorbar(aes(ymin = lower_SPH_age, ymax = upper_SPH_age, group=zone),
                width = 0.2,
                linetype = "solid",
                position = position_dodge(width = 1),
                color="black", size=1)+
  labs(x = "Age", y =  "Back-Calculated Length (mm)")-> SPH_age2
cowplot::plot_grid(SPH2, BPH2,SPH_age2,   align = "hv", nrow = 3, ncol=1) 
ggsave("figs/length.png", dpi = 500, height = 8, width = 8, units = "in")

# Linear models by zone (SPH age)
lm_data_zone %>% 
  do(A1 = lm(radcap ~ lencap +agecap, data = lm_data_zone[lm_data_zone$zone=="A1",]),
     A2 = lm(radcap ~ lencap +agecap, data = lm_data_zone[lm_data_zone$zone=="A2",]),
     A3 = lm(radcap ~ lencap +agecap, data = lm_data_zone[lm_data_zone$zone=="A3",]),
     A4 = lm(radcap ~ lencap +agecap, data = lm_data_zone[lm_data_zone$zone=="A4",]),
     A6 = lm(radcap ~ lencap +agecap, data = lm_data_zone[lm_data_zone$zone=="A6",])) -> lm_out_SPH_age

lm_out_SPH_age %>% 
  tidy(A1) %>% 
  mutate(zone = "A1") -> A1
lm_out_SPH_age %>% 
  tidy(A2) %>% 
  mutate(zone = "A2") -> A2
lm_out_SPH_age %>% 
  tidy(A3) %>% 
  mutate(zone = "A3") -> A3
lm_out_SPH_age %>% 
  tidy(A4) %>% 
  mutate(zone = "A4") -> A4
lm_out_SPH_age %>% 
  tidy(A6) %>% 
  mutate(zone = "A6") -> A6

x<- rbind(A1, A2) #combine data for all zones
x<- rbind(x, A3)
x<- rbind(x, A4)
x<- rbind(x, A6)
write_csv(x, "output/SPH_age_lm.csv")

lm_out_SPH_age %>% 
  glance(A1) %>% 
  mutate(zone = "A1") -> A1
lm_out_SPH_age %>% 
  glance(A2) %>% 
  mutate(zone = "A2") -> A2
lm_out_SPH_age %>% 
  glance(A3) %>% 
  mutate(zone = "A3") -> A3
lm_out_SPH_age %>% 
  glance(A4) %>% 
  mutate(zone = "A4") -> A4
lm_out_SPH_age %>% 
  glance(A6) %>% 
  mutate(zone = "A6") -> A6

x<- rbind(A1, A2) #combine data for all zones
x<- rbind(x, A3)
x<- rbind(x, A4)
x<- rbind(x, A6)
write_csv(x, "output/SPH_age_lm_R2.csv")

# BPH and SPH regression fig
lm_out_SPH %>% 
  augment(A1) %>% 
  mutate(fit_SPH = (.fitted)) %>% 
  as.data.frame()-> lm_out_SPH1

lm_out_BPH %>% 
  augment(A1) %>% 
  mutate(fit_BPH = (.fitted)) %>% 
  as.data.frame()-> lm_out_BPH1

ggplot(aes(x = radcap, y = lencap), data=lm_out_BPH1) +
  geom_point(color ="grey50")+geom_line(aes(x=radcap, y=fit_BPH), color = "black", size = 1) + 
  geom_line(data =lm_out_SPH1, aes(x = fit_SPH, y = lencap), color = "black", size=2) + 
  geom_segment(aes(x = 5.32, y = 41, xend = 5.32, yend = 213), size=1, colour="grey80") + #Sc
  geom_segment(aes(x = 0, y = 213, xend = 5.32, yend = 213), size=1, colour="grey80") + #Lc
  geom_segment(aes(x = 0, y = 118.52, xend = 2.35, yend = 118.52), size=1, colour="grey80") + #SPH
  geom_segment(aes(x = 2.35, y = 41, xend = 2.35, yend = 118.52), size=1, colour="grey80") +     #SPH
  geom_segment(aes(x = 5.32, y = 213, xend = 0, yend = 43.77), size=2, colour="black", lty=2) +     #SPH
  geom_segment(aes(x = 0, y = 143.43, xend = 2.35, yend = 143.43), size=1, colour="grey80") +        #BPH
  geom_segment(aes(x = 2.35, y = 41, xend = 2.35, yend = 143.43), size=1, colour="grey80") +         #BPH
  geom_segment(aes(x = 5.32, y = 213, xend = 0, yend = 88.38), size=1, colour="black", lty=2)+  #BPH
  annotate("text", x = 0.3, y=146, label="L1 (BPH)", family="Times New Roman")+
  annotate("text", x = 0.3, y=122, label="L1 (SPH)", family="Times New Roman")+
  annotate("text", x = 5.32, y=217, label="P", family="Times New Roman")+
  annotate("text", x = 0.1, y=217, label="Lc", family="Times New Roman")+
  annotate("text", x = 5.32, y=40, label="Sc", family="Times New Roman")+
  annotate("text", x = 2.35, y=40, label="S1", family="Times New Roman")+
  scale_y_continuous(breaks = c(40, 80, 120, 160, 200, 240), limits = c(40,240))+
  #scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7), limits = c(0,7))+
  scale_x_continuous(expand = c(0, 0), limits=c(0,6.5))+
  labs(x = "Scale Radius (mm)", y =  "Capture Length (mm)")-> A1 
ggsave("figs/BPH_SPH_regression_fit.png", dpi = 500, height = 6, width =8, units = "in")

my.formula <- y ~ x
my.formula1 <- x ~ y
  ggplot(aes(x = radcap, y = lencap), data = lm_out_BPH1) +
    stat_poly_eq(formula = my.formula,
                 eq.with.lhs = "italic(hat(y))~`=`~",
                 aes(label = paste(..eq.label.., ..rr.label.., sep = "*plain(\",\")~")), 
                 parse = TRUE, size = 3, vjust=2) +
    stat_poly_eq(formula = my.formula1,
                 eq.with.lhs = "italic(hat(y))~`=`~",
                 aes(label = paste(..eq.label.., ..rr.label.., sep = "*plain(\",\")~")), 
                 parse = TRUE, vjust=3, size=3) +
  geom_point(color ="grey50")+geom_line(aes(x=radcap, y=fit_BPH), color = "black", lty=2, size = 1) + 
  geom_line(data =lm_out_SPH1, aes(x = fit_SPH, y = lencap), color = "black", size=1) + 
  annotate("text", x = 4.5, y=250, label="A1", size=5, family="Times New Roman")+
  scale_y_continuous(breaks = c(150, 200, 250), limits = c(150,250))+
  scale_x_continuous(breaks = c(2,3,4,5,6,7), limits = c(2,7))+
  labs(x = "Scale Radius (mm)", y =  "Capture Length (mm)")-> A1 

lm_out_SPH %>% 
  augment(A2) %>% 
  mutate(fit_SPH = (.fitted)) -> lm_out_SPH2
lm_out_BPH %>% 
  augment(A2) %>% 
  mutate(fit_BPH = (.fitted)) %>% 
  ggplot(aes(x = radcap, y = lencap)) +
  stat_poly_eq(formula = my.formula,
               eq.with.lhs = "italic(hat(y))~`=`~",
               aes(label = paste(..eq.label.., ..rr.label.., sep = "*plain(\",\")~")), 
               parse = TRUE, size = 3, vjust=2) +
  stat_poly_eq(formula = my.formula1,
               eq.with.lhs = "italic(hat(y))~`=`~",
               aes(label = paste(..eq.label.., ..rr.label.., sep = "*plain(\",\")~")), 
               parse = TRUE, vjust=3, size=3) +
  geom_point(color ="grey50")+geom_line(aes(x=radcap, y=fit_BPH), color = "black", lty=2, size = 1) + 
  geom_line(data =lm_out_SPH2, aes(x = fit_SPH, y = lencap), color = "black", size=1, lty=1) + 
  annotate("text", x = 4.5, y=250, label="A2", size=5,family="Times New Roman")+
  scale_y_continuous(breaks = c(150, 200, 250), limits = c(150,250))+
  scale_x_continuous(breaks = c(2,3,4,5,6,7), limits = c(2,7))+
  labs(x = "Scale Radius (mm)", y =  "Capture Length (mm)")-> A2 

lm_out_SPH %>% 
  augment(A3) %>% 
  mutate(fit_SPH = (.fitted)) -> lm_out_SPH3
lm_out_BPH %>% 
  augment(A3) %>% 
  mutate(fit_BPH = (.fitted)) %>% 
  ggplot(aes(x = radcap, y = lencap)) +
  stat_poly_eq(formula = my.formula,
               eq.with.lhs = "italic(hat(y))~`=`~",
               aes(label = paste(..eq.label.., ..rr.label.., sep = "*plain(\",\")~")), 
               parse = TRUE, size = 3, vjust=2) +
  stat_poly_eq(formula = my.formula1,
               eq.with.lhs = "italic(hat(y))~`=`~",
               aes(label = paste(..eq.label.., ..rr.label.., sep = "*plain(\",\")~")), 
               parse = TRUE, vjust=3, size=3) +
  geom_point(color ="grey50")+geom_line(aes(x=radcap, y=fit_BPH), color = "black", lty=2, size = 1) + 
  geom_line(data =lm_out_SPH3, aes(x = fit_SPH, y = lencap), color = "black", size=1, lty=1) + 
  annotate("text", x = 4.75, y=250, label="A3", size=5, family="Times New Roman")+
  scale_y_continuous(breaks = c(150, 200, 250), limits = c(150,250))+
  labs(x = "Scale Radius (mm)", y =  "Capture Length (mm)")-> A3 

lm_out_SPH %>%
  augment(A4) %>% 
  mutate(fit_SPH = (.fitted))  -> lm_out_SPH4
lm_out_BPH %>% 
  augment(A4) %>% 
  mutate(fit_BPH = (.fitted)) %>% 
  ggplot(aes(x = radcap, y = lencap)) +
  stat_poly_eq(formula = my.formula,
               eq.with.lhs = "italic(hat(y))~`=`~",
               aes(label = paste(..eq.label.., ..rr.label.., sep = "*plain(\",\")~")), 
               parse = TRUE, size = 3, vjust=2) +
  stat_poly_eq(formula = my.formula1,
               eq.with.lhs = "italic(hat(y))~`=`~",
               aes(label = paste(..eq.label.., ..rr.label.., sep = "*plain(\",\")~")), 
               parse = TRUE, vjust=3, size=3) +
  geom_point(color ="grey50")+geom_line(aes(x=radcap, y=fit_BPH), color = "black", size = 1, lty=2) + 
  geom_line(data =lm_out_SPH4, aes(x = fit_SPH, y = lencap), color = "black", size=1, lty=1) + 
  annotate("text", x = 3.75, y=250, label="A4", size=5, family="Times New Roman")+
  scale_y_continuous(breaks = c(150, 200, 250), limits = c(150,250))+
  labs(x = "Scale Radius (mm)", y =  "Capture Length (mm)")-> A4 

lm_out_SPH %>%  
  augment(A6) %>% 
  mutate(fit_SPH = (.fitted)) -> lm_out_SPH6
lm_out_BPH %>% 
  augment(A6) %>% 
  mutate(fit_BPH = (.fitted)) %>% 
  ggplot(aes(x = radcap, y = lencap)) +
  stat_poly_eq(formula = my.formula,
               eq.with.lhs = "italic(hat(y))~`=`~",
               aes(label = paste(..eq.label.., ..rr.label.., sep = "*plain(\",\")~")), 
               parse = TRUE, size = 3, vjust=2) +
  stat_poly_eq(formula = my.formula1,
               eq.with.lhs = "italic(hat(y))~`=`~",
               aes(label = paste(..eq.label.., ..rr.label.., sep = "*plain(\",\")~")), 
               parse = TRUE, vjust=3, size=3) +
  geom_point(color ="grey50")+geom_line(aes(x=radcap, y=fit_BPH), color = "black", size = 1, lty=2) + 
  geom_line(data =lm_out_SPH6, aes(x = fit_SPH, y = lencap), color = "black", size=1, lty=1) + 
  annotate("text", x = 4, y=250, label="A6", size=5, family="Times New Roman")+
  scale_y_continuous(breaks = c(150, 200, 250), limits = c(150,250))+
  labs(x = "Scale Radius (mm)", y =  "Capture Length (mm)")-> A6 
cowplot::plot_grid(A1, A2, A3, A4, A6,  align = "vh", nrow = 2, ncol=3)
ggsave("figs/BPH_SPH_regression.png", dpi = 500, height = 6, width =8, units = "in") 

# ANCOVA Models (SPH)
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

#Test #5: Maximum difference in BPH and SPH 
