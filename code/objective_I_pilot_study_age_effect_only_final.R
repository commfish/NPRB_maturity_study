# notes---
# author: Sara E Miller 
# contact: sara.miller@alaska.gov; 907-465-4245
# Last edited: July, 2020
# Much of the regression code was orginally created by Ben Williams, and then adapted by Sara Miller.

# load libraries----
# devtools::install_github("commfish/fngr")
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
# calculate e, f, g for each zone
sample1 %>%
  filter(zone == "A1" & agei == 1)-> lm_data
sample1 %>%
  filter(agei == 1) %>% 
  as.data.frame() -> lm_data_zone

lm.age <- lm(radcap~lencap+agecap,data=lm_data)
e <- coef(lm.age)[[1]] 
f <- coef(lm.age)[[2]]
g <- coef(lm.age)[[3]]
h<-e/f
i<-g/f

sample1 %>%
    filter(zone == "A1") %>%
    mutate(SPH_age=((-h)+(lencap+h+i*agecap)*(radi/radcap)-(i*agei)))-> sample_A1 #age effect (Morita & Matsuishi 2001)

sample1 %>%
  filter(zone == "A2" & agei == 1)-> lm_data
lm.age <- lm(radcap~lencap+agecap,data=lm_data)
e <- coef(lm.age)[[1]] 
f <- coef(lm.age)[[2]]
g <- coef(lm.age)[[3]]
h<-e/f
i<-g/f

sample1 %>%
  filter(zone == "A2") %>%
  mutate(SPH_age=((-h)+(lencap+h+i*agecap)*(radi/radcap)-(i*agei)))-> sample_A2 #age effect (Morita & Matsuishi 2001)

sample1 %>%
  filter(zone == "A3" & agei == 1)-> lm_data #filter by agei so only one sample/region/fish
lm.age <- lm(radcap~lencap+agecap,data=lm_data)
e <- coef(lm.age)[[1]] 
f <- coef(lm.age)[[2]]
g <- coef(lm.age)[[3]]
h<-e/f
i<-g/f

sample1 %>%
  filter(zone == "A3") %>%
  mutate(SPH_age=((-h)+(lencap+h+i*agecap)*(radi/radcap)-(i*agei)))-> sample_A3 #age effect (Morita & Matsuishi 2001)

sample1 %>%
  filter(zone == "A4" & agei == 1)-> lm_data
lm.age <- lm(radcap~lencap+agecap,data=lm_data)
e <- coef(lm.age)[[1]] 
f <- coef(lm.age)[[2]]
g <- coef(lm.age)[[3]]
h<-e/f
i<-g/f

sample1 %>%
  filter(zone == "A4") %>%
  mutate(SPH_age=((-h)+(lencap+h+i*agecap)*(radi/radcap)-(i*agei)))-> sample_A4 #age effect (Morita & Matsuishi 2001)

sample1 %>%
  filter(zone == "A6" & agei == 1)-> lm_data
lm.age <- lm(radcap~lencap+agecap,data=lm_data)
e <- coef(lm.age)[[1]] 
f <- coef(lm.age)[[2]]
g <- coef(lm.age)[[3]]
h<-e/f
i<-g/f

sample1 %>%
  filter(zone == "A6") %>%
  mutate(SPH_age=((-h)+(lencap+h+i*agecap)*(radi/radcap)-(i*agei)))-> sample_A6 #age effect (Morita & Matsuishi 2001)

x<- rbind(sample_A1, sample_A2) #combine data for all zones
x<- rbind(x, sample_A3)
x<- rbind(x, sample_A4)
x<- rbind(x, sample_A6)
write.csv(x, "output/test.csv") # separate constants used for each region

# determine mean radi by region and annulus
x %>% 
  dplyr::select(fish, agei, zone, lencap, agecap, radi) %>%
  spread(key = zone, value = radi) %>%
  as.data.frame() -> data_wide_sph_age_radi # this does same as code above with merge
write.csv(data_wide_sph_age_radi, "output/data_wide_sph_age_radi.csv") 

alpha = 0.05
data_wide_sph_age_radi %>% 
  dplyr::select(fish, agei,A1, A2, A3, A4, A6) %>%  
  gather(variable, value, -agei, -fish) %>% 
  group_by(agei, variable) %>% 
  summarise(n = n(),
            mean_SPH_age=mean(value, na.rm=T),
            sd_SPH_age=sd(value, na.rm=T),
            n_SPH_age=n(),
            cv_SPH_age = sd_SPH_age / mean_SPH_age,
            se_SPH_age = sd(value, na.rm=T)/sqrt(n()),
            t = qt((1-alpha)/2 + .5, df = n-1), #https://www.r-graph-gallery.com/4-barplot-with-error-bar.html
            ci = t * se_SPH_age,
            upper_SPH_age = mean_SPH_age + ci,
            lower_SPH_age = mean_SPH_age - ci) %>% 
  mutate(age = as.factor(agei),
         zone=as.factor(variable)) -> data_wide_sph_age_radi

# create back calculated lengths by region
x %>% 
  dplyr::select(fish, agei, zone, SPH_age,lencap, agecap) %>%
  spread(key = zone, value = SPH_age) %>%
  as.data.frame() -> data_wide_sph_age 
write.csv(data_wide_sph_age, "output/data_wide_sph_age.csv") 

# figure 4 output found in data_wide_sph_age

# calculate difference in back-calculated zones in percents from zone A1
# source: http://www2.phy.ilstu.edu/~wenning/slh/Percent%20Difference%20Error.pdf

# SCALE PROPORTIONAL HYPOTHESIS
# AGE EFFECT SPH
# delete final annulus as it is the same for all regions
target <- c(11,22,33,44,55,66,77,88,99)
cols_to_concat <- c("agei", "agecap")
data_wide_sph_age %>% 
  dplyr::select(fish, agei, agecap, lencap, A1, A2, A3, A4, A6) %>%
  mutate(new_col = do.call(paste0, .[cols_to_concat])) %>%
  filter(!new_col %in% target) -> data_wide_sph_age
write.csv(data_wide_sph_age, "output/data_wide_sph_age_test.csv") 

alpha = 0.05
data_wide_sph_age %>% 
  mutate(A2= abs((A1-A2)/A1)*100,
         A3= abs((A1-A3)/A1)*100,
         A4= abs((A1-A4)/A1)*100,
         A6= abs((A1-A6)/A1)*100) %>% 
  dplyr::select(fish, agei,A2, A3, A4, A6) %>%  
  gather(variable, value, -agei, -fish) %>% 
  group_by(agei, variable) %>% 
  summarise(n = n(),
            mean_SPH_age=mean(value, na.rm=T),
            sd_SPH_age=sd(value, na.rm=T),
            n_SPH_age=n(),
            cv_SPH_age = sd_SPH_age / mean_SPH_age,
            se_SPH_age = sd(value, na.rm=T)/sqrt(n()),
            t = qt((1-alpha)/2 + .5, df = n-1), #https://www.r-graph-gallery.com/4-barplot-with-error-bar.html
            ci = t * se_SPH_age,
            upper_SPH_age = mean_SPH_age + ci,
            lower_SPH_age = mean_SPH_age - ci) %>% 
  mutate(age = as.factor(agei),
         zone=as.factor(variable)) -> data_wide_sph_age 
write.csv(data_wide_sph_age, "output/data_wide_sph_age_ouput.csv") 
# plots by age effect SPH
data_wide_sph_age %>% 
  group_by(age) %>% 
  summarise(labels = mean(n_SPH_age)) -> labels

data_wide_sph_age %>% 
  ggplot(aes(age, mean_SPH_age)) +
  geom_bar(aes(fill=zone), stat="identity", position="dodge", alpha=0.9) +
  scale_y_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60), limits = c(0,60))+
  scale_fill_grey(start = 0, end = .8, name="")+
  theme(legend.position=c(.9,.80)) +
  geom_errorbar(aes(ymin = lower_SPH_age, ymax = upper_SPH_age, group=zone),
                width = 0.2,
                linetype = "solid",
                position = position_dodge(width = 1),
                color="black", size = 0.75)+
  #annotate("text", x = 0.65, y=60, label="B)", family="Times New Roman") +
  geom_text(data = labels, aes(age, y = 60, label=paste ("n = " ,labels, sep =""), group=age),
            size = 3) +
  xlab("Age (annulus)") +
  ylab("Mean Percent Error") -> SPH_age
cowplot::plot_grid(SPH_age, align = "hv", nrow = 1, ncol=1) 
ggsave("figs/length_diff.png", dpi = 500, height = 4, width = 8, units = "in")

# radial figure
data_wide_sph_age_radi %>% 
  group_by(age) %>% 
  summarise(labels = mean(n_SPH_age)) -> labels

data_wide_sph_age_radi %>% 
  ggplot(aes(age, mean_SPH_age)) +
  geom_bar(aes(fill=zone), stat="identity", position="dodge", alpha=0.9) +
  scale_y_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6,7,8), limits = c(0,8))+
  scale_fill_grey(start = 0, end = .8,name ="") + 
  theme(legend.position=c(.08,.76)) +
  geom_errorbar(aes(ymin = lower_SPH_age, ymax = upper_SPH_age, group=zone),
                width = 0.2,
                linetype = "solid",
                position = position_dodge(width = 1),
                color="black", size = 0.75)+
  #annotate("text", x = 0.65, y=6, label="A)", family="Times New Roman") +
  geom_text(data = labels, aes(age, y = 8, label=paste ("n = " ,labels, sep =""), group=age),
            size = 3) +
  xlab("Age (annulus)") +
  ylab("Radial measurement (mm)") -> SPH_age2
cowplot::plot_grid(SPH_age2, align = "hv", nrow = 1, ncol=1) 
ggsave("figs/radi_diff.png", dpi = 500, height = 4, width = 8, units = "in")

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

# compare back calculated lengths for each region against region A1
#x %>%
#  filter(zone == "A1") %>%
#  dplyr::select(fish, agecap, lencap, radcap, agei, radi, zone, SPH_age) %>%
#  rename(zone_A1 = zone) %>%
#  rename(SPH_age_A1 = SPH_age) -> region_A1

#x %>%
#  dplyr::select(fish, agecap, lencap, radcap, agei, radi, zone, SPH_age) %>%
#  filter(zone == "A2") -> region_A2
#  
#merge(region_A1, region_A2, by = c('fish', 'agei', 'agecap', 'lencap')) -> A2_compare


#x %>%
#  dplyr::select(fish, agecap, lencap, radcap, agei, radi, zone, SPH_age) %>%
#  filter(zone == "A3") -> region_A3

#merge(region_A1, region_A3, by = c('fish', 'agei', 'agecap', 'lencap')) -> A3_compare

#x %>%
#  dplyr::select(fish, agecap, lencap, radcap, agei, radi, zone, SPH_age) %>%
#  filter(zone == "A4") -> region_A4

#merge(region_A1, region_A4, by = c('fish', 'agei', 'agecap', 'lencap')) -> A4_compare

#x %>%
#  dplyr::select(fish, agecap, lencap, radcap, agei, radi, zone, SPH_age) %>%
#  filter(zone == "A6") -> region_A6

#merge(region_A1, region_A6, by = c('fish', 'agei', 'agecap', 'lencap')) -> A6_compare  

#x %>%
#  group_by(agei, zone) %>%
#  summarize(n = n(),
#             mn_radcap = mean(radcap),
#             sd_radcap = sd(radcap),
#             ss_radcap = sum(radcap^2),
#             cv_radcap = sd_radcap / mn_radcap,
#             mn_SPH_age = mean(SPH_age),
#             sd_SPH_age = sd(SPH_age), 
#             ss_SPH_age = sum(SPH_age^2),
#             cv_SPH_age = sd_SPH_age / mn_SPH_age,
#             error_SPH_age = qt(0.975,df = n - 1) * sd_SPH_age / sqrt(n),
#             upper_SPH_age = mn_SPH_age + error_SPH_age,
#            lower_SPH_age = mn_SPH_age - error_SPH_age) %>%
#  as.data.frame() -> sample2
#write.csv(sample2, "output/table_summary_obj1.csv") # summary output for table 2