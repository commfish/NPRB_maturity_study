#notes---
#author: Sara E Miller 
#contact: sara.miller@alaska.gov; 907-465-4245
#Last edited: July, 2019

# load libraries----
# devtools::install_github("ben-williams/FNGr")
source("code/functions.R")
library(extrafont)
loadfonts(device="win")
library(nlme)
library(tidyverse)
library(FNGr)
library(broom)
library(cowplot)
library(mixtools)
#library(mixdist)
#library(fitdistrplus)
#library(gam)
library(car)
#library(lmtest)
#library(mgcv)
#library(visreg)
#library(boot)
#library(AICcmodavg) #AICc
#library(rms)
library(modEvA) #pseudo-R squared
library(MASS)
library(digest)
theme_set(theme_sleek())

# load data ---- 
data <- read.csv("data/data_11_26_2018.csv", check.names = FALSE) 
increment <- read.csv("data/increment_measurements.csv", check.names = FALSE) 

#data clean ----
data %>% 
  filter(!(maturation_status_histology %in% c("", NA, "no slide", "no score", "4-1", "3-4", "2-3", "1-2")))-> data_clean

# check for duplicate image names in the increment file
increment %>% 
  group_by(image_name) %>% 
  filter(n()>1) %>% 
  summarize(n=n())-> x

# match increments to awl data
merge1<-merge(increment, data_clean, all.y=T)
merge1 %>% 
  filter(!(Increment1 %in% c("", NA))) %>% 
  filter(!(scale_region%in% c("OOA (OUT OF AREA)", "A", "C", "D", "E", "G", "H")))-> clean_dataset #n=296 useable samples

#sample sizes for study 
clean_dataset %>% 
  dplyr::select(age,maturity_state_histology, scale_region, sex_histology) %>% 
  group_by(maturity_state_histology, age, sex_histology) %>% 
  summarize(n=n())-> table1 

clean_dataset %>% 
  dplyr::select(maturity_state_histology, maturation_status_histology) %>% 
  group_by(maturity_state_histology,maturation_status_histology) %>% 
  summarize(n=n())-> table2 

clean_dataset %>% 
  dplyr::select(scale_region) %>% 
  group_by(scale_region) %>% 
  summarize(n=n())-> table3 

#calculate proportion of outer scale ring
clean_dataset %>% 
  mutate(mature = ifelse(maturity_state_histology == 'I', 'immature', 'mature'),
         anu1 = ifelse(is.na(Increment1), 0 , Increment1),
         anu2 = ifelse(is.na(Increment2), 0 , Increment2),
         anu3 = ifelse(is.na(Increment3), 0 , Increment3),
         anu4 = ifelse(is.na(Increment4), 0 , Increment4),
         anu5 = ifelse(is.na(Increment5), 0 , Increment5),
         anu6 = ifelse(is.na(Increment6), 0 , Increment6),
         anu7 = ifelse(is.na(Increment7), 0 , Increment7)) %>% 
  rowwise() %>% 
  mutate(radcap = sum(anu1, anu2, anu3, anu4, anu5, anu6, anu7)) %>% 
  mutate(rad1 = anu1,#calculate radial units (ie. focus to each annulus)
         rad2 = sum(anu1, anu2),
         rad3 = sum(anu1, anu2, anu3),
         rad4 = sum(anu1, anu2, anu3, anu4),
         rad5 = sum(anu1, anu2, anu3, anu4, anu5),
         rad6 = sum(anu1, anu2, anu3, anu4, anu5, anu6),
         rad7 = sum(anu1, anu2, anu3, anu4, anu5, anu6, anu7),
         prop1 = rad1 / radcap,
         prop2 = rad2 / radcap,
         prop3 = rad3 / radcap,
         prop4 = rad4 / radcap,
         prop5 = rad5 / radcap,
         prop6 = rad6 / radcap,
         prop7 = rad7 / radcap)%>%
  mutate(prop1_adj = prop1,
         prop2_adj = ifelse(prop2==1, NA, prop2),
         prop3_adj = ifelse(prop3==1, NA, prop3),
         prop4_adj = ifelse(prop4==1, NA, prop4),
         prop5_adj = ifelse(prop5==1, NA, prop5),
         prop6_adj = ifelse(prop6==1, NA, prop6),
         prop7_adj = ifelse(prop7==1, NA, prop7)) -> sample1
 
sample1 %>%
  gather(value, variable, prop1_adj:prop7_adj) %>% 
  group_by(image_name) %>% 
  filter(!is.na(variable)) %>%
  summarize(max=max(variable))%>%
  mutate(outer_prop=1-max)-> sample2 #outer prop

merge<- merge(x = sample1, y= sample2, by=c("image_name"), all.x  = T)

#anu_adj is measurement of outer ring
#aprop is arcsine transform of outer_prop
merge %>%
  mutate(aprop = asintrans(outer_prop),
         anu_adj = ifelse(anu1>0 & anu2==0 & anu3 ==0 & anu4 ==0 & anu5==0 & anu6==0 & anu7==0, anu1,
                          ifelse(anu1>0 & anu2>0& anu3 ==0 & anu4 ==0 & anu5==0 & anu6==0 & anu7==0, anu2,
                                 ifelse(anu1>0 & anu2>0 & anu3>0 & anu4 ==0 & anu5==0 & anu6==0 & anu7==0, anu3,
                                        ifelse(anu1>0 & anu2>0 & anu3>0 & anu4>0 & anu5==0 & anu6==0 & anu7==0, anu4,
                                               ifelse(anu1>0 & anu2>0 & anu3>0 & anu4>0 & anu5>0 & anu6==0 & anu7==0, anu5,
                                                      ifelse(anu1>0 & anu2>0 & anu3>0 & anu4>0 & anu5>0 & anu6>0 & anu7==0, anu6,anu7))))))) %>%
  mutate(maturity = ifelse(mature=='mature', 1, 0)) %>%
  dplyr::select(image_name,year, age, sex_histology, maturation_status_histology, mature, max, outer_prop, aprop, anu_adj, maturity) -> merge
merge %>% 
  mutate(age=as.factor(age),
         maturity = ifelse(mature == "mature", 1, 0)) -> merge

merge %>%
  filter(age == 3) -> dataset

#Exploratory Plots
#A) Histograms of outer ring 
merge %>%
  filter(age == 2 & mature == "mature")-> age2mature
merge %>%
  filter(age == 2 & mature == "immature")-> age2immature
merge %>%
  filter(age == 3 & mature == "mature")-> age3mature
merge %>%
  filter(age == 3 & mature == "immature")-> age3immature
merge %>%
  filter(age == 4 & mature == "mature")-> age4mature
merge %>%
  filter(age == 4 & mature == "immature")-> age4immature
merge %>%
  filter(age == 5 & mature == "mature")-> age5mature
merge %>%
  filter(age == 5 & mature == "immature")-> age5immature
merge %>%
  filter(age == 6 & mature == "mature")-> age6mature
merge %>%
  filter(age == 6 & mature == "immature")-> age6immature

#test each age/maturity for normality
eda.norm(as.numeric(age2immature$aprop))
#eda.norm(as.numeric(age2mature$aprop))
eda.norm(as.numeric(age3immature$aprop))
eda.norm(as.numeric(age3mature$aprop))
#eda.norm(as.numeric(age4immature$aprop))
eda.norm(as.numeric(age4mature$aprop))
eda.norm(as.numeric(age5mature$aprop)) #normal
eda.norm(as.numeric(age6mature$aprop)) #normal

merge %>%
  filter(age == 2) %>%
ggplot(., aes(x=outer_prop, color=mature, fill=mature)) + geom_histogram(alpha=0.5, position = 'identity') +
  ylab("Frequency")+ xlab("Outer increment proportion") +
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6), limits = c(0,0.6))+
  ggtitle("Age 2; n=126") + theme(legend.position="none") +
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_fill_manual(values=c("#999999", "#E69F00" )) -> plot1

merge %>%
  filter(age == 2) %>%
ggplot(., aes(x=outer_prop, color=mature, fill=mature)) +
  geom_density(alpha=0.5, adjust=1) +scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_fill_manual(values=c("#999999", "#E69F00" )) + theme(legend.title=element_blank(), legend.position=c(.85,.85)) + 
  ylab("Density")+ xlab("Outer increment proportion") +
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6), limits = c(0,0.6))-> plot2

merge %>%
  filter(age == 3) %>%
ggplot(., aes(x=outer_prop, color=mature, fill=mature)) + geom_histogram(alpha=0.5, position = 'identity') +
  ylab("Frequency")+ xlab("Outer increment proportion") +
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6), limits = c(0,0.6))+
  ggtitle("Age 3; n=72") + theme(legend.position="none") +
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_fill_manual(values=c("#999999", "#E69F00" )) -> plot3

merge %>%
  filter(age == 3) %>%
  ggplot(., aes(x=outer_prop, color=mature, fill=mature)) +
  geom_density(alpha=0.5, adjust=1) +scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_fill_manual(values=c("#999999", "#E69F00" )) + theme(legend.title=element_blank(), legend.position=c(.85,.85)) +
  ylab("Density")+ xlab("Outer increment proportion") +
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6), limits = c(0,0.6))-> plot4

merge %>%
  filter(age == 4) %>%
  ggplot(., aes(x=outer_prop, color=mature, fill=mature)) + geom_histogram(alpha=0.5, position = 'identity') +
  ylab("Frequency")+ xlab("Outer increment proportion") +
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6), limits = c(0,0.6)) +
  ggtitle("Age 4; n=41") + theme(legend.position="none") +
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_fill_manual(values=c("#999999", "#E69F00" )) -> plot5

merge %>%
  filter(age == 4) %>%
  ggplot(., aes(x=outer_prop, color=mature, fill=mature)) +
  geom_density(alpha=0.5, adjust=1) +scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_fill_manual(values=c("#999999", "#E69F00" )) + theme(legend.title=element_blank(), legend.position=c(.85,.85)) +
  ylab("Density")+ xlab("Outer increment proportion") +
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6), limits = c(0,0.6))-> plot6

merge %>%
  filter(age == 5) %>%
  ggplot(., aes(x=outer_prop, color=mature, fill=mature)) + geom_histogram(alpha=0.5, position = 'identity') +
  ylab("Frequency")+ xlab("Outer increment proportion") +
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6), limits = c(0,0.6)) +
  ggtitle("Age 5; n=24") + theme(legend.position="none") +
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_fill_manual(values=c("#999999", "#E69F00" )) -> plot7

merge %>%
  filter(age == 5) %>%
  ggplot(., aes(x=outer_prop, color=mature, fill=mature)) +
  geom_density(alpha=0.5, adjust=1) +scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_fill_manual(values=c("#999999", "#E69F00" )) + theme(legend.title=element_blank(), legend.position=c(.85,.85)) +
  ylab("Density")+ xlab("Outer increment proportion") +
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6), limits = c(0,0.6)) -> plot8

merge %>%
  filter(age == 6) %>%
  ggplot(., aes(x=outer_prop, color=mature, fill=mature)) + geom_histogram(alpha=0.5, position = 'identity') +
  ylab("Frequency")+ xlab("Outer increment proportion") +
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6), limits = c(0,0.6)) +
  ggtitle("Age 6; n=33") + theme(legend.position="none") +
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_fill_manual(values=c("#999999", "#E69F00" )) -> plot9

merge %>%
  filter(age == 6) %>%
  ggplot(., aes(x=outer_prop, color=mature, fill=mature)) +
  geom_density(alpha=0.5, adjust=1) +scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_fill_manual(values=c("#999999", "#E69F00" )) + theme(legend.title=element_blank(), legend.position=c(.85,.85)) +
  ylab("Density")+ xlab("Outer increment proportion") +
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6), limits = c(0,0.6))-> plot10
cowplot::plot_grid(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, plot9, plot10, align = "vh", nrow = 5, ncol=2)
ggsave("figs/histogram.png", dpi = 500, height = 10, width =8, units = "in")

#B) Fit gaussian mixture models
#datasets by ages
merge %>%
  filter(age == 2) -> age2
merge %>%
  filter(age == 3) -> age3
merge %>%
  filter(age == 4) -> age4
merge %>%
  filter(age == 5) -> age5
merge %>%
  filter(age == 6) -> age6  

png(file='figs/Mixture_Models.png', res=200, width=13, height=10, units ="in") 
par(mfrow=c(3,2)) 
age2<-age2[,c(8)]
mixmdl2 <- normalmixEM(age2)
x3<-plot(mixmdl2, which=2, cex.axis=1, cex.lab=1, cex.main=1.1,
         main2="Age 2", xlab2="Scale increment Proportions", xlim=c(0,1))
x3<-lines(density(age2), lty=2, lwd=2)
x3<-legend("topright", c("mature", "immature"), col = c(2,3), lty = c(1,1), lwd=2,
           merge = TRUE, bty = "n")
summary(mixmdl2)

age3<-age3[,c(8)]
mixmdl3 <- normalmixEM(age3)
x3<-plot(mixmdl3, which=2, cex.axis=1, cex.lab=1, cex.main=1.1,
         main2="Age 3", xlab2="Scale increment Proportions", xlim=c(0,1))
x3<-lines(density(age3), lty=2, lwd=2)
x3<-legend("topright", c("mature", "immature"), col = c(2,3), lty = c(1,1), lwd=2,
           merge = TRUE, bty = "n")
summary(mixmdl3)

age4<-age4[,c(8)]
mixmdl4 <- normalmixEM(age4)
x3<-plot(mixmdl4, which=2, cex.axis=1, cex.lab=1, cex.main=1.1,
         main2="Age 4", xlab2="Scale increment Proportions", xlim=c(0,1))
x3<-lines(density(age4), lty=2, lwd=2)
x3<-legend("topright", c("mature", "immature"), col = c(2,3), lty = c(1,1), lwd=2,
           merge = TRUE, bty = "n")
summary(mixmdl4)

age5<-age5[,c(8)]
mixmdl5 <- normalmixEM(age5)
x3<-plot(mixmdl5, which=2, cex.axis=1, cex.lab=1, cex.main=1.1,
         main2="Age 5", xlab2="Scale increment Proportions", xlim=c(0,1))
x3<-lines(density(age5), lty=2, lwd=2)
x3<-legend("topright", c("mature", "immature"), col = c(2,3), lty = c(1,1), lwd=2,
           merge = TRUE, bty = "n")
summary(mixmdl5)

age6<-age6[,c(8)]
mixmdl6 <- normalmixEM(age6)
x3<-plot(mixmdl6, which=2, cex.axis=1, cex.lab=1, cex.main=1.1,
         main2="Age 6", xlab2="Scale increment Proportions", xlim=c(0,1))
x3<-lines(density(age6), lty=2, lwd=2)
x3<-legend("topright", c("mature", "immature"), col = c(2,3), lty = c(1,1), lwd=2,
           merge = TRUE, bty = "n")
summary(mixmdl6)
dev.off()

merge %>%
  filter(age == 2) -> age2
merge %>%
  filter(age == 3) -> age3
merge %>%
  filter(age == 4) -> age4
merge %>%
  filter(age == 5) -> age5
merge %>%
  filter(age == 6) -> age6 

png(file='figs/Mixture_Models_Transform.png', res=200, width=13, height=10, units ="in") 
par(mfrow=c(3,2)) 
age2<-age2[,c(9)]
mixmdl2 <- normalmixEM(age2)
x3<-plot(mixmdl2, which=2, cex.axis=1, cex.lab=1, cex.main=1.1,
         main2="Age 2", xlab2="Scale increment Proportions (T)", xlim=c(0,1))
x3<-lines(density(age2), lty=2, lwd=2)
x3<-legend("topright", c("mature", "immature"), col = c(2,3), lty = c(1,1), lwd=2,
           merge = TRUE, bty = "n")
summary(mixmdl2)

age3<-age3[,c(9)]
mixmdl3 <- normalmixEM(age3)
x3<-plot(mixmdl3, which=2, cex.axis=1, cex.lab=1, cex.main=1.1,
         main2="Age 3", xlab2="Scale increment Proportions (T)", xlim=c(0,1))
x3<-lines(density(age3), lty=2, lwd=2)
x3<-legend("topright", c("mature", "immature"), col = c(2,3), lty = c(1,1), lwd=2,
           merge = TRUE, bty = "n")
summary(mixmdl3)

age4<-age4[,c(9)]
mixmdl4 <- normalmixEM(age4)
x3<-plot(mixmdl4, which=2, cex.axis=1, cex.lab=1, cex.main=1.1,
         main2="Age 4", xlab2="Scale increment Proportions (T)", xlim=c(0,1))
x3<-lines(density(age4), lty=2, lwd=2)
x3<-legend("topright", c("mature", "immature"), col = c(2,3), lty = c(1,1), lwd=2,
           merge = TRUE, bty = "n")
summary(mixmdl4)

age5<-age5[,c(9)]
mixmdl5 <- normalmixEM(age5)
x3<-plot(mixmdl5, which=2, cex.axis=1, cex.lab=1, cex.main=1.1,
         main2="Age 5", xlab2="Scale increment Proportions (T)", xlim=c(0,1))
x3<-lines(density(age5), lty=2, lwd=2)
x3<-legend("topright", c("mature", "immature"), col = c(2,3), lty = c(1,1), lwd=2,
           merge = TRUE, bty = "n")
summary(mixmdl5)

age6<-age6[,c(9)]
mixmdl6 <- normalmixEM(age6)
x3<-plot(mixmdl6, which=2, cex.axis=1, cex.lab=1, cex.main=1.1,
         main2="Age 6", xlab2="Scale increment Proportions (T)", xlim=c(0,1))
x3<-lines(density(age6), lty=2, lwd=2)
x3<-legend("topright", c("mature", "immature"), col = c(2,3), lty = c(1,1), lwd=2,
           merge = TRUE, bty = "n")
summary(mixmdl6)
dev.off()

merge %>%
  filter(age == 2) -> age2
merge %>%
  filter(age == 3) -> age3
merge %>%
  filter(age == 4) -> age4
merge %>%
  filter(age == 5) -> age5
merge %>%
  filter(age == 6) -> age6

png(file='figs/Mixture_Models_Measurement.png', res=200, width=13, height=10, units ="in") 
par(mfrow=c(3,2)) 
age2<-age2[,c(10)]
mixmdl2 <- normalmixEM(age2)
x3<-plot(mixmdl2, which=2, cex.axis=1, cex.lab=1, cex.main=1.1,
         main2="Age 2", xlab2="Scale increment (mm)", xlim=c(0,4))
x3<-lines(density(age2), lty=2, lwd=2)
x3<-legend("topright", c("mature", "immature"), col = c(2,3), lty = c(1,1), lwd=2,
           merge = TRUE, bty = "n")
summary(mixmdl2)

age3<-age3[,c(10)]
mixmdl3 <- normalmixEM(age3)
x3<-plot(mixmdl3, which=2, cex.axis=1, cex.lab=1, cex.main=1.1,
         main2="Age 3", xlab2="Scale increment (mm)", xlim=c(0,4))
x3<-lines(density(age3), lty=2, lwd=2)
x3<-legend("topright", c("mature", "immature"), col = c(2,3), lty = c(1,1), lwd=2,
           merge = TRUE, bty = "n")
summary(mixmdl3)

age4<-age4[,c(10)]
mixmdl4 <- normalmixEM(age4)
x3<-plot(mixmdl4, which=2, cex.axis=1, cex.lab=1, cex.main=1.1,
         main2="Age 4", xlab2="Scale increment (Proportions (mm)", xlim=c(0,4))
x3<-lines(density(age4), lty=2, lwd=2)
x3<-legend("topright", c("mature", "immature"), col = c(2,3), lty = c(1,1), lwd=2,
           merge = TRUE, bty = "n")
summary(mixmdl4)

age5<-age5[,c(10)]
mixmdl5 <- normalmixEM(age5)
x3<-plot(mixmdl5, which=2, cex.axis=1, cex.lab=1, cex.main=1.1,
         main2="Age 5", xlab2="Scale increment (Proportions (mm)", xlim=c(0,4))
x3<-lines(density(age5), lty=2, lwd=2)
x3<-legend("topright", c("mature", "immature"), col = c(2,3), lty = c(1,1), lwd=2,
           merge = TRUE, bty = "n")
summary(mixmdl5)

age6<-age6[,c(10)]
mixmdl6 <- normalmixEM(age6)
x3<-plot(mixmdl6, which=2, cex.axis=1, cex.lab=1, cex.main=1.1,
         main2="Age 6", xlab2="Scale increment (mm)", xlim=c(0,4))
x3<-lines(density(age6), lty=2, lwd=2)
x3<-legend("topright", c("mature", "immature"), col = c(2,3), lty = c(1,1), lwd=2,
           merge = TRUE, bty = "n")
summary(mixmdl6)
dev.off()

#C) Boxplot Figures 
merge %>% 
  ggplot(aes(x = age, y = outer_prop, color = mature)) +
  geom_jitter(size = 1) + labs(x = "Age",y = "Outer Increment Proportion") + 
  theme(legend.title=element_blank(), legend.position=c(.8,.9)) +
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_fill_manual(values=c("#999999", "#E69F00" ))-> plot1

merge %>% 
 ggplot(aes(x = age, y = outer_prop, color = mature)) + labs(x = "Age",
                                                              y = "Outer Increment Proportion")  +
  geom_boxplot() + theme(legend.position= "none") +
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_fill_manual(values=c("#999999", "#E69F00" ))-> plot2

dataset %>% 
  ggplot(aes(x = mature, y = outer_prop, color=mature)) + labs(x = "",
                                                               y = "Outer Increment Proportion (Age 3)")  +
  geom_boxplot() + theme(legend.position= "none") +
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_fill_manual(values=c("#999999", "#E69F00" ))-> plot3

cowplot::plot_grid(plot1, plot3, align = "vh", nrow = 1, ncol=2)
ggsave("figs/boxplot.png", dpi = 500, height = 6, width = 8, units = "in")


#test for difference in proportions for mature/immature outer ring
res <- t.test(age3mature$aprop, age3immature$aprop) 

#Generalized Linear models (age-3 only)
#merge %>% 
#  do(A1 = glm(maturity ~ outer_prop, data = merge, family = binomial(link=logit)),
#     A2 = glm(maturity ~ outer_prop, data = merge[merge$age==3,], family= binomial(link=logit))) -> lm_out

fit3 <- glm(maturity ~ (outer_prop) , family = binomial, data = dataset) 
AICc(fit3)
Anova(fit3)
RsqGLM(fit3)#peudo R2 
summary(fit3)

dataset %>% 
  do(A2 = glm(maturity ~ outer_prop, data = dataset, family = binomial(link=logit))) -> lm_out
head(augment(lm_out,dataset, type.residuals="pearson"))
outlierTest(fit3) #Bonferroni p-values
residualPlots(fit3) #lack-of fit curvature test
marginalModelPlots(fit3) #marginal model plots
mmp(fit3, dataset$outer_prop, xlab="outer proportion", ylab="maturity" , family="A")
#remove datapoint 92 to determine differnce in coefficients
fit3outlier<-update(fit3, subset=-c(92))
compareCoefs(fit3, fit3outlier)

lm_out %>% 
  tidy(A2) %>% 
  mutate(model = "fit_age3") -> A2
write_csv(A2, "output/lm.csv")

lm_out %>% 
  glance(A2) %>% 
  mutate(model = "fit_age3") -> A2
write_csv(A2, "output/lm_R2.csv")

lm_out %>% # Pearson residuals against covariate
  augment(A2) %>% 
  mutate(resid = (.resid))%>% 
  ggplot(aes(x = outer_prop, y = resid)) +
  geom_hline(yintercept = 0, lty=2) + 
  geom_point(color ="grey50") + 
  scale_x_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4,0.5), limits = c(0, 0.5))+
  scale_y_continuous(breaks = c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2), limits = c(-2,2)) +
  geom_smooth(aes(colour = outer_prop, fill = outer_prop)) +
  labs(y = "Pearson residuals", x =  "Outer proportion")-> plot1

lm_out %>% #Pearson residuals against fitted
  augment(A2) %>% 
  mutate(resid = (.resid),
         fit = (.fitted)) %>% 
  ggplot(aes(x = fit, y = resid)) +
  geom_point(color ="grey50") + 
  scale_x_continuous(breaks = c(-1.0, -0.5, 0, 0.5), limits = c(-1, 0.5))+
  scale_y_continuous(breaks = c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2), limits = c(-2,2)) +
  geom_smooth(aes(colour = fit, fill = fit)) +
  geom_hline(yintercept = 0, lty=2) + 
  labs(y = "Pearson residuals", x =  "Fitted values")-> plot2


outer.data <- data.frame(outer_prop = seq(0, 1, 0.1))# Create a temporary data frame of hypothetical values
predicted.data <- as.data.frame(predict(fit3, newdata = outer.data, # Predict the fitted values given the model and hypothetical data
                                        type="link", se=TRUE))
new.data <- cbind(outer.data, predicted.data)# Combine the hypothetical data and predicted values
std <- qnorm(0.95 / 2 + 0.5)# Calculate confidence intervals
new.data$ymin <- fit3$family$linkinv(new.data$fit - std * new.data$se)
new.data$ymax <- fit3$family$linkinv(new.data$fit + std * new.data$se)
new.data$fit <- fit3$family$linkinv(new.data$fit)  # Rescale to 0-1

dataset %>% #outer verus maturity
ggplot(aes(x=outer_prop, y=maturity)) +
  geom_point() + 
  geom_ribbon(data=new.data, aes(y=fit, ymin=ymin, ymax=ymax), alpha=0.5) + 
  geom_line(data=new.data, aes(y=fit)) + 
  labs(x="Outer proportions", y="Maturity") + 
  scale_x_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4,0.5), limits = c(0, 0.5)) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8,1.0), limits = c(0, 1.0)) -> plot3


lm_out %>% #Cook's distance plot
  augment(A2) %>% 
  mutate(cooksd = (.cooksd),
         count = 1:72,
         name= ifelse(cooksd >0.05, count, ""))%>% 
  ggplot(aes(x = count, y = cooksd, label=name)) +
  geom_bar(stat = "identity", colour = "grey50", 
           fill = "lightgrey",alpha=.7,
           width = 0.8, position = position_dodge(width = 0.2)) + 
  geom_text(size = 3, position = position_stack(vjust = 1.1)) + 
  scale_y_continuous(breaks = c(0, 0.05, 0.1, 0.15, .2), limits = c(0,0.2)) +
  labs(y = "Cook's distance", x =  "Index") -> plot4

lm_out %>% #leverage plot
  augment(A2) %>% 
  mutate(hat= (.hat),
         count = 1:72,
         name= ifelse(hat >0.05, count, "")) %>% 
  ggplot(aes(x = count, y = hat, label=name)) +
  geom_bar(stat = "identity", colour = "grey50", 
           fill = "lightgrey",alpha=.7,
           width = 0.8, position = position_dodge(width = 0.2)) + 
  geom_text(size = 3, position = position_stack(vjust = 1.1)) + 
  scale_y_continuous(breaks = c(0, 0.05, 0.1, 0.15, .2), limits = c(0,0.2)) +
  labs(y = "hat-values", x =  "Index") -> plot5

lm_out %>% #Pearson by index
  augment(A2) %>% 
  mutate(resid = (.resid),
         count = 1:72) %>% 
  ggplot(aes(x = count, y = resid)) +
  geom_bar(stat = "identity", colour = "grey50", 
           fill = "lightgrey",alpha=.7,
           width = 0.8, position = position_dodge(width = 0.2)) + 
  scale_y_continuous(breaks = c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2), limits = c(-2,2)) +
  labs(y = "Pearson residuals", x =  "Index") -> plot6


cowplot::plot_grid(plot1, plot2, plot3, plot4, plot5, plot6,  align = "vh", nrow = 3, ncol=2)
ggsave("figs/glm_diagnostics.png", dpi = 500, height = 6, width = 8, units = "in")

par(mfrow=c(1,1))
png(filename="figs/glm_diagnostics1.png", width = 6, height = 6, units = 'in', res = 400)
windowsFonts(A = windowsFont("Times New Roman"))
influencePlot(fit3, xlab="Hat-values", ylab="Studentized residuals", family="A")
dev.off()
#png(filename="figs/glm_diagnostics2.png", width = 8, height = 6, units = 'in', res = 400)
#influenceIndexPlot(fitall, family="A", main="")
#dev.off()

#quantile quantile plots
par(mfrow = c(1,1)) #qqplot
plot(fit3, 2, family="Times")
tau_est <- fit3$summary.hyperpar[1,4]
nu_est <- fit3$summary.hyperpar[2,4]

#Generalized Linear models (ages 2-7)
fit0<- glm(maturity ~ 1 , family = binomial, data=merge) 
fit1<- glm(maturity ~  age, data=merge, family = binomial) 
fit2 <- glm(maturity ~ (outer_prop) , family = binomial, data=merge)
fit3 <- glm(maturity ~ (outer_prop) + age, family = binomial, data=merge) 
AICc(fit3)
Anova(fit3)
RsqGLM(fit3)#peudo R2 
summary(fit3)

#Generalized Additive Models
#fit3 <- glm(maturity ~ (outer_prop) , family = binomial, data = merge[merge$age==3,]) 
#fit5 <- gam(maturity ~ s(outer_prop, by = age) , family = binomial, data=dataset) 
#fit6 <- gam(maturity ~ s(outer_prop, by = age) + age, family = binomial, data=dataset) 
AIC(fit1, fit2, fit3, fit4)

#diagnostics of best model
coef <- coefficients(fit3) # coefficients
resid <- residuals(fit3, type="pearson")
resid<-as.data.frame(resid)# residuals
pred <- predict(fit3) # fitted values
rsq <- summary(fit3)$r.squared # R-sq for the fit
se <- summary(fit)$sigma # se of the fit
plot.gam(fit3)
gam.check(fit3)
shapiro.test(resid)
durbinWatsonTest(resid)
dwtest(fit3)
data = dataset[dataset$age==3,]
data <-cbind(data, resid)
xnumeracy <- seq (0.1, 0.4, 0.01)
ynumeracy <- predict(fit3, list(outer_prop=xnumeracy),type="response")
#dataset<-cbind(xnumeracy, ynumeracy)
plot(data$outer_prop, data$maturity, pch = 16, xlab = "NUMERACY SCORE", ylab = "ADMISSION")
lines(xnumeracy, ynumeracy, col = "red", lwd = 2)

#Durbin-Watson Test
#2 is no autocorrelation.
#0 to <2 is positive autocorrelation (common in time series data).
#>2 to 4 is negative autocorrelation (less common in time series data).
#A rule of thumb is that test statistic values in the range of 1.5 to 2.5 are relatively normal
#Field(2009) suggests that values under 1 or more than 3 are a definite cause for concern.
#Field, A.P. (2009). Discovering statistics using SPSS: and sex and drugs and rock ‘n’ roll (3rd edition). London:Sage.

#If the model distributional assumptions are met then usually these plots
#should be close to a straight line (although discrete data can yield marked random departures from this line). 
par(mfrow=c(2,2))
qqnorm(residuals(fit3),pch=19,cex=.3)
qq.gam(pif,pch=19,cex=.3)
qq.gam(pif,rep=100,level=.9)
qq.gam(pif,rep=100,level=1,type="pearson",pch=19,cex=.2)
dev.off()
#cook's distance plot [Zuur et al. (2013): A beginner's Guide to GLM and GLMM with R]
#plot(cooks.distance(fit3), ylim=c(0,1), ylab="Cook distance values", type="h")
#qqline plot
#Pearson residuals vs. continous covariate

#Residual plots
fit3.diag <- glm.diag(fit3)
glm.diag.plots(fit3, fit3.diag)

# Residual vs. fitted
E2 <- resid(fit3, type="pearson")
F2 <- fitted(fit3, type="response")
plot(x=F2, y=E2, xlab="fitted values", ylab="Pearson residuals")
abline(h=0, lty=2)

# Pearson residuals vs. continous covariate
plot(x=df$SYNTAXSCORE, y=E2, xlab="SYNTAXSCORE", ylab="Pearson residuals")
abline(h=0, lty=2)
#Zuur et al. (2013): A beginner's Guide to GLM and GLMM with R.
