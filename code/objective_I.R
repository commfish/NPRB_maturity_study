#notes---
#author: Sara E Miller 
#contact: sara.miller@alaska.gov; 907-465-4245
#Last edited: May, 2019

# load libraries----
# devtools::install_github("ben-williams/FNGr")
library(extrafont)
loadfonts(device="win")
library(tidyverse)
library(FNGr)
library(broom)
library(cowplot)
library(mixtools)
library(mixdist)
library(fitdistrplus)
library(gam)
library(car)
library(lmtest)
library(mgcv)
#library(visreg)

theme_set(theme_sleek())

asintrans <- function(p) {asin(sqrt(p))} # arctransformation function

eda.norm <- function(x, ...) #normality test function
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
  filter(!(scale_region%in% c("OOA (OUT OF AREA)")))-> clean_dataset #n=568 useable samples

clean_dataset %>% 
  dplyr::select(age,maturity_state_histology, scale_region, sex_histology) %>% 
  group_by(maturity_state_histology, age, sex_histology) %>% 
  summarize(n=n())-> table1 #sample sizes for study

clean_dataset %>% 
  dplyr::select(maturity_state_histology, maturation_status_histology) %>% 
  group_by(maturity_state_histology,maturation_status_histology) %>% 
  summarize(n=n())-> table2 #sample sizes for study

clean_dataset %>% 
  dplyr::select(scale_region) %>% 
  group_by(scale_region) %>% 
  summarize(n=n())-> table3 #sample sizes for study

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
  dplyr::select(image_name,year, age, sex_histology, maturation_status_histology, mature, max, outer_prop, aprop, anu_adj) -> merge
         
#Histograms of outer ring 
#datasets by ages
merge %>%
  filter(age == 2) %>%
  filter(mature == "mature")-> age2mature
merge %>%
  filter(age == 2) %>%
  filter(mature == "immature")-> age2immature
merge %>%
  filter(age == 3) %>%
  filter(mature == "mature")-> age3mature
merge %>%
  filter(age == 3) %>%
  filter(mature == "immature")-> age3immature
merge %>%
  filter(age == 4) %>%
  filter(mature == "mature")-> age4mature
merge %>%
  filter(age == 4) %>%
  filter(mature == "immature")-> age4immature
merge %>%
  filter(age == 5) %>%
  filter(mature == "mature")-> age5mature
merge %>%
  filter(age == 5) %>%
  filter(mature == "immature")-> age5immature
merge %>%
  filter(age == 6) %>%
  filter(mature == "mature")-> age6mature
merge %>%
  filter(age == 6) %>%
  filter(mature == "immature")-> age6immature

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
  filter(age == 2) -> age2
merge %>%
  filter(age == 3) -> age3
merge %>%
  filter(age == 4) -> age4
merge %>%
  filter(age == 5) -> age5
merge %>%
  filter(age == 6) -> age6 
 
ggplot(age2, aes(x=outer_prop, color=mature, fill=mature)) + geom_histogram(alpha=0.5, position = 'identity') +
  ylab("Frequency")+ xlab("Outer increment proportion") +
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6), limits = c(0,0.6))+
  ggtitle("Age 2; n=215") + theme(legend.position="none") +
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_fill_manual(values=c("#999999", "#E69F00" )) -> plot1

ggplot(age2, aes(x=outer_prop, color=mature, fill=mature)) +
  geom_density(alpha=0.5, adjust=1) +scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_fill_manual(values=c("#999999", "#E69F00" )) + theme(legend.title=element_blank(), legend.position=c(.85,.85)) + 
  ylab("Density")+ xlab("Outer increment proportion") +
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6), limits = c(0,0.6))-> plot2

ggplot(age3, aes(x=outer_prop, color=mature, fill=mature)) + geom_histogram(alpha=0.5, position = 'identity') +
  ylab("Frequency")+ xlab("Outer increment proportion") +
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6), limits = c(0,0.6))+
  ggtitle("Age 3; n=243") + theme(legend.position="none") +
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_fill_manual(values=c("#999999", "#E69F00" )) -> plot3

ggplot(age3, aes(x=outer_prop, color=mature, fill=mature)) +
  geom_density(alpha=0.5, adjust=1) +scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_fill_manual(values=c("#999999", "#E69F00" )) + theme(legend.title=element_blank(), legend.position=c(.85,.85)) +
  ylab("Density")+ xlab("Outer increment proportion") +
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6), limits = c(0,0.6))-> plot4

ggplot(age4, aes(x=outer_prop, color=mature, fill=mature)) + geom_histogram(alpha=0.5, position = 'identity') +
  ylab("Frequency")+ xlab("Outer increment proportion") +
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6), limits = c(0,0.6)) +
  ggtitle("Age 4; n=93") + theme(legend.position="none") +
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_fill_manual(values=c("#999999", "#E69F00" )) -> plot5

ggplot(age4, aes(x=outer_prop, color=mature, fill=mature)) +
  geom_density(alpha=0.5, adjust=1) +scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_fill_manual(values=c("#999999", "#E69F00" )) + theme(legend.title=element_blank(), legend.position=c(.85,.85)) +
  ylab("Density")+ xlab("Outer increment proportion") +
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6), limits = c(0,0.6))-> plot6

ggplot(age5, aes(x=outer_prop, color=mature, fill=mature)) + geom_histogram(alpha=0.5, position = 'identity') +
  ylab("Frequency")+ xlab("Outer increment proportion") +
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6), limits = c(0,0.6)) +
  ggtitle("Age 5; n=47") + theme(legend.position="none") +
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_fill_manual(values=c("#999999", "#E69F00" )) -> plot7

ggplot(age5, aes(x=outer_prop, color=mature, fill=mature)) +
  geom_density(alpha=0.5, adjust=1) +scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_fill_manual(values=c("#999999", "#E69F00" )) + theme(legend.title=element_blank(), legend.position=c(.85,.85)) +
  ylab("Density")+ xlab("Outer increment proportion") +
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6), limits = c(0,0.6)) -> plot8

ggplot(age6, aes(x=outer_prop, color=mature, fill=mature)) + geom_histogram(alpha=0.5, position = 'identity') +
  ylab("Frequency")+ xlab("Outer increment proportion") +
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6), limits = c(0,0.6)) +
  ggtitle("Age 6; n=68") + theme(legend.position="none") +
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_fill_manual(values=c("#999999", "#E69F00" )) -> plot9

ggplot(age6, aes(x=outer_prop, color=mature, fill=mature)) +
  geom_density(alpha=0.5, adjust=1) +scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_fill_manual(values=c("#999999", "#E69F00" )) + theme(legend.title=element_blank(), legend.position=c(.85,.85)) +
  ylab("Density")+ xlab("Outer increment proportion") +
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6), limits = c(0,0.6))-> plot10
cowplot::plot_grid(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, plot9, plot10, align = "vh", nrow = 5, ncol=2)
ggsave("figs/histogram.png", dpi = 500, height = 10, width =8, units = "in")

#Fit gaussian mixture models
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


#test for difference in proportions for mature/immature outer ring
res <- t.test(age3mature$aprop, age3immature$aprop) 

#GAM 
merge %>% 
  mutate(age=as.factor(age),
         maturity = ifelse(mature == "mature", 1, 0)) -> dataset
fit <- glm(maturity ~ (outer_prop) , family = binomial, data=dataset)
fit3 <- glm(maturity ~ (outer_prop) , family = binomial, data = dataset[dataset$age==3,]) 
fit1 <- glm(maturity ~ (outer_prop) + age, family = binomial, data=dataset) 
fit2 <- gam(maturity ~ s(outer_prop, by = age) , family = binomial, data=dataset) 
fit3 <- gam(maturity ~ s(outer_prop, by = age) + age, family = binomial, data=dataset) 
fit4 <- gam(maturity ~  age, family = binomial, data=dataset) 
summary(fit)
summary(fit1)
summary(fit2)
summary(fit3)
summary(fit4)
plot(fit1)
plot(fit2)
plot(fit3)
AIC(fit, fit1, fit2, fit3, fit4)

fit5 <- gam(outer_prop ~ age, data = dataset)
summary(fit5)

#diagnostics of best model
coef <- coefficients(fit3) # coefficients
resid <- residuals(fit3) # residuals
pred <- predict(fit3) # fitted values
rsq <- summary(fit3)$r.squared # R-sq for the fit
se <- summary(fit)$sigma # se of the fit
plot.gam(fit3)
gam.check(fit3)
shapiro.test(resid)
durbinWatsonTest(resid)
dwtest(fit3)
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
#Generalized Linear models 
dataset %>% 
  do(A1 = glm(maturity ~ outer_prop, data = dataset, family = binomial),
     A2 = glm(maturity ~ outer_prop, data = dataset[dataset$age==3,], family= binomial)) -> lm_out

lm_out %>% 
  tidy(A1) %>% 
  mutate(model = "fit") -> A1
lm_out %>% 
  tidy(A2) %>% 
  mutate(model = "fit_age3") -> A2
x<- rbind(A1, A2) #combine data for all zones
write_csv(x, "output/lm.csv")

lm_out %>% 
  glance(A1) %>% 
  mutate(model = "fit") -> A1
lm_out %>% 
  glance(A2) %>% 
  mutate(model = "fit_age3") -> A2
x<- rbind(A1, A2) #combine data for all zones
write_csv(x, "output/lm_R2.csv")

lm_out %>% 
  augment(A1) %>% 
  mutate(fit = (.fitted)) %>% 
  ggplot(aes(x = outer_prop, y = maturity)) +
  geom_point(color ="grey50")+geom_line(aes(maturity, fit), color = "black") + 
  annotate("text", x = 200, y=7, label="A1", family="Times New Roman")+
  scale_x_continuous(breaks = c(0, .25,0.5,.75, 1), limits = c(0,1))+
  labs(y = "Maturity", x =  "Outer Proportion")-> A1

lm_out %>% 
  augment(A2) %>% 
  mutate(fit = (.fitted)) %>% 
  ggplot(aes(x = outer_prop, y = maturity)) +
  geom_point(color ="grey50")+geom_line(aes(maturity, fit), color = "black") + 
  #annotate("text", x = 0, y=2, label="Age 3", family="Times New Roman")+
  scale_x_continuous(breaks = c(0, .25,0.5,.75, 1), limits = c(0,1))+
  labs(y = "Maturity", x =  "Outer Proportion")-> plot1

lm_out %>% 
  augment(A2) %>% 
  mutate(resid = (.resid),
         fit = (.fitted)) %>% 
  ggplot(aes(x = fit, y = resid)) +
  geom_point(color ="grey50") + 
  #annotate("text", x = 0, y=2, label="Age 3", family="Times New Roman")+
  scale_x_continuous(breaks = c(0, .25,0.5,.75, 1), limits = c(0,1)) +
  geom_hline(yintercept = 0, lty=2) + geom_smooth(method="auto") +
  labs(y = "Residuals", x =  "Fitted Values")-> plot2

#Residual plots
op <- par(oma=c(4,1,1,1))
fit3.diag <- glm.diag(fit3)
glm.diag.plots(fit3, fit3.diag)
par(op)

ggplot(data=dataset, aes(x = fitted(fit3), y = residuals(fit3))) + geom_point() +
  labs(y = "Residuals", x =  "Fitted Values") +
abline(h=0, lty=2) +
lines(smooth.spline(fitted(fit3), residuals(fit3)))#plot of the fitted vs. residuals (upper right)
 
'#Boxplot Figures 
dataset %>% 
  ggplot(aes(x = age, y = outer_prop, color = mature)) +
  geom_jitter(size = 1) + labs(x = "Age",y = "Outer Increment Proportion") + 
  theme(legend.title=element_blank(), legend.position=c(.85,.9)) +
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_fill_manual(values=c("#999999", "#E69F00" ))-> plot1

dataset %>% 
  ggplot(aes(x = age, y = outer_prop, color = mature)) + labs(x = "Age",
  y = "Outer Increment Proportion")  +
  geom_boxplot() + theme(legend.position= "none") +
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_fill_manual(values=c("#999999", "#E69F00" ))-> plot2

cowplot::plot_grid(plot1, plot2, align = "vh", nrow = 1, ncol=2)
ggsave("figs/boxplot.png", dpi = 500, height = 4, width = 6, units = "in")
