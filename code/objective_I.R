#author: Sara E Miller 
#contact: sara.miller@alaska.gov; 907-465-4245
#Last edited: February, 2020
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4885900/
#https://stats.stackexchange.com/questions/82105/mcfaddens-pseudo-r2-interpretation
#S:\Region1Shared-DCF\Research\Herring-Dive Fisheries\Herring\NPRB-Herring Maturity Study\datasheets\datasheets\MTA Lab Datasheets\fall field study-MTA lab measurements\remeasurements_corrected data\ MTA Lab data_length_weight_Jan_2019
#matched to histology data sheet (S:\Region1Shared-DCF\Research\Herring-Dive Fisheries\Herring\NPRB-Herring Maturity Study\datasheets\datasheets\Histology-Pathologist report\Alaska Department of Fish and Game 10_17_2018_corrected)
#the following code then takes the histology dataset data_11_26_2018 and matches it to the MTA lab increment data by image name

# load libraries----
# devtools::install_github("ben-williams/FNGr")
source("code/functions.R")
library(extrafont)
loadfonts(device="win")
library(nlme)
library(tidyverse)
library(fngr)
library(broom)
library(cowplot)
library(mixtools)
library(car)
library(modEvA) #pseudo-R squared
library(MASS)
library(digest)
library(ResourceSelection)
library("ggpubr")
#library(lmtest)
#library(mgcv)
#library(visreg)
#library(boot)
#library(AICcmodavg) #AICc
#library(rms)
#library(mixdist)
#library(fitdistrplus)
#library(gam)
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
  filter(!(scale_region%in% c("OOA (OUT OF AREA)", "A", "C", "D", "E", "G", "H"))) %>% 
  mutate(maturity = ifelse(maturation_status_histology == 1, 'immature', 'mature')) -> clean_dataset #n=296 useable samples

#sample sizes for study 
clean_dataset %>% 
  dplyr::select(age,maturity, scale_region, sex_histology) %>% 
  group_by(maturity, age, sex_histology) %>% 
  summarize(n=n())-> table1 

clean_dataset %>% 
  dplyr::select(maturity, maturation_status_histology) %>% 
  group_by(maturity,maturation_status_histology) %>% 
  summarize(n=n())-> table2 

clean_dataset %>% 
  dplyr::select(scale_region) %>% 
  group_by(scale_region) %>% 
  summarize(n=n())-> table3 

#calculate proportion of outer scale ring
clean_dataset %>% 
  mutate(anu1 = ifelse(is.na(Increment1), 0 , Increment1),
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
         rad7 = sum(anu1, anu2, anu3, anu4, anu5, anu6, anu7))-> sample1
 
# anu_adj is measurement of outer ring
# radcap is cumulative rowth
sample1 %>%
  mutate(anu_adj = ifelse(anu1>0 & anu2==0 & anu3 ==0 & anu4 ==0 & anu5==0 & anu6==0 & anu7==0, anu1,
                          ifelse(anu1>0 & anu2>0& anu3 ==0 & anu4 ==0 & anu5==0 & anu6==0 & anu7==0, anu2,
                                 ifelse(anu1>0 & anu2>0 & anu3>0 & anu4 ==0 & anu5==0 & anu6==0 & anu7==0, anu3,
                                        ifelse(anu1>0 & anu2>0 & anu3>0 & anu4>0 & anu5==0 & anu6==0 & anu7==0, anu4,
                                               ifelse(anu1>0 & anu2>0 & anu3>0 & anu4>0 & anu5>0 & anu6==0 & anu7==0, anu5,
                                                      ifelse(anu1>0 & anu2>0 & anu3>0 & anu4>0 & anu5>0 & anu6>0 & anu7==0, anu6,anu7))))))) %>%
  filter(sex_histology == "Female") %>% 
  dplyr::select(image_name,year, age, sex_histology, maturation_status_histology, maturity, anu_adj, radcap, length_mm) -> merge
write.csv(merge, "data/cpue_new_test.csv") 

# Exploratory Plots----
# Histograms of outer ring---- 
merge %>%
  filter(age == 2 & maturity == "mature")-> age2mature
merge %>%
  filter(age == 2 & maturity == "immature")-> age2immature
merge %>%
  filter(age == 3 & maturity == "mature")-> age3mature
merge %>%
  filter(age == 3 & maturity == "immature")-> age3immature
merge %>%
  filter(age == 4 & maturity == "mature")-> age4mature
merge %>%
  filter(age == 4 & maturity == "immature")-> age4immature
merge %>%
  filter(age == 5 & maturity == "mature")-> age5mature
merge %>%
  filter(age == 5 & maturity == "immature")-> age5immature
merge %>%
  filter(age == 6 & maturity == "mature")-> age6mature
merge %>%
  filter(age == 6 & maturity == "immature")-> age6immature

# test each age/maturity for normality
eda.norm(as.numeric(age2immature$aprop))
eda.norm(as.numeric(age2mature$aprop))
eda.norm(as.numeric(age3immature$aprop))
eda.norm(as.numeric(age3mature$aprop))
eda.norm(as.numeric(age4immature$aprop))
eda.norm(as.numeric(age4mature$aprop))
eda.norm(as.numeric(age5mature$aprop)) #normal
eda.norm(as.numeric(age6mature$aprop)) #normal

# Frequencies by age----
merge %>%
  filter(age == 2) %>%
ggplot(., aes(x=anu_adj, color=maturity, fill=maturity)) + geom_histogram(alpha=0.5, position = 'identity') +
  ylab("Frequency")+ xlab("Outer increment measurement (mm)") +
  scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2), limits = c(0,2))+
  ggtitle("Age 2; n = 64") + theme(legend.position="none") +
  scale_color_manual(values=c("#999999", "black")) +theme(legend.title=element_blank(), legend.position=c(.15,.85)) +
  scale_fill_manual(values=c("#999999", "black" )) -> plot1

merge %>%
  filter(age == 3) %>%
  ggplot(., aes(x=anu_adj, color=maturity, fill=maturity)) + geom_histogram(alpha=0.5, position = 'identity') +
  ylab("Frequency")+ xlab("Outer increment measurement (mm)") +
  scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2), limits = c(0,2))+
  ggtitle("Age 3; n = 51") + theme(legend.position="none") +
  scale_color_manual(values=c("#999999", "black")) +
  scale_fill_manual(values=c("#999999", "black" )) -> plot2

merge %>%
  filter(age == 4) %>%
  ggplot(., aes(x=anu_adj, color=maturity, fill=maturity)) + geom_histogram(alpha=0.5, position = 'identity') +
  ylab("Frequency")+ xlab("Outer increment measurement (mm)") +
  scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2), limits = c(0,2))+
  ggtitle("Age 4; n = 40") + theme(legend.position="none") +
  scale_color_manual(values=c("#999999", "black")) + 
  scale_fill_manual(values=c("#999999", "black" )) -> plot3

merge %>%
  filter(age == 5) %>%
  ggplot(., aes(x=anu_adj, color=maturity, fill=maturity)) + geom_histogram(alpha=0.5, position = 'identity') +
  ylab("Frequency")+ xlab("Outer increment measurement (mm)") +
  scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2), limits = c(0,2))+
  ggtitle("Age 5; n = 23") + theme(legend.position="none") +
  scale_color_manual(values=c("black")) + 
  scale_fill_manual(values=c("black" )) -> plot4

merge %>%
  filter(age == 6) %>%
  ggplot(., aes(x=anu_adj, color=maturity, fill=maturity)) + geom_histogram(alpha=0.5, position = 'identity') +
  ylab("Frequency")+ xlab("Outer increment measurement (mm)") +
  scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2), limits = c(0,2))+
  ggtitle("Age 6; n = 33") + theme(legend.position="none") +
  scale_color_manual(values=c("#999999", "black")) + 
  scale_fill_manual(values=c("#999999", "black" )) -> plot5

cowplot::plot_grid(plot1, plot2, plot3, plot4, plot5, align = "vh", nrow = 2, ncol=3)
ggsave("figs/frequency_obj1.png", dpi = 500, height = 8, width =10, units = "in")

# Histogreams by age----
merge %>%
  filter(age == 2) %>%
  ggplot(., aes(x=anu_adj, color=maturity, fill=maturity)) +ggtitle("Age 2; n = 64")+
  geom_density(alpha=0.5, adjust=1) +scale_color_manual(values=c("#999999", "black")) +
  scale_fill_manual(values=c("#999999", "black" )) + theme(legend.title=element_blank(), legend.position=c(.19,.88)) + 
  ylab("Density")+ xlab("Outer increment measurement (mm)") +
  scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2), limits = c(0,2))-> plot1

merge %>%
  filter(age == 3) %>%
  ggplot(., aes(x=anu_adj, color=maturity, fill=maturity)) +ggtitle("Age 3; n = 51") +
  geom_density(alpha=0.5, adjust=1) +scale_color_manual(values=c("#999999", "black")) +
  scale_fill_manual(values=c("#999999", "black" )) + theme(legend.position="none") +
  ylab("Density")+ xlab("Outer increment measurement (mm)") +
  scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2), limits = c(0,2))-> plot2

merge %>%
  filter(age == 4) %>%
  ggplot(., aes(x=anu_adj, color=maturity, fill=maturity)) +  ggtitle("Age 4; n = 40") + 
  geom_density(alpha=0.5, adjust=1) +scale_color_manual(values=c("#999999", "black")) +
  scale_fill_manual(values=c("#999999", "black" )) + theme(legend.position="none") + 
  ylab("Density")+ xlab("Outer increment measurement (mm)") +
  scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2), limits = c(0,2))-> plot3

merge %>%
  filter(age == 5) %>%
  ggplot(., aes(x=anu_adj, color=maturity, fill=maturity)) +ggtitle("Age 5; n = 23") +
  geom_density(alpha=0.5, adjust=1) +scale_color_manual(values=c("black")) +
  scale_fill_manual(values=c("black" )) + theme(legend.position="none") +
  ylab("Density")+ xlab("Outer increment measurement (mm)") +
  scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2), limits = c(0,2))-> plot4

merge %>%
  filter(age == 6) %>%
  ggplot(., aes(x=anu_adj, color=maturity, fill=maturity)) +
  geom_density(alpha=0.5, adjust=1) +scale_color_manual(values=c("#999999", "black")) +
  scale_fill_manual(values=c("#999999", "black" )) +  ggtitle("Age 6; n = 33")+
  ylab("Density")+ xlab("Outer increment measurement (mm)") +theme(legend.position="none")+
  scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2), limits = c(0,2))-> plot5

cowplot::plot_grid(plot1, plot2, plot3, plot4, plot5, align = "vh", nrow = 2, ncol=3)
ggsave("figs/hist_obj1.png", dpi = 500, height = 8, width =10, units = "in")

# Scatterplot Figures---- 
merge %>% 
  ggplot(data=., aes(x = age, y = anu_adj, color = maturity, shape = maturity)) +
  geom_jitter(size = 1) +
  labs(x = "Age",y = "Outer increment measurement (mm)") + 
  theme(legend.title=element_blank(), legend.position=c(.8,.9), legend.text=element_text(size=12)) +
  scale_color_manual(values=c("#999999", "black")) +
  scale_fill_manual(values=c("#999999", "black" )) +
  scale_shape_manual(values=c(8, 16)) +
  annotate("text",x = 2, y=2, label="A)", family="Arial" ,size=4) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0), limits = c(0, 2.0)) +
  scale_x_continuous(breaks = c(2,3,4,5,6), limits = c(2,6))-> plot1

merge %>%
  mutate(age=as.factor(age))%>% 
  ggplot(data=.,aes(x = age, y = anu_adj)) + labs(x = "Age",y = "Outer increment measurement (mm)") + 
  geom_boxplot(aes(fill = maturity)) +
  scale_color_manual(values=c("#999999", "black")) + theme(legend.title=element_blank(), legend.position=c(.8,.9), legend.text=element_text(size=12)) +
  scale_fill_manual(values=c("#999999", "black" )) + 
  scale_y_continuous(breaks = c(0, 0.5, 1, 1.5, 2), limits = c(0, 2)) +
  annotate("text",x = 0.75, y=2, label="B)", family="Arial" ,size=4)-> plot2

merge %>% 
  ggplot(data=., aes(x = age, y = radcap, color = maturity, shape = maturity)) +
  geom_jitter(size = 1) +
  labs(x = "Age",y = "Cumulative growth (mm)") + 
  theme(legend.title=element_blank(), legend.position=c(.8,.9), legend.text=element_text(size=12)) +
  scale_color_manual(values=c("#999999", "black")) +
  scale_fill_manual(values=c("#999999", "black" )) +
  scale_shape_manual(values=c(8, 16)) +
  annotate("text",x = 2, y=10, label="B)", family="Arial", colour="black", size=4) +
  scale_y_continuous(breaks = c(0, 1,2,3,4,5,6,7,8,9,10), limits = c(0, 10)) +
  scale_x_continuous(breaks = c(2,3,4,5,6), limits = c(2,6))-> plot3

cowplot::plot_grid(plot1, plot2,  align = "vh", nrow = 1, ncol=2)
ggsave("figs/scatterplot_obj1.png", dpi = 500, height = 4, width = 8, units = "in")

# Correlation Analysis---
# age versus length at capture
cor.test(merge$age, merge$length_mm, method=c("pearson", "kendall", "spearman"))

ggscatter(merge, x = "age", y = "length_mm",  add="loess",
          conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Age", ylab = "Length at capture (mm)") -> plot1
#ggscatter(merge, x = "radcap", y = "length_mm", 
#          add = "reg.line", conf.int = TRUE, 
#          cor.coef = TRUE, cor.method = "spearman",
#          xlab = "Total Scale Length", ylab = "Length at capture (mm)") -> plot2
ggscatter(merge, x = "age", y = "radcap", 
          conf.int = TRUE, add="loess", 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Age", ylab = "Total Scale Length")+ annotate("text",x = 6, y=7, label="B)", family="Arial" ,size=4) -> plot3
cowplot::plot_grid(plot1,   align = "vh", nrow = 1, ncol=1)
ggsave("figs/correlation.png", dpi = 500, height = 4, width = 8, units = "in")

par(mfrow=c(1,1))
with(merge,interaction.plot(age, maturity, anu_adj, type="b", pch=c(1,16), ylab ="outer incremental length (mm)",
                            xlab="Age"))
# Generalized Linear models----
## Outer Increment Measurement
merge %>%
  mutate(age = as.factor(age)) %>%
  mutate(maturity = ifelse(maturity == "mature", 1 , 0))-> merge
fit <- glm(maturity ~ (anu_adj) , family = binomial, data = merge) 
fit1 <- glm(maturity ~ (age*anu_adj) , family = binomial, data = merge) 
fit2 <- glm(maturity ~ (age) , family = binomial, data = merge) 
vif(fit)
Anova(fit)
RsqGLM(fit)#peudo R2 
summary(fit)
hoslem.test(merge$maturity, fitted(fit)) #goodness of fit test (ResourceSelection package); https://www.theanalysisfactor.com/r-glm-model-fit/

merge %>% 
  do(A2 = glm(maturity ~ anu_adj *age, data = ., family = binomial(link=logit))) -> lm_out
head(augment(lm_out, merge, type.residuals="pearson"))
outlierTest(fit) #Bonferroni p-values
residualPlot(fit, variable = "fitted", type = "pearson",
             plot = TRUE, quadratic = FALSE, smooth=TRUE) #curvature test
marginalModelPlots(fit) #marginal model plots
mmp(fit, merge$anu_adj, xlab="outer proportion", ylab="maturity" , family="A")


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
  ggplot(aes(x = anu_adj, y = resid)) +
  geom_hline(yintercept = 0, lty=2) + 
  geom_point(color ="grey50") + 
  scale_x_continuous(breaks = c(0, 0.5,1,1.5,2), limits = c(0, 2))+
  scale_y_continuous(breaks = c(-2, -1.5, -1,-0.5, 0, 0.5, 1, 1.5, 2), limits = c(-2, 2)) +
  geom_smooth(aes(colour = anu_adj, fill = anu_adj), colour="black") +
  labs(y = "Pearson residuals", x =  "Outer measurement") +
  geom_text(aes(x = 0, y = 1.9, label="a)"),family="Times New Roman", colour="black", size=5)-> plot1

lm_out %>% #Pearson residuals against fitted
  augment(A2) %>% 
  mutate(resid = (.resid),
         fit = (.fitted)) %>% 
  ggplot(aes(x = fit, y = resid)) +
  geom_point(color ="grey50") + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2))+
  scale_y_continuous(breaks = c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2), limits = c(-2,2)) +
  geom_smooth(aes(colour = fit, fill = fit),colour="black") +
  geom_hline(yintercept = 0, lty=2) + 
  labs(y = "Pearson residuals", x =  "Fitted values") +
  geom_text(aes(x = 0, y = 2, label="b)", hjust = 1),family="Times New Roman", colour="black", size=5)-> plot2

#pf(0.49, 3, 17) # to figure out influence of cook's distance; relate D to the F(p, n-p) distribution and accertain the corresponding percentile value; see Neter et al. book
qf(0.5, 3+1, 211-3-1) # see https://rpubs.com/mpfoley73/460943 and Das and Gogoi 2015 pg 80 
qf(0.5, 4, 207) # p+1; n-p-1


lm_out %>% #Cook's distance plot
  augment(A2) %>% 
  mutate(cooksd = (.cooksd),
         count = 1:211,
         name= ifelse(cooksd >0.84, count, ""))%>% 
  ggplot(aes(x = count, y = cooksd, label=name)) +
  geom_bar(stat = "identity", colour = "grey50", 
           fill = "lightgrey",alpha=.7,
           width = 0.8, position = position_dodge(width = 0.2)) + 
  geom_text(size = 3, position = position_stack(vjust = 1.1)) + 
  scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5), limits = c(0,0.5)) +
  labs(y = "Cook's distance", x =  "Index") +
  geom_text(aes(x = 0, y = 0.5, label="c)"),family="Times New Roman", colour="black", size=5)-> plot3

# greater than 2 or three times p/n, it may be a concern (where p is the number of parameters (i.e., 2) and n is the number of observations (i.e., 51); p/n = 0.04 for this study; Dobson 2002). 
# 3/211 = 0.01; 0.01 *3 =0.03
lm_out %>% #leverage plot
  augment(A2) %>% 
  mutate(hat= (.hat),
         count = 1:211,
         name= ifelse(hat >0.03, count, "")) -> x%>% 
  ggplot(aes(x = count, y = hat, label=name)) +
  geom_bar(stat = "identity", colour = "grey50", 
           fill = "lightgrey",alpha=.7,
           width = 0.8, position = position_dodge(width = 0.2)) + 
  geom_text(size = 3, position = position_stack(vjust = 1.2)) + 
  scale_y_continuous(breaks = c(0, 0.05, 0.1, 0.15), limits = c(0,0.15)) +
  labs(y = "hat-values", x =  "Index") +
  geom_text(aes(x = 1, y = 0.15, label="d)"),family="Times New Roman", colour="black", size=5)-> plot4

lm_out %>% #Pearson by index
  augment(A2) %>% 
  mutate(resid = (.resid),
         count = 1:211) %>% 
  ggplot(aes(x = count, y = resid)) +
  geom_bar(stat = "identity", colour = "grey50", 
           fill = "lightgrey",alpha=.7,
           width = 0.8, position = position_dodge(width = 0.2)) + 
  scale_y_continuous(breaks = c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2), limits = c(-2,2)) +
  labs(y = "Pearson residuals", x =  "Index") +
  geom_text(aes(x = 0, y = 1.9, label="e)"),family="Times New Roman", colour="black", size=5)-> plot5

cowplot::plot_grid(plot1, plot2, plot3, plot4, plot5,  align = "vh", nrow = 3, ncol=2)
ggsave("figs/glm_diagnostics_outer_increment.png", dpi = 500, height = 6, width = 8, units = "in")


# Generalized Linear models----
## Cumulative Growth
fit <- glm(maturity ~ (anu_adj +radcap) , family = binomial, data = merge) 
vif(fit)
Anova(fit)
RsqGLM(fit)#peudo R2 
summary(fit)
hoslem.test(merge$maturity, fitted(fit)) #goodness of fit test (ResourceSelection package); https://www.theanalysisfactor.com/r-glm-model-fit/

merge %>% 
  do(A2 = glm(maturity ~ radcap +age, data = ., family = binomial(link=logit))) -> lm_out
head(augment(lm_out, merge, type.residuals="pearson"))
outlierTest(fit) #Bonferroni p-values
residualPlot(fit, variable = "fitted", type = "pearson",
             plot = TRUE, quadratic = FALSE, smooth=TRUE) #curvature test
marginalModelPlots(fit) #marginal model plots
mmp(fit, merge$radcap, xlab="outer proportion", ylab="maturity" , family="A")


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
  ggplot(aes(x = radcap, y = resid)) +
  geom_hline(yintercept = 0, lty=2) + 
  geom_point(color ="grey50") + 
  scale_x_continuous(breaks = c(0, 0.5,1,1.5,2), limits = c(0, 2))+
  scale_y_continuous(breaks = c(-2, -1.5, -1,-0.5, 0, 0.5, 1, 1.5, 2), limits = c(-2, 2)) +
  geom_smooth(aes(colour = radcap, fill = radcap), colour="black") +
  labs(y = "Pearson residuals", x =  "Outer measurement") +
  geom_text(aes(x = 0, y = 1.9, label="a)"),family="Times New Roman", colour="black", size=5)-> plot1

lm_out %>% #Pearson residuals against fitted
  augment(A2) %>% 
  mutate(resid = (.resid),
         fit = (.fitted)) %>% 
  ggplot(aes(x = fit, y = resid)) +
  geom_point(color ="grey50") + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2))+
  scale_y_continuous(breaks = c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2), limits = c(-2,2)) +
  geom_smooth(aes(colour = fit, fill = fit),colour="black") +
  geom_hline(yintercept = 0, lty=2) + 
  labs(y = "Pearson residuals", x =  "Fitted values") +
  geom_text(aes(x = 0, y = 2, label="b)", hjust = 1),family="Times New Roman", colour="black", size=5)-> plot2

#pf(0.49, 3, 17) # to figure out influence of cook's distance; relate D to the F(p, n-p) distribution and accertain the corresponding percentile value; see Neter et al. book
qf(0.5, 3+1, 211-3-1) # see https://rpubs.com/mpfoley73/460943 and Das and Gogoi 2015 pg 80 
qf(0.5, 4, 207) # p+1; n-p-1


lm_out %>% #Cook's distance plot
  augment(A2) %>% 
  mutate(cooksd = (.cooksd),
         count = 1:211,
         name= ifelse(cooksd >0.84, count, ""))%>% 
  ggplot(aes(x = count, y = cooksd, label=name)) +
  geom_bar(stat = "identity", colour = "grey50", 
           fill = "lightgrey",alpha=.7,
           width = 0.8, position = position_dodge(width = 0.2)) + 
  geom_text(size = 3, position = position_stack(vjust = 1.1)) + 
  scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5), limits = c(0,0.5)) +
  labs(y = "Cook's distance", x =  "Index") +
  geom_text(aes(x = 0, y = 0.5, label="c)"),family="Times New Roman", colour="black", size=5)-> plot3

# greater than 2 or three times p/n, it may be a concern (where p is the number of parameters (i.e., 2) and n is the number of observations (i.e., 51); p/n = 0.04 for this study; Dobson 2002). 
# 3/211 = 0.01; 0.01 *3 =0.03
lm_out %>% #leverage plot
  augment(A2) %>% 
  mutate(hat= (.hat),
         count = 1:211,
         name= ifelse(hat >0.03, count, "")) -> x%>% 
  ggplot(aes(x = count, y = hat, label=name)) +
  geom_bar(stat = "identity", colour = "grey50", 
           fill = "lightgrey",alpha=.7,
           width = 0.8, position = position_dodge(width = 0.2)) + 
  geom_text(size = 3, position = position_stack(vjust = 1.2)) + 
  scale_y_continuous(breaks = c(0, 0.05, 0.1, 0.15), limits = c(0,0.15)) +
  labs(y = "hat-values", x =  "Index") +
  geom_text(aes(x = 1, y = 0.15, label="d)"),family="Times New Roman", colour="black", size=5)-> plot4

lm_out %>% #Pearson by index
  augment(A2) %>% 
  mutate(resid = (.resid),
         count = 1:211) %>% 
  ggplot(aes(x = count, y = resid)) +
  geom_bar(stat = "identity", colour = "grey50", 
           fill = "lightgrey",alpha=.7,
           width = 0.8, position = position_dodge(width = 0.2)) + 
  scale_y_continuous(breaks = c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2), limits = c(-2,2)) +
  labs(y = "Pearson residuals", x =  "Index") +
  geom_text(aes(x = 0, y = 1.9, label="e)"),family="Times New Roman", colour="black", size=5)-> plot5

cowplot::plot_grid(plot1, plot2, plot3, plot4, plot5,  align = "vh", nrow = 3, ncol=2)
ggsave("figs/glm_diagnostics_cumulative_growth.png", dpi = 500, height = 6, width = 8, units = "in")



















