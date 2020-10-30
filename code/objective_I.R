#author: Sara E Miller 
#contact: sara.miller@alaska.gov; 907-465-4245
#Last edited: February, 2020
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4885900/
#https://stats.stackexchange.com/questions/82105/mcfaddens-pseudo-r2-interpretation
#S:\Region1Shared-DCF\Research\Herring-Dive Fisheries\Herring\NPRB-Herring Maturity Study\datasheets\datasheets\MTA Lab Datasheets\fall field study-MTA lab measurements\remeasurements_corrected data\ MTA Lab data_length_weight_Jan_2019
#matched to histology data sheet (S:\Region1Shared-DCF\Research\Herring-Dive Fisheries\Herring\NPRB-Herring Maturity Study\datasheets\datasheets\Histology-Pathologist report\Alaska Department of Fish and Game 10_17_2018_corrected)
#the following code then takes the histology dataset data_11_26_2018 and matches it to the MTA lab increment data by image name
#https://stackoverflow.com/questions/44993844/how-to-report-likelihood-ratio-test-results
#https://stats.stackexchange.com/questions/59879/logistic-regression-anova-chi-square-test-vs-significance-of-coefficients-ano
#http://ww2.coastal.edu/kingw/statistics/R-tutorials/logistic.html

# load libraries----
# devtools::install_github("ben-williams/FNGr")
source("code/functions.R")
library(extrafont)
loadfonts(device="win")
#library(nlme)
library(tidyverse)
library(fngr)
library(broom)
library(cowplot)
#library(mixtools)
library(car)
library(modEvA) #pseudo-R squared
#library(MASS)
#library(digest)
library(ResourceSelection)
library("ggpubr")
library("devtools")
devtools::install_github("commfish/fngr")
library(fngr)
theme_set(theme_report(base_size = 12))

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

# delete a or b at the end of the image name (this just depicts what image was used)
increment$image_name = gsub("a","",increment$image_name)

# match increments to awl data
merge <-merge(increment, data_clean, all.y=T)
write.csv(merge, "data/check.csv") 
merge %>% 
  filter(!(Increment1 %in% c("", NA)))%>% 
  filter(!(scale_region%in% c("OOA (OUT OF AREA)", "A", "C", "D", "E", "G", "H"))) %>% 
  filter(!(sex_histology %in% c("Male"))) %>% 
  mutate(maturity = ifelse(maturation_status_histology == 1, 'immature', 'mature')) -> clean_dataset #n=211 useable samples

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
# radcap is cumulative growth
sample1 %>%
  mutate(anu_adj = ifelse(anu1>0 & anu2==0 & anu3 ==0 & anu4 ==0 & anu5==0 & anu6==0 & anu7==0, anu1,
                          ifelse(anu1>0 & anu2>0& anu3 ==0 & anu4 ==0 & anu5==0 & anu6==0 & anu7==0, anu2,
                                 ifelse(anu1>0 & anu2>0 & anu3>0 & anu4 ==0 & anu5==0 & anu6==0 & anu7==0, anu3,
                                        ifelse(anu1>0 & anu2>0 & anu3>0 & anu4>0 & anu5==0 & anu6==0 & anu7==0, anu4,
                                               ifelse(anu1>0 & anu2>0 & anu3>0 & anu4>0 & anu5>0 & anu6==0 & anu7==0, anu5,
                                                      ifelse(anu1>0 & anu2>0 & anu3>0 & anu4>0 & anu5>0 & anu6>0 & anu7==0, anu6,anu7))))))) %>%
  filter(sex_histology == "Female") %>% 
  #mutate(age = as.factor(age)) %>% 
  #filter(sample_no != 743) %>% 
  dplyr::select(sample_no, image_name,year, age, sex_histology, maturation_status_histology, maturity, anu_adj, radcap, length_mm, Increment1, Increment2, 
                Increment3, Increment4, Increment5, Increment6, Increment7, anu1, anu2, anu3, anu4, anu5, anu6, anu7, histology_notes) -> merge
write.csv(merge, "data/obj1_data.csv") 

merge %>% 
  mutate (scale_length_prior = radcap - anu_adj) %>% 
dplyr::select(sample_no, age, maturation_status_histology, maturity, anu_adj, radcap, length_mm, scale_length_prior) -> merge
write.csv(merge, "data/obj1_data2.csv") 

merge %>% 
  dplyr::select(maturity, age, length_mm, anu_adj, scale_length_prior) %>% 
  group_by(maturity,age) %>% 
  summarize(n=n(),
            mean_anu_adj =mean(anu_adj),
            stdv_anu_adj = sd(anu_adj),
            min_anu_adj =min(anu_adj),
            max_anu_adj = max(anu_adj),
            mean_length =mean(length_mm),
            stdv_length = sd(length_mm),
            min_length =min(length_mm),
            max_length = max(length_mm),
            mean_length_prior =mean(scale_length_prior),
            stdv_length_prior = sd(scale_length_prior),
            min_length_prior =min(scale_length_prior),
            max_length_prior = max(scale_length_prior)) -> summary_table 


# Generalized Linear models (age two)----
## Outer Increment Measurement
## https://stackoverflow.com/questions/23826621/plotting-glm-interactions-newdata-structure-in-predict-function
merge %>%
  mutate(age = as.factor(age)) %>%
  dplyr::mutate(maturity = ifelse(maturity == "mature", 1 , 0)) -> merge_dataset

merge_dataset %>%
  filter(age == 2) %>%
  as.data.frame()-> merge_dataset_two

merge_dataset %>%
  filter(age == 3) -> merge_dataset_three

fit <- glm(maturity ~ (anu_adj * scale_length_prior) , family = binomial(link=logit), data = merge_dataset_two) 
summary(fit)
vif(fit) # car package
anova(fit, test = "Chisq") #http://ww2.coastal.edu/kingw/statistics/R-tutorials/logistic.html
Anova(fit)
1 - pchisq(63.376-60.324, df=3)
RsqGLM(fit)#peudo R2 
summary(fit)
hoslem.test(merge_dataset_two$maturity, fitted(fit)) #goodness of fit test (ResourceSelection package); https://www.theanalysisfactor.com/r-glm-model-fit/
exp(coef(fit)) #https://stats.idre.ucla.edu/r/dae/logit-regression/ ; odds ratio
exp(cbind(OR = coef(fit), confint(fit))) # odds ratio and 95% confidence intervals

merge_dataset_two %>% 
  do(A2 = glm(maturity ~ anu_adj * scale_length_prior, data = merge_dataset_two, family = binomial(link=logit))) -> lm_out
head(augment(lm_out, merge_dataset_two, type.residuals="pearson"))
outlierTest(fit) #Bonferroni p-values
residualPlot(fit, variable = "fitted", type = "pearson",
             plot = TRUE, quadratic = FALSE, smooth=TRUE) #curvature test
marginalModelPlots(fit) #marginal model plots
mmp(fit, merge_dataset_two$anu_adj, xlab="outer proportion", ylab="maturity" , family="A")


lm_out %>% 
  tidy(A2) %>% 
  mutate(model = "fit_age2") -> A2
write_csv(A2, "output/lm.csv")

lm_out %>% 
  glance(A2) %>% 
  mutate(model = "fit_age2") -> A2
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
  scale_x_continuous(breaks = c(-4, -3, -2, -1, 0, 1), limits = c(-4, 1))+
  scale_y_continuous(breaks = c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2), limits = c(-2,2)) +
  geom_smooth(aes(colour = fit, fill = fit),colour="black") +
  geom_hline(yintercept = 0, lty=2) + 
  labs(y = "Pearson residuals", x =  "Fitted values") +
  geom_text(aes(x = -4, y = 2, label="b)", hjust = 1),family="Times New Roman", colour="black", size=5)-> plot2

#pf(0.49, 3, 17) # to figure out influence of cook's distance; relate D to the F(p, n-p) distribution and accertain the corresponding percentile value; see Neter et al. book
qf(0.5, 3, 64-3-1) # see https://rpubs.com/mpfoley73/460943 and Das and Gogoi 2015 pg 80 
qf(0.5, 3, 200) # p+1; n-p-1


lm_out %>% #Cook's distance plot
  augment(A2) %>% 
  mutate(cooksd = (.cooksd),
         count = 1:68,
         name= ifelse(cooksd >0.79, count, "")) %>% 
  ggplot(aes(x = count, y = cooksd, label=name)) +
  geom_bar(stat = "identity", colour = "grey50", 
           fill = "lightgrey",alpha=.7,
           width = 0.8, position = position_dodge(width = 0.2)) + 
  geom_text(size = 3, position = position_stack(vjust = 1.1)) + 
  scale_y_continuous(breaks = c(0,0.05, 0.10, 0.15,0.20), limits = c(0,0.20)) +
  labs(y = "Cook's distance", x =  "Index") +
  geom_text(aes(x = 0, y = 0.19, label="c)"),family="Times New Roman", colour="black", size=5)-> plot3

# greater than 2 or three times p/n, it may be a concern (where p is the number of parameters (i.e., 2) and n is the number of observations (i.e., 51); p/n = 0.04 for this study; Dobson 2002). 
# 3/64 = 0.05; 0.05 *3 =0.14
lm_out %>% #leverage plot
  augment(A2) %>% 
  mutate(hat= (.hat),
         count = 1:68,
         name= ifelse(hat >0.14, count, "")) %>% 
  ggplot(aes(x = count, y = hat, label=name)) +
  geom_bar(stat = "identity", colour = "grey50", 
           fill = "lightgrey",alpha=.7,
           width = 0.8, position = position_dodge(width = 0.2)) + 
  geom_text(size = 3, position = position_stack(vjust = 1.2)) + 
  scale_y_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25), limits = c(0, 0.25)) +
  labs(y = "hat-values", x =  "Index") +
  geom_text(aes(x = 1, y = 0.25, label="d)"),family="Times New Roman", colour="black", size=5)-> plot4

lm_out %>% #Pearson by index
  augment(A2) %>% 
  mutate(resid = (.resid),
         count = 1:68) %>% 
  ggplot(aes(x = count, y = resid)) +
  geom_bar(stat = "identity", colour = "grey50", 
           fill = "lightgrey",alpha=.7,
           width = 0.8, position = position_dodge(width = 0.2)) + 
  scale_y_continuous(breaks = c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5), limits = c(-2, 2.5)) +
  labs(y = "Pearson residuals", x =  "Index") +
  geom_text(aes(x = 0, y = 2.5, label="e)"),family="Times New Roman", colour="black", size=5)-> plot5

cowplot::plot_grid(plot1, plot2, plot3, plot4, plot5,  align = "vh", nrow = 3, ncol=2)
ggsave("figs/glm_diagnostics_outer_increment_age2.png", dpi = 500, height = 6, width = 8, units = "in")

# prediction plot; example: https://www.theanalysisfactor.com/generalized-linear-models-r-graphs/
fit_anu_adj <- glm(maturity ~ (anu_adj) , family = binomial, data = merge_dataset_two) 
fit_scale_length_prior <- glm(maturity ~ (scale_length_prior) , family = binomial, data = merge_dataset_two) 
range(merge_dataset_two$anu_adj)
range(merge_dataset_two$scale_length_prior)
xanu_adj <-seq(0,2,0.01)
yanu_adj <- predict(fit_anu_adj, list(anu_adj=xanu_adj),type="response")
plot(merge_dataset_two$anu_adj, merge_dataset_two$maturity, pch = 16, xlab = "ANU_ADJ", ylab = "MATURITY")
lines(xanu_adj, yanu_adj, col = "red", lwd = 2)

xscale_length_prior <-seq(0,3,0.01)
yscale_length_prior <- predict(fit_scale_length_prior, list(scale_length_prior=xscale_length_prior),type="response")
plot(merge_dataset_two$scale_length_prior, merge_dataset_two$maturity, pch = 16, xlab = "SCALE_LENGTH_PRIOR", ylab = "MATURITY")
lines(xscale_length_prior, yscale_length_prior, col = "red", lwd = 2)

# Generalized Linear models (age three)----
## Outer Increment Measurement
fit <- glm(maturity ~ (anu_adj * scale_length_prior) , family = binomial(link=logit), data = merge_dataset_three) 
vif(fit) # car package
summary(fit)
anova(fit, test = "Chisq") #http://ww2.coastal.edu/kingw/statistics/R-tutorials/logistic.html
Anova(fit)
74.786-74.084
1 - pchisq(74.786-74.084, df=3)
1 - pchisq(0.41+0.26+0.03, df=3)
RsqGLM(fit)#peudo R2 
summary(fit)
hoslem.test(merge_dataset_three$maturity, fitted(fit)) #goodness of fit test (ResourceSelection package); https://www.theanalysisfactor.com/r-glm-model-fit/
exp(coef(fit)) #https://stats.idre.ucla.edu/r/dae/logit-regression/ ; odds ratio
exp(cbind(OR = coef(fit), confint(fit))) # odds ratio and 95% confidence intervals

merge_dataset_three %>% 
  do(A2 = glm(maturity ~ anu_adj * scale_length_prior, data = merge_dataset_three, family = binomial(link=logit))) -> lm_out
head(augment(lm_out, merge_dataset_three, type.residuals="pearson"))
outlierTest(fit) #Bonferroni p-values
residualPlot(fit, variable = "fitted", type = "pearson",
             plot = TRUE, quadratic = FALSE, smooth=TRUE) #curvature test
marginalModelPlots(fit) #marginal model plots
mmp(fit, merge_dataset_three$anu_adj, xlab="outer proportion", ylab="maturity" , family="A")


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
  scale_x_continuous(breaks = c(-1, 0, 1), limits = c(-1, 1))+
  scale_y_continuous(breaks = c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2), limits = c(-2,2)) +
  geom_smooth(aes(colour = fit, fill = fit),colour="black") +
  geom_hline(yintercept = 0, lty=2) + 
  labs(y = "Pearson residuals", x =  "Fitted values") +
  geom_text(aes(x = -4, y = 2, label="b)", hjust = 1),family="Times New Roman", colour="black", size=5)-> plot2

#pf(0.49, 3, 17) # to figure out influence of cook's distance; relate D to the F(p, n-p) distribution and accertain the corresponding percentile value; see Neter et al. book
qf(0.5, 3, 64-3-1) # see https://rpubs.com/mpfoley73/460943 and Das and Gogoi 2015 pg 80 
qf(0.5, 3, 200) # p+1; n-p-1


lm_out %>% #Cook's distance plot
  augment(A2) %>% 
  mutate(cooksd = (.cooksd),
         count = 1:54,
         name= ifelse(cooksd >0.79, count, "")) %>% 
  ggplot(aes(x = count, y = cooksd, label=name)) +
  geom_bar(stat = "identity", colour = "grey50", 
           fill = "lightgrey",alpha=.7,
           width = 0.8, position = position_dodge(width = 0.2)) + 
  geom_text(size = 3, position = position_stack(vjust = 1.1)) + 
  scale_y_continuous(breaks = c(0,0.05, 0.10, 0.15,0.20), limits = c(0,0.20)) +
  labs(y = "Cook's distance", x =  "Index") +
  geom_text(aes(x = 0, y = 0.19, label="c)"),family="Times New Roman", colour="black", size=5)-> plot3

# greater than 2 or three times p/n, it may be a concern (where p is the number of parameters (i.e., 2) and n is the number of observations (i.e., 51); p/n = 0.04 for this study; Dobson 2002). 
# 3/64 = 0.05; 0.05 *3 =0.14
lm_out %>% #leverage plot
  augment(A2) %>% 
  mutate(hat= (.hat),
         count = 1:54,
         name= ifelse(hat >0.14, count, "")) %>% 
  ggplot(aes(x = count, y = hat, label=name)) +
  geom_bar(stat = "identity", colour = "grey50", 
           fill = "lightgrey",alpha=.7,
           width = 0.8, position = position_dodge(width = 0.2)) + 
  geom_text(size = 3, position = position_stack(vjust = 1.2)) + 
  scale_y_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25), limits = c(0, 0.25)) +
  labs(y = "hat-values", x =  "Index") +
  geom_text(aes(x = 1, y = 0.25, label="d)"),family="Times New Roman", colour="black", size=5)-> plot4

lm_out %>% #Pearson by index
  augment(A2) %>% 
  mutate(resid = (.resid),
         count = 1:54) %>% 
  ggplot(aes(x = count, y = resid)) +
  geom_bar(stat = "identity", colour = "grey50", 
           fill = "lightgrey",alpha=.7,
           width = 0.8, position = position_dodge(width = 0.2)) + 
  scale_y_continuous(breaks = c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5), limits = c(-2, 2.5)) +
  labs(y = "Pearson residuals", x =  "Index") +
  geom_text(aes(x = 0, y = 2.5, label="e)"),family="Times New Roman", colour="black", size=5)-> plot5

cowplot::plot_grid(plot1, plot2, plot3, plot4, plot5,  align = "vh", nrow = 3, ncol=2)
ggsave("figs/glm_diagnostics_outer_increment_age3.png", dpi = 500, height = 6, width = 8, units = "in")

# prediction plot; example: https://www.theanalysisfactor.com/generalized-linear-models-r-graphs/
fit_anu_adj <- glm(maturity ~ (anu_adj) , family = binomial, data = merge_dataset_three) 
fit_scale_length_prior <- glm(maturity ~ (scale_length_prior) , family = binomial, data = merge_dataset_three) 
range(merge_dataset_three$anu_adj)
range(merge_dataset_three$scale_length_prior)
xanu_adj <-seq(0,1,0.01)
yanu_adj <- predict(fit_anu_adj, list(anu_adj=xanu_adj),type="response")
plot(merge_dataset_three$anu_adj, merge_dataset_three$maturity, pch = 16, xlab = "ANU_ADJ", ylab = "MATURITY")
lines(xanu_adj, yanu_adj, col = "red", lwd = 2)

xscale_length_prior <-seq(0,5,0.01)
yscale_length_prior <- predict(fit_scale_length_prior, list(scale_length_prior=xscale_length_prior),type="response")
plot(merge_dataset_three$scale_length_prior, merge_dataset_three$maturity, pch = 16, xlab = "SCALE_LENGTH_PRIOR", ylab = "MATURITY")
lines(xscale_length_prior, yscale_length_prior, col = "red", lwd = 2)


# Exploratory Plots----
# NOTE: order of variables is wrong
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
eda.norm(as.numeric(age2immature$anu_adj))
eda.norm(as.numeric(age2mature$anu_adj))
eda.norm(as.numeric(age3immature$anu_adj))
eda.norm(as.numeric(age3mature$anu_adj))

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

# Histograms by age----
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
## outer incremment growth
merge %>% 
  mutate(age=as.numeric(age))%>% 
  ggplot(data=., aes(x = age, y = anu_adj, color = maturity, shape = maturity)) +
  geom_jitter(size = 1) +
  labs(x = "Age",y = "Outer increment measurement (mm)") + 
  theme(legend.title=element_blank(), legend.position=c(.8,.9), legend.text=element_text(size=12)) +
  scale_color_manual(values=c("#999999", "black")) +
  scale_fill_manual(values=c("#999999", "black" )) +
  scale_shape_manual(values=c(8, 16)) +
  scale_x_continuous(breaks = c(2, 3, 4, 5, 6), limits = c(1.5,6.5)) +
  annotate("text",x = 1.5, y=2, label="A)", family="Times" ,size=4) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0), limits = c(0, 2.0))-> plot1

merge %>%
  as.data.frame()%>% 
  mutate(age=as.factor(age)) %>% 
  ggplot(data=.,aes(x = age, y = anu_adj)) + labs(x = "Age",y = "Outer increment measurement (mm)") + 
  geom_boxplot(aes(fill = maturity)) +
  scale_color_manual(values=c("#999999", "black")) + theme(legend.title=element_blank(), legend.position=c(.8,.9), legend.text=element_text(size=12)) +
  scale_fill_manual(values=c("#999999", "black" )) + 
  scale_y_continuous(breaks = c(0, 0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_x_continuous(breaks = c(2, 3, 4, 5, 6), limits = c(1.5,6.5)) +
  annotate("text",x = 0.75, y=2, label="B)", family="Times" ,size=4)-> plot2

merge %>% 
  ggplot(data=., aes(x = age, y = radcap, color = maturity, shape = maturity)) +
  geom_jitter(size = 1) +
  labs(x = "Age",y = "Cumulative growth (mm)") + 
  theme(legend.title=element_blank(), legend.position=c(.8,.9), legend.text=element_text(size=12)) +
  scale_color_manual(values=c("#999999", "black")) +
  scale_fill_manual(values=c("#999999", "black" )) +
  scale_shape_manual(values=c(8, 16)) +
  annotate("text",x = 2, y=10, label="B)", family="Arial", colour="black", size=4) +
  scale_y_continuous(breaks = c(0, 1,2,3,4,5,6,7,8,9,10), limits = c(0, 10))-> plot3

cowplot::plot_grid(plot1, plot2,  align = "vh", nrow = 1, ncol=2)
ggsave("figs/scatterplot_obj1.png", dpi = 500, height = 4, width = 8, units = "in")

# Scatterplot Figures---- 
## growth up into the summer 
merge %>% 
  mutate(age=as.numeric(age))%>% 
  ggplot(data=., aes(x = age, y = scale_length_prior, color = maturity, shape = maturity)) +
  geom_jitter(size = 1) +
  labs(x = "Age",y = "Scale measurement up to the last increment (mm)") + 
  theme(legend.title=element_blank(), legend.position=c(.86,.2), legend.text=element_text(size=12)) +
  scale_color_manual(values=c("#999999", "black")) +
  scale_fill_manual(values=c("#999999", "black" )) +
  scale_shape_manual(values=c(8, 16)) +
  scale_x_continuous(breaks = c(2, 3, 4, 5, 6), limits = c(1.5,6.5)) +
  annotate("text",x = 1.5, y=6, label="A)", family="Times" ,size=4) +
  scale_y_continuous(breaks = c(0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0), limits = c(0, 6.0)) -> plot1

merge %>%
  as.data.frame()%>% 
  mutate(age=as.factor(age)) %>% 
  ggplot(data=.,aes(x = age, y = scale_length_prior)) + labs(x = "Age",y = "Scale measurement up to the last increment (mm)") + 
  geom_boxplot(aes(fill = maturity)) +
  scale_color_manual(values=c("#999999", "black")) + theme(legend.title=element_blank(), legend.position=c(.86,.2), legend.text=element_text(size=12)) +
  scale_fill_manual(values=c("#999999", "black" )) + 
  scale_y_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6), limits = c(0, 6)) +
  #scale_x_continuous(breaks = c(2, 3, 4, 5, 6), limits = c(1.5,6.5)) +
  annotate("text",x = 0.75, y=6, label="B)", family="Times" ,size=4) -> plot2

cowplot::plot_grid(plot1, plot2,  align = "vh", nrow = 1, ncol=2)
ggsave("figs/scatterplot__scale_length.png", dpi = 500, height = 4, width = 8, units = "in")


merge %>%
  as.data.frame()%>% 
  mutate(age=as.factor(age)) %>% 
  ggplot(data=.,aes(x = age, y = length_mm)) + labs(x = "Age",y = "Length at capture (mm)") + 
  geom_boxplot(aes(fill = maturity)) +
  scale_color_manual(values=c("#999999", "black")) + theme(legend.title=element_blank(), legend.position=c(.15,.15), legend.text=element_text(size=12)) +
  scale_fill_manual(values=c("#999999", "black" )) + 
  scale_y_continuous(breaks = c(0,50,100,150,200,250), limits = c(0, 250))-> plot1
cowplot::plot_grid(plot1,  align = "vh", nrow = 1, ncol=1)
ggsave("figs/length_at_capture.png", dpi = 500, height = 4, width = 8, units = "in")

# Correlation Analysis---
# age versus length at capture
merge$age<-as.numeric(merge$age)
ggscatter(merge, x = "age", y = "length_mm",  
          conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Age", ylab = "Length at capture (mm)") -> plot1
ggscatter(merge, y = "length_mm", x = "radcap", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Total Scale Length", ylab = "Length at capture (mm)") -> plot2
ggscatter(merge, x = "age", y = "radcap", 
          conf.int = TRUE, add="loess", 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Age", ylab = "Total Scale Length")+ annotate("text",x = 6, y=7, label="B)", family="Arial" ,size=4) -> plot3
cowplot::plot_grid(plot1,   align = "vh", nrow = 1, ncol=1)
ggsave("figs/correlation.png", dpi = 500, height = 4, width = 8, units = "in")

fit <- lm(length_mm ~ (radcap), data = merge) 
fit1 <- lm(length_mm ~ poly(radcap,3), data = merge) 

merge %>% #GSI data not in data set 
  filter (age =="3") %>%
  ggscatter(., x = "anu_adj", y = "GSI", shape= "maturity", color = "maturity",
            conf.int = TRUE, add="loess",  palette=c("#999999", "black"), 
            cor.coef = FALSE, cor.method = "spearman",
            xlab = "Outer increment measurement (mm)", ylab = "GSI")-> plot1
cowplot::plot_grid(plot1,   align = "vh", nrow = 1, ncol=1)
#ggsave("figs/GSI_corr.png", dpi = 500, height = 4, width = 8, units = "in")

merge %>% #GSI data not in data set
  filter (age =="3") %>%
  filter (maturity =="mature") %>%
  ggscatter(., x = "anu_adj", y = "GSI", shape= "maturity", color = "maturity",
            conf.int = TRUE, add="loess",  palette=c("#999999", "black"), 
            cor.coef = TRUE, cor.method = "spearman",
            xlab = "Outer increment measurement (mm)", ylab = "GSI")-> plot1
cowplot::plot_grid(plot1,   align = "vh", nrow = 1, ncol=1)
#ggsave("figs/GSI_corr.png", dpi = 500, height = 4, width = 8, units = "in")

par(mfrow=c(1,1))
with(merge,interaction.plot(age, maturity, anu_adj, type="b", pch=c(1,16), ylab ="outer incremental length (mm)",
                            xlab="Age"))
