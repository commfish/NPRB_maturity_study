#notes---
#author: Sara E Miller 
#contact: sara.miller@alaska.gov; 907-465-4245
#Last edited: November 26, 2018

# load ---
library(extrafont)
library(tidyverse)
library(PropCIs)
library(MASS)

# load data--- 
data <- read.csv("data/data_11_26_2018.csv", check.names = FALSE) 

#confidence intervals for a difference in proportions (package PropCIs) (see package for refs)
#diffscoreci(x1, n1, x2, n2, conf.level)
#wald2ci(x1, n1, x2, n2, conf.level, adjust)
#prop.test (MASS package)
prop.test(x=c(122,99), n=c(600,400), correct=F)

diffscoreci(7, 21, 13, 17, conf.level=0.95)

wald2ci(7, 21, 13, 17, conf.level=0.95, adjust = "Wald")