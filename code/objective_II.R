# notes---
# author: Sara E Miller 
# contact: sara.miller@alaska.gov; 907-465-4245
# Last edited: December, 2018

# load ----
source("code/helper.r")

# load data ---- 
read.csv("data/data_11_26_2018.csv", check.names = FALSE) %>% 
  filter(!(maturation_status_histology %in% c("", "NA", "no slide", "no score", "4-1", "3-4", "2-3", "1-2")))  -> data_clean

data_clean %>%
  group_by(age, maturation_status_histology) %>%
  summarise(count = n()) -> table1 

data_clean %>%
  group_by(maturity_state_GSI, maturation_status_histology) %>%
  summarise(count = n()) -> table2 


data_clean %>%
  group_by(maturity_state_field, maturation_status_histology) %>%
  summarise(count = n()) -> table3 

data %>% 
  filter(!(sex_histology %in% c("", "NA", "no slide"))) -> data_clean

data_clean %>%
  mutate (date = mdy(catch_date)) %>%
  group_by(maturity_state_GSI, date, sex_histology) %>%
  summarise(count = n(), 
            mean = mean (GSI)) -> table4 

#figures ----


#bubble plot of age versus histology
bb <- c(201,100,50,20,10,5) # define breaks.
ll <- c("200+","100","50","20","10", "5") # labels.

ggplot(table1) + 
  geom_point(aes(age, maturation_status_histology, size = count), shape=16, alpha=0.80) +
  labs(x="Age", y = "Maturation Status (histology)\n") +  
  scale_size_continuous(name = "",
                        breaks = bb,
                        limits = c(0, 250),
                        labels = ll,
                        range = c(0, 10))

ggsave(filename = 'figs/bubbleplot.png', dpi =200, width=6, height=8, units = "in")



# confidence intervals for a difference in proportions (package PropCIs) (see package for refs)
# diffscoreci(x1, n1, x2, n2, conf.level)
# wald2ci(x1, n1, x2, n2, conf.level, adjust)
# prop.test (MASS package)
prop.test(x=c(122,99), n=c(600,400), correct=F)

diffscoreci(7, 21, 13, 17, conf.level=0.95)

wald2ci(7, 21, 13, 17, conf.level=0.95, adjust = "Wald")
