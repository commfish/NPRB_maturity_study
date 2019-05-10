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

theme_set(theme_sleek())

asintrans <- function(p) {asin(sqrt(p))} # arctransformation function

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
  summarize(n=n())-> x #sample sizes for study

#calculate proportion of outer scale ring
clean_dataset %>% 
  dplyr::select(age,maturity_state_histology, scale_region, sex_histology, Increment1, Increment2, Increment3, 
                Increment4, Increment5, Increment6, Increment7, image_name, length_mm, weight_grams, 
                scale_condition_code, scale_region, scale_shape, sex_histology, maturation_status_histology) %>%
  gather(value, variable, Increment1:Increment7) %>%
  filter(!is.na(variable)) %>%
  spread (key = value, value = variable)->s1

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
         rad7 = sum(anu1, anu2, anu3, anu4, anu5, anu6, anu7),
         prop1 = rad1 / radcap,
         prop2 = rad2 / radcap,
         prop3 = rad3 / radcap,
         prop4 = rad4 / radcap,
         prop5 = rad5 / radcap,
         prop6 = rad6 / radcap,
         prop7 = rad7 / radcap) %>%
  mutate(prop1_adj = prop1,
         prop2_adj = ifelse(prop2==1, NA, prop2),
         prop3_adj = ifelse(prop3==1, NA, prop3),
         prop4_adj = ifelse(prop4==1, NA, prop4),
         prop5_adj = ifelse(prop5==1, NA, prop5),
         prop6_adj = ifelse(prop6==1, NA, prop6),
         prop7_adj = ifelse(prop7==1, NA, prop7)) %>%
  rowwise() %>% 
  mutate(max=max(prop1_adj, prop2_adj, prop3_adj, prop4_adj, prop5_adj, prop6_adj, prop7_adj))->sample1

sample1 %>%
  gather(value, variable, prop1_adj:prop7_adj) %>% 
  group_by(image_name) %>% 
  filter(!is.na(variable)) %>%
  summarize(max=max(variable))%>%
  mutate(outer_prop=1-max)-> sample2 #outer prop

merge<- merge(x = sample1, y= sample2, by=c("image_name"), all.x  = T)

#Histograms of outer ring 
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
merge %>%
  filter(age == 7) -> age7

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

eda.norm(as.numeric(age2$outer_prop))
eda.norm(as.numeric(age3$outer_prop))
eda.norm(as.numeric(age4$outer_prop))
eda.norm(as.numeric(age5$outer_prop))
eda.norm(as.numeric(age6$outer_prop))
eda.norm(as.numeric(age7$outer_prop))

tickr_length <- data.frame(outer_prop = 0:1)
axisb <- tickr(tickr_length, outer_prop, 0.2)
 
ggplot(age2, aes(x=outer_prop)) + geom_histogram(alpha=0.5, position = 'identity') +
  ylab("Frequency")+ xlab("Outer increment proportion") +
  scale_x_continuous(breaks = axisb$breaks, labels = axisb$labels) +
  ggtitle("Age 2; n=215") +
  facet_wrap(~maturity_state_histology)

x3<-ggplot(Age_3, aes(x=value, color=variable, fill=variable)) +
  geom_density(alpha=0.5, adjust=1)+
  labs(x="Scale increment length (mm)", y="Density")
x3<-x3+ggtitle("Age 3")+theme_set(theme_bw(base_size=14,base_family=
                                             'Times New Roman')+
                                    theme(panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank()))
x3<-x3 + theme(legend.position="none")+scale_x_continuous(breaks=seq(0,max(Age_3$value+0.1, na.rm=TRUE),0.2))


Age_4$Year<-as.factor(Age_4$Year)
cdat <- ddply(Age_4, "Year", summarise, mean.Year=mean(value))
y<-ggplot(Age_4, aes(x=value)) + geom_histogram(alpha=0.5, position = 'identity')+
  ylab("Frequency")+ xlab("Scale increment length (mm)")+
  scale_x_continuous(breaks=seq(0,max(Age_4$value, na.rm=TRUE),0.2) )                       
y<-y+geom_vline(data=cdat, aes(xintercept=mean.Year, colour=Year),   # Ignore NA values for mean
                linetype="dashed", size=1)
y<-y+ggtitle("Age 4; n=1229")+theme_set(theme_bw(base_size=14,base_family=
                                                   'Times New Roman')+
                                          theme(panel.grid.major = element_blank(),
                                                panel.grid.minor = element_blank()))

x4<-ggplot(Age_4, aes(x=value, color=variable, fill=variable)) +
  geom_density(alpha=0.5, adjust=1)+
  labs(x="Scale increment length (mm)", y="Density")
x4<-x4+ggtitle("Age 4")+theme_set(theme_bw(base_size=14,base_family=
                                             'Times New Roman')+
                                    theme(panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank()))
x4<-x4 + theme(legend.position="none")+scale_x_continuous(breaks=seq(0,max(Age_4$value+0.1, na.rm=TRUE),0.2))



Age_5$Year<-as.factor(Age_5$Year)
cdat <- ddply(Age_5, "Year", summarise, mean.Year=mean(value))
z<-ggplot(Age_5, aes(value)) + geom_histogram(alpha=0.5, position = 'identity')+
  ylab("Frequency")+ xlab("Scale increment length (mm)")+
  scale_x_continuous(breaks=seq(0,max(Age_5$value, na.rm=TRUE),0.1) )                       
z<-z+geom_vline(data=cdat, aes(xintercept=mean.Year, colour=Year),   # Ignore NA values for mean
                linetype="dashed", size=1)
z<-z+ggtitle("Age 5; n=469")+theme_set(theme_bw(base_size=14,base_family=
                                                  'Times New Roman')+
                                         theme(panel.grid.major = element_blank(),
                                               panel.grid.minor = element_blank()))
x5<-ggplot(Age_5, aes(x=value, color=variable, fill=variable)) +
  geom_density(alpha=0.5, adjust=1)+
  labs(x="Scale increment length (mm)", y="Density")
x5<-x5+ggtitle("Age 5")+theme_set(theme_bw(base_size=14,base_family=
                                             'Times New Roman')+
                                    theme(panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank()))
x5<-x5 + theme(legend.position="none")+scale_x_continuous(breaks=seq(0,max(Age_6$value+0.1, na.rm=TRUE),0.2))

Age_6$Year<-as.factor(Age_6$Year)
cdat <- ddply(Age_6, "Year", summarise, mean.Year=mean(value))
a<-ggplot(Age_6, aes(value)) + geom_histogram(alpha=0.5, position = 'identity')+
  ylab("Frequency")+ xlab("Scale increment length (mm)")+
  scale_x_continuous(breaks=seq(0,max(Age_6$value, na.rm=TRUE),0.1) )                       
a<-a+geom_vline(data=cdat, aes(xintercept=mean.Year, colour=Year),   # Ignore NA values for mean
                linetype="dashed", size=1)
a<-a+ggtitle("Age 6; n=447")+theme_set(theme_bw(base_size=14,base_family=
                                                  'Times New Roman')+
                                         theme(panel.grid.major = element_blank(),
                                               panel.grid.minor = element_blank()))
x6<-ggplot(Age_6, aes(x=value, color=variable, fill=variable)) +
  geom_density(alpha=0.5, adjust=1)+
  labs(x="Scale increment length (mm)", y="Density")
x6<-x6+ggtitle("Age 6")+theme_set(theme_bw(base_size=14,base_family=
                                             'Times New Roman')+
                                    theme(panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank()))
x6<-x6 + theme(legend.position="none")+scale_x_continuous(breaks=seq(0,max(Age_6$value+0.1, na.rm=TRUE),0.2))


Age_7$Year<-as.factor(Age_7$Year)
cdat <- ddply(Age_7, "Year", summarise, mean.Year=mean(value))
b<-ggplot(Age_7, aes(value)) + geom_histogram(alpha=0.5, position = 'identity')+
  ylab("Frequency")+ xlab("Scale increment length (mm)")+
  scale_x_continuous(breaks=seq(0,max(Age_7$value, na.rm=TRUE),0.1) ) 
b<-b+geom_vline(data=cdat, aes(xintercept=mean.Year, colour=Year),   # Ignore NA values for mean
                linetype="dashed", size=1)
b<-b+ggtitle("Age 7; n=446")+theme_set(theme_bw(base_size=14,base_family=
                                                  'Times New Roman')+
                                         theme(panel.grid.major = element_blank(),
                                               panel.grid.minor = element_blank()))
x7<-ggplot(Age_7, aes(x=value, color=variable, fill=variable)) +
  geom_density(alpha=0.5, adjust=1)+
  labs(x="Scale increment length (mm)", y="Density")
x7<-x7+ggtitle("Age 7")+theme_set(theme_bw(base_size=14,base_family=
                                             'Times New Roman')+
                                    theme(panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank()))
x7<-x7 + theme(legend.position="none")+scale_x_continuous(breaks=seq(0,max(Age_7$value+0.1, na.rm=TRUE),0.2))

library(grid)
g6<-grid.arrange(x,x3,y,x4,z,x5,a,x6,b,x7, ncol=2)
dev.off()

#FIGURE 4 (matching age-3 2011 with age-8 2016) #mixture models
Data <- read.csv("scale.csv", sep=",", header = TRUE, check.names = FALSE)  #Read in csv file, subset the spawning stock
Dataset <- melt(Data, id=c("Year", "AGE", "SEX","LENGTH"), na.rm=TRUE)
#create dataset for 2011 with the last ring
#AGE 3 DATASET
Age_3LRa<-subset(Dataset, Dataset$AGE==3& Dataset$Year==2011&Dataset$variable=='CSW3') 
Age_8_3RING<-subset(Dataset, Dataset$AGE==8& Dataset$Year==2016&Dataset$variable=='CSW3') 
Age_3LRa["Ring"] <- "Last Ring"
Age_8_3RING["Ring"] <- "Age Ring"
Dataset_Age3<- rbind(Age_3LRa,Age_8_3RING) 
#AGE 4 DATASET
Age_4LRa<-subset(Dataset, Dataset$AGE==4& Dataset$Year==2011&Dataset$variable=='CSW4') 
Age_4LRb<-subset(Dataset, Dataset$AGE==4& Dataset$Year==2012&Dataset$variable=='CSW4') 
Age_8_4RINGa<-subset(Dataset, Dataset$AGE==8& Dataset$Year==2015&Dataset$variable=='CSW4')
Age_8_4RINGb<-subset(Dataset, Dataset$AGE==8& Dataset$Year==2016&Dataset$variable=='CSW4')
Age_9_4RINGa<-subset(Dataset, Dataset$AGE==9& Dataset$Year==2016&Dataset$variable=='CSW4')
Age_4LRa["Ring"] <- "Last Ring"
Age_4LRb["Ring"] <- "Last Ring"
Age_8_4RINGa["Ring"]<-"Age Ring"
Age_8_4RINGb["Ring"]<-"Age Ring"
Age_9_4RINGa["Ring"]<-"Age Ring"
Dataset_Age4<- rbind(Age_4LRa,Age_4LRb,Age_8_4RINGa, Age_8_4RINGb, 
                     Age_9_4RINGa) 
#AGE 5 DATASET
Age_5LRa<-subset(Dataset, Dataset$AGE==5& Dataset$Year==2011&Dataset$variable=='CSW5') 
Age_5LRb<-subset(Dataset, Dataset$AGE==5& Dataset$Year==2012&Dataset$variable=='CSW5') 
Age_5LRc<-subset(Dataset, Dataset$AGE==5& Dataset$Year==2013&Dataset$variable=='CSW5') 
Age_8_5RINGa<-subset(Dataset, Dataset$AGE==8& Dataset$Year==2014&Dataset$variable=='CSW5') 
Age_8_5RINGb<-subset(Dataset, Dataset$AGE==8& Dataset$Year==2015&Dataset$variable=='CSW5')
Age_8_5RINGc<-subset(Dataset, Dataset$AGE==8& Dataset$Year==2016&Dataset$variable=='CSW5')
Age_9_5RINGa<-subset(Dataset, Dataset$AGE==9& Dataset$Year==2015&Dataset$variable=='CSW5')
Age_9_5RINGb<-subset(Dataset, Dataset$AGE==9& Dataset$Year==2016&Dataset$variable=='CSW5')
Age_10_5RINGa<-subset(Dataset, Dataset$AGE==10& Dataset$Year==2016&Dataset$variable=='CSW5')
Age_5LRa["Ring"] <- "Last Ring"
Age_5LRb["Ring"] <- "Last Ring"
Age_5LRc["Ring"] <- "Last Ring"
Age_8_5RINGa["Ring"]<-"Age Ring"
Age_8_5RINGb["Ring"]<-"Age Ring"
Age_8_5RINGc["Ring"]<-"Age Ring"
Age_9_5RINGa["Ring"]<-"Age Ring"
Age_9_5RINGb["Ring"]<-"Age Ring"
Age_10_5RINGa["Ring"]<-"Age Ring"
Dataset_Age5<- rbind(Age_5LRa,Age_5LRb,Age_5LRc,
                     Age_8_5RINGa,Age_8_5RINGb,Age_8_5RINGc,Age_9_5RINGa,Age_9_5RINGb,Age_10_5RINGa)
#AGE 6 DATASET
Age_6LRa<-subset(Dataset, Dataset$AGE==6& Dataset$Year==2011&Dataset$variable=='CSW6') 
Age_6LRb<-subset(Dataset, Dataset$AGE==6& Dataset$Year==2012&Dataset$variable=='CSW6') 
Age_6LRc<-subset(Dataset, Dataset$AGE==6& Dataset$Year==2013&Dataset$variable=='CSW6') 
Age_6LRd<-subset(Dataset, Dataset$AGE==6& Dataset$Year==2014&Dataset$variable=='CSW6') 
Age_8_6RINGa<-subset(Dataset, Dataset$AGE==8& Dataset$Year==2013&Dataset$variable=='CSW6') 
Age_8_6RINGb<-subset(Dataset, Dataset$AGE==8& Dataset$Year==2014&Dataset$variable=='CSW6') 
Age_8_6RINGc<-subset(Dataset, Dataset$AGE==8& Dataset$Year==2015&Dataset$variable=='CSW6') 
Age_8_6RINGd<-subset(Dataset, Dataset$AGE==8& Dataset$Year==2016&Dataset$variable=='CSW6') 
Age_9_6RINGa<-subset(Dataset, Dataset$AGE==9& Dataset$Year==2014&Dataset$variable=='CSW6') 
Age_9_6RINGb<-subset(Dataset, Dataset$AGE==9& Dataset$Year==2015&Dataset$variable=='CSW6') 
Age_9_6RINGc<-subset(Dataset, Dataset$AGE==9& Dataset$Year==2016&Dataset$variable=='CSW6') 
Age_10_6RINGa<-subset(Dataset, Dataset$AGE==10& Dataset$Year==2015&Dataset$variable=='CSW6') 
Age_10_6RINGb<-subset(Dataset, Dataset$AGE==10& Dataset$Year==2016&Dataset$variable=='CSW6') 
Age_11_6RINGa<-subset(Dataset, Dataset$AGE==11& Dataset$Year==2016&Dataset$variable=='CSW6') 
Age_6LRa["Ring"] <- "Last Ring"
Age_6LRb["Ring"] <- "Last Ring"
Age_6LRc["Ring"] <- "Last Ring"
Age_6LRd["Ring"] <- "Last Ring"
Age_8_6RINGa["Ring"]<-"Age Ring"
Age_8_6RINGb["Ring"]<-"Age Ring"
Age_8_6RINGc["Ring"]<-"Age Ring"
Age_8_6RINGd["Ring"]<-"Age Ring"
Age_9_6RINGa["Ring"]<-"Age Ring"
Age_9_6RINGb["Ring"]<-"Age Ring"
Age_9_6RINGc["Ring"]<-"Age Ring"
Age_10_6RINGa["Ring"]<-"Age Ring"
Age_10_6RINGb["Ring"]<-"Age Ring"
Age_11_6RINGa["Ring"]<-"Age Ring"
Dataset_Age6<- rbind(Age_6LRa,Age_6LRb,Age_6LRc,Age_6LRd,
                     Age_8_6RINGa,Age_8_6RINGb, Age_8_6RINGc,Age_8_6RINGd,Age_9_6RINGa,Age_9_6RINGb,
                     Age_9_6RINGc,Age_10_6RINGa,Age_10_6RINGb,Age_11_6RINGa)

#AGE 7 DATASET
Age_7LRa<-subset(Dataset, Dataset$AGE==7& Dataset$Year==2011&Dataset$variable=='CSW7') 
Age_7LRb<-subset(Dataset, Dataset$AGE==7& Dataset$Year==2012&Dataset$variable=='CSW7') 
Age_7LRc<-subset(Dataset, Dataset$AGE==7& Dataset$Year==2013&Dataset$variable=='CSW7') 
Age_7LRd<-subset(Dataset, Dataset$AGE==7& Dataset$Year==2014&Dataset$variable=='CSW7') 
Age_7LRe<-subset(Dataset, Dataset$AGE==7& Dataset$Year==2015&Dataset$variable=='CSW7') 
Age_8_7RINGa<-subset(Dataset, Dataset$AGE==8& Dataset$Year==2012&Dataset$variable=='CSW7') 
Age_8_7RINGb<-subset(Dataset, Dataset$AGE==8& Dataset$Year==2013&Dataset$variable=='CSW7') 
Age_8_7RINGc<-subset(Dataset, Dataset$AGE==8& Dataset$Year==2014&Dataset$variable=='CSW7') 
Age_8_7RINGd<-subset(Dataset, Dataset$AGE==8& Dataset$Year==2015&Dataset$variable=='CSW7') 
Age_8_7RINGe<-subset(Dataset, Dataset$AGE==8& Dataset$Year==2016&Dataset$variable=='CSW7') 
Age_9_7RINGa<-subset(Dataset, Dataset$AGE==9& Dataset$Year==2013&Dataset$variable=='CSW7') 
Age_9_7RINGb<-subset(Dataset, Dataset$AGE==9& Dataset$Year==2014&Dataset$variable=='CSW7') 
Age_9_7RINGc<-subset(Dataset, Dataset$AGE==9& Dataset$Year==2015&Dataset$variable=='CSW7') 
Age_9_7RINGd<-subset(Dataset, Dataset$AGE==9& Dataset$Year==2016&Dataset$variable=='CSW7') 
Age_10_7RINGa<-subset(Dataset, Dataset$AGE==10& Dataset$Year==2014&Dataset$variable=='CSW7') 
Age_10_7RINGb<-subset(Dataset, Dataset$AGE==10& Dataset$Year==2015&Dataset$variable=='CSW7') 
Age_10_7RINGc<-subset(Dataset, Dataset$AGE==10& Dataset$Year==2016&Dataset$variable=='CSW7') 
Age_11_7RINGa<-subset(Dataset, Dataset$AGE==11& Dataset$Year==2015&Dataset$variable=='CSW7') 
Age_11_7RINGb<-subset(Dataset, Dataset$AGE==11& Dataset$Year==2016&Dataset$variable=='CSW7') 
Age_7LRa["Ring"] <- "Last Ring"
Age_7LRb["Ring"] <- "Last Ring"
Age_7LRc["Ring"] <- "Last Ring"
Age_7LRd["Ring"] <- "Last Ring"
Age_7LRe["Ring"] <- "Last Ring"
Age_8_7RINGa["Ring"]<-"Age Ring"
Age_8_7RINGb["Ring"]<-"Age Ring"
Age_8_7RINGc["Ring"]<-"Age Ring"
Age_8_7RINGd["Ring"]<-"Age Ring"
Age_8_7RINGe["Ring"]<-"Age Ring"
Age_9_7RINGa["Ring"]<-"Age Ring"
Age_9_7RINGb["Ring"]<-"Age Ring"
Age_9_7RINGc["Ring"]<-"Age Ring"
Age_9_7RINGd["Ring"]<-"Age Ring"
Age_10_7RINGa["Ring"]<-"Age Ring"
Age_10_7RINGb["Ring"]<-"Age Ring"
Age_10_7RINGc["Ring"]<-"Age Ring"
Age_11_7RINGa["Ring"]<-"Age Ring"
Age_11_7RINGb["Ring"]<-"Age Ring"
Dataset_Age7<- rbind(Age_7LRa,Age_7LRb,Age_7LRc,Age_7LRd,Age_7LRe,Age_8_7RINGa,Age_8_7RINGb,Age_8_7RINGc,Age_8_7RINGd,
                     Age_8_7RINGe,Age_9_7RINGa,Age_9_7RINGb,Age_9_7RINGc,Age_9_7RINGd,Age_10_7RINGa,Age_10_7RINGb,Age_10_7RINGc,
                     Age_11_7RINGa,Age_11_7RINGb)


png(file='Histograms.png', res=200, width=13, height=10, units ="in") #default adjust=1 for density curve
x3<-ggplot(Dataset_Age3, aes(x=value, color=Ring, fill=Ring)) +
  geom_histogram(position="identity", alpha=0.2)+
  geom_density(alpha=0.5, adjust=1)+
  scale_color_manual(values=c("#999999", "#E69F00"))+
  scale_fill_manual(values=c("#999999", "#E69F00" ))+
  labs(title="Age 3",x="Scale increment length (mm)", y="Count")+
  theme_classic()
x3<-x3+ggtitle("Age 3")+theme_set(theme_bw(base_size=14,base_family=
                                             'Times New Roman')+
                                    theme(panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank()))
x3<-x3 + theme(legend.position="none")+scale_x_continuous(breaks=seq(0,max(Dataset_Age3$value+0.1, na.rm=TRUE),0.1))


x4<-ggplot(Dataset_Age4, aes(x=value, color=Ring, fill=Ring)) +
  geom_histogram(position="identity", alpha=0.2)+
  geom_density(alpha=0.5,adjust=1)+
  scale_color_manual(values=c("#999999", "#E69F00"))+
  scale_fill_manual(values=c("#999999", "#E69F00" ))+
  labs(title="Age 4",x="Scale increment length (mm)", y="Count")+
  theme_classic()
x4<-x4+ggtitle("Age 4")+theme_set(theme_bw(base_size=14,base_family=
                                             'Times New Roman')+
                                    theme(panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank()))
x4<-x4 + theme(legend.position="none")+scale_x_continuous(breaks=seq(0,max(Dataset_Age4$value+0.1, na.rm=TRUE),0.1)) 

x5<-ggplot(Dataset_Age5, aes(x=value, color=Ring, fill=Ring)) +
  geom_histogram(position="identity", alpha=0.2)+
  geom_density(alpha=0.5,adjust=1)+
  scale_color_manual(values=c("#999999", "#E69F00"))+
  scale_fill_manual(values=c("#999999", "#E69F00" ))+
  labs(title="Age 5",x="Scale increment length (mm)", y="Count")+
  theme_classic()
x5<-x5+ggtitle("Age 5")+theme_set(theme_bw(base_size=14,base_family=
                                             'Times New Roman')+
                                    theme(panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank()))
x5<-x5 + theme(legend.position="none")+scale_x_continuous(breaks=seq(0,max(Dataset_Age5$value+0.1, na.rm=TRUE),0.1)) 

x6<-ggplot(Dataset_Age6, aes(x=value, color=Ring, fill=Ring)) +
  geom_histogram(position="identity", alpha=0.2)+
  geom_density(alpha=0.5,adjust=1)+
  scale_color_manual(values=c("#999999", "#E69F00"))+
  scale_fill_manual(values=c("#999999", "#E69F00" ))+
  labs(title="Age 6",x="Scale increment length (mm)", y="Count")+
  theme_classic()
x6<-x6+ggtitle("Age 6")+theme_set(theme_bw(base_size=14,base_family=
                                             'Times New Roman')+
                                    theme(panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank()))
x6<-x6 + theme(legend.position="none")+scale_x_continuous(breaks=seq(0,max(Dataset_Age6$value+0.1, na.rm=TRUE),0.1)) 

x7<-ggplot(Dataset_Age7, aes(x=value, color=Ring, fill=Ring)) +
  geom_histogram(position="identity", alpha=0.2)+
  geom_density(alpha=0.5,adjust=1)+
  scale_color_manual(values=c("#999999", "#E69F00"))+
  scale_fill_manual(values=c("#999999", "#E69F00" ))+
  labs(title="Age 7",x="Scale increment length (mm)", y="Count")+
  theme_classic()
x7<-x7+ggtitle("Age 7")+theme_set(theme_bw(base_size=14,base_family=
                                             'Times New Roman')+
                                    theme(panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank()))
x7<-x7 + theme(legend.position="none")+scale_x_continuous(breaks=seq(0,max(Dataset_Age7$value+0.1, na.rm=TRUE),0.1))
x7<-ggplot(Dataset_Age7, aes(value, fill = Ring)) + geom_histogram(alpha=0.5, position = 'identity')+
  ylab("Frequency")+ xlab("Scale increment length (mm)")+geom_line(stat="density")+
  scale_x_continuous(breaks=seq(0,max(Age_7$value, na.rm=TRUE),0.5) )                       
x7<-x7 + theme(legend.title=element_blank())+ggtitle("Age 7")
library(grid)
g6<-grid.arrange(x3,x4,x5,x6,x7, ncol=2)
dev.off()
#######################################################################################################
png(file='Density.png', res=200, width=13, height=10, units ="in") #default adjust=1 for density curve
x3<-ggplot(Dataset_Age3, aes(x=value, color=Ring, fill=Ring)) +
  geom_density(alpha=0.5, adjust=1)+
  scale_color_manual(values=c("#999999", "#E69F00"))+
  scale_fill_manual(values=c("#999999", "#E69F00" ))+
  labs(title="Age 3",x="Scale increment length (mm)",y="Density")+
  theme_classic()
x3<-x3+ggtitle("Age 3")+theme_set(theme_bw(base_size=14,base_family=
                                             'Times New Roman')+
                                    theme(panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank()))
x3<-x3 + theme(legend.position="none")+scale_x_continuous(breaks=seq(0,max(Dataset_Age3$value+0.1, na.rm=TRUE),0.1))


x4<-ggplot(Dataset_Age4, aes(x=value, color=Ring, fill=Ring)) +
  geom_density(alpha=0.5,adjust=1)+
  scale_color_manual(values=c("#999999", "#E69F00"))+
  scale_fill_manual(values=c("#999999", "#E69F00" ))+
  labs(title="Age 4",x="Scale increment length (mm)", y="Density")+
  theme_classic()
x4<-x4+ggtitle("Age 4")+theme_set(theme_bw(base_size=14,base_family=
                                             'Times New Roman')+
                                    theme(panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank()))
x4<-x4 + theme(legend.position="none")+scale_x_continuous(breaks=seq(0,max(Dataset_Age4$value+0.1, na.rm=TRUE),0.1)) 

x5<-ggplot(Dataset_Age5, aes(x=value, color=Ring, fill=Ring)) +
  geom_density(alpha=0.5,adjust=1)+
  scale_color_manual(values=c("#999999", "#E69F00"))+
  scale_fill_manual(values=c("#999999", "#E69F00" ))+
  labs(title="Age 5",x="Scale increment length (mm)", y="Density")+
  theme_classic()
x5<-x5+ggtitle("Age 5")+theme_set(theme_bw(base_size=14,base_family=
                                             'Times New Roman')+
                                    theme(panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank()))
x5<-x5 + theme(legend.position="none")+scale_x_continuous(breaks=seq(0,max(Dataset_Age5$value+0.1, na.rm=TRUE),0.1)) 

x6<-ggplot(Dataset_Age6, aes(x=value, color=Ring, fill=Ring)) +
  geom_density(alpha=0.5,adjust=1)+
  scale_color_manual(values=c("#999999", "#E69F00"))+
  scale_fill_manual(values=c("#999999", "#E69F00" ))+
  labs(title="Age 6",x="Scale increment length (mm)", y="Density")+
  theme_classic()
x6<-x6+ggtitle("Age 6")+theme_set(theme_bw(base_size=14,base_family=
                                             'Times New Roman')+
                                    theme(panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank()))
x6<-x6 + theme(legend.position="none")+scale_x_continuous(breaks=seq(0,max(Dataset_Age6$value+0.1, na.rm=TRUE),0.1)) 

x7<-ggplot(Dataset_Age7, aes(x=value, color=Ring, fill=Ring)) +
  geom_density(alpha=0.5,adjust=1)+
  scale_color_manual(values=c("#999999", "#E69F00"))+
  scale_fill_manual(values=c("#999999", "#E69F00" ))+
  labs(title="Age 7",x="Scale increment length (mm)", y="Density")+
  theme_classic()
x7<-x7+ggtitle("Age 7")+theme_set(theme_bw(base_size=14,base_family=
                                             'Times New Roman')+
                                    theme(panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank()))
x7<-x7 + theme(legend.position="none")+scale_x_continuous(breaks=seq(0,max(Dataset_Age7$value+0.1, na.rm=TRUE),0.1))


g6<-grid.arrange(x3,x4,x5,x6,x7, ncol=2)
dev.off()
####################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
###################################################################################
#Fit gaussian mixture models to age ring data and plot (AGE 3)
#mixtools package
png(file='Mixture_Models.png', res=200, width=13, height=10, units ="in") 
par(mfrow=c(5,2)) 
Data_Age3<- subset(Dataset_Age3, Dataset_Age3$Ring=='Age Ring') 
Data_Age3<-Data_Age3[,c(6)]
mixmdl3 <- normalmixEM(Data_Age3)
x3<-plot(mixmdl3, which=2, cex.axis=1, cex.lab=1, cex.main=1.1,
         main2="Age 3 (Inner Rings;n=6)", xlab2="Scale increments (mm)", xlim=c(0.2,1.5))
x3<-lines(density(Data_Age3), lty=2, lwd=2)
x3<-legend("topright", c("mature", "immature"), col = c(2,3), lty = c(1,1), lwd=2,
           merge = TRUE, bty = "n")
summary(mixmdl3)
## Fit normal distr to Last Ring data (AGE 3)
#fitdistrplus package
Data_Age3<- subset(Dataset_Age3, Dataset_Age3$Ring=='Last Ring') 
Data_Age3<-as.data.frame(Data_Age3[,c(6)])
f <- apply(Data_Age3, 2, fitdist, "norm")
colnames(Data_Age3) <- "Measurement"
Data_Age3$Measurement<-as.numeric(Data_Age3$Measurement)
hx <- dnorm(Data_Age3$Measurement)
Measurement<-Data_Age3$Measurement
x<-fitdist(Measurement, "norm")
x3a<-denscomp(x, main="Age 3 (Outer Ring)", xlab="Scale increments (mm)", xlim = c(0.2,1.5),addlegend=T,
              legendtext=c("mature (n=37)"),lwd=2)
#Fit gaussian mixture models to age ring data and plot (AGE 4)
Data_Age4<- subset(Dataset_Age4, Dataset_Age4$Ring=='Age Ring') 
Data_Age4<-Data_Age4[,c(6)]
mixmdl4 <- normalmixEM(Data_Age4)
x4<-plot(mixmdl4, which=2, cex.axis=1, cex.lab=1, cex.main=1.1,xlim = c(0.2,1.0),
         main2="Age 4 (Inner Rings;n=50)", xlab2="Scale increments (mm)")
x4<-lines(density(Data_Age4), lty=2, lwd=2)
x4<-legend("topright", c("mature", "immature"), col = c(2,3), lty = c(1,1), lwd=2,
           merge = TRUE, bty = "n")
summary(mixmdl4)
#Age 4 did not converge
## Fit normal distr to Last Ring data (AGE 4)
Data_Age4<- subset(Dataset_Age4, Dataset_Age4$Ring=='Last Ring') 
Data_Age4<-as.data.frame(Data_Age4[,c(6)])
f <- apply(Data_Age4, 2, fitdist, "norm")
colnames(Data_Age4) <- "Measurement"
Data_Age4$Measurement<-as.numeric(Data_Age4$Measurement)
hx <- dnorm(Data_Age4$Measurement)
Measurement<-Data_Age4$Measurement
x<-fitdist(Measurement, "norm")
x4a<-denscomp(x, main="Age 4 (Outer Ring)", xlab="Scale increments (mm)", xlim = c(0.2,1.0),addlegend=T,
              legendtext=c("mature (n=230)"), lwd=2)
#Fit gaussian mixture models to age ring data and plot (AGE 5)
Data_Age5<- subset(Dataset_Age5, Dataset_Age5$Ring=='Age Ring') 
Data_Age5<-Data_Age5[,c(6)]
mixmdl5 <- normalmixEM(Data_Age5)
x5<-plot(mixmdl5, which=2, cex.axis=1, cex.lab=1, cex.main=1.1,xlim = c(0.1,0.8),
         main2="Age 5 (Inner Rings; n=159)", xlab2="Scale increments (mm)")
x5<-legend("topright", c("mature", "immature"), col = c(2,3), lty = c(1,1), lwd=2,
           merge = TRUE, bty = "n")
x5<-lines(density(Data_Age5), lty=2, lwd=2)
summary(mixmdl5)
## Fit normal distr to Last Ring data (AGE 5)
Data_Age5<- subset(Dataset_Age5, Dataset_Age5$Ring=='Last Ring') 
Data_Age5<-as.data.frame(Data_Age5[,c(6)])
f <- apply(Data_Age5, 2, fitdist, "norm")
colnames(Data_Age5) <- "Measurement"
Data_Age5$Measurement<-as.numeric(Data_Age5$Measurement)
hx <- dnorm(Data_Age5$Measurement)
Measurement<-Data_Age5$Measurement
x<-fitdist(Measurement, "norm")
x5a<-denscomp(x, main="Age 5 (Outer Ring)", xlab="Scale increments (mm)", xlim = c(0.1,0.8),addlegend=T,
              legendtext=c("mature (n=252)"),lwd=2)
#Fit gaussian mixture models to age ring data and plot (AGE 6)
Data_Age6<- subset(Dataset_Age6, Dataset_Age6$Ring=='Age Ring') 
Data_Age6<-Data_Age6[,c(6)]
mixmdl6 <- normalmixEM(Data_Age6)
x6<-plot(mixmdl6, which=2, cex.axis=1, cex.lab=1, cex.main=1.1,xlim = c(0,0.7),
         main2="Age 6 (Inner Rings; n=357)", xlab2="Scale increments (mm)")
x6<-lines(density(Data_Age6), lty=2, lwd=2)
x6<-legend("topright", c("mature", "immature"), col = c(2,3), lty = c(1,1), lwd=2,
           merge = TRUE, bty = "n")
summary(mixmdl6)
## Fit normal distr to Last Ring data (AGE 6)
Data_Age6<- subset(Dataset_Age6, Dataset_Age6$Ring=='Last Ring') 
Data_Age6<-as.data.frame(Data_Age6[,c(6)])
f <- apply(Data_Age6, 2, fitdist, "norm")
colnames(Data_Age6) <- "Measurement"
Data_Age6$Measurement<-as.numeric(Data_Age6$Measurement)
hx <- dnorm(Data_Age6$Measurement)
Measurement<-Data_Age6$Measurement
x<-fitdist(Measurement, "norm")
x6a<-denscomp(x, main="Age 6 (Outer Ring)", xlab="Scale increments (mm)", addlegend=T,
              legendtext=c("mature (n=375)"),lwd=2,xlim = c(0,0.7))

#Fit gaussian mixture models to age ring data and plot (AGE 7)
Data_Age7<- subset(Dataset_Age7, Dataset_Age7$Ring=='Age Ring') 
Data_Age7<-Data_Age7[,c(6)]
mixmdl7 <- normalmixEM(Data_Age7)
x7<-plot(mixmdl7, which=2, cex.axis=1, cex.lab=1, cex.main=1.1,xlim = c(0,0.7),
         main2="Age 7 (Inner Rings; n=477)", xlab2="Scale increments (mm)")
x7<-lines(density(Data_Age7), lty=2, lwd=2)
x7<-legend("topright", c("mature", "immature"), col = c(2,3), lty = c(1,1), lwd=2,
           merge = TRUE, bty = "n")
summary(mixmdl7)
## Fit normal distr to Last Ring data (AGE 6)
Data_Age7<- subset(Dataset_Age7, Dataset_Age7$Ring=='Last Ring') 
Data_Age7<-as.data.frame(Data_Age7[,c(6)])
f <- apply(Data_Age7, 2, fitdist, "norm")
colnames(Data_Age7) <- "Measurement"
Data_Age7$Measurement<-as.numeric(Data_Age7$Measurement)
hx <- dnorm(Data_Age7$Measurement)
Measurement<-Data_Age7$Measurement
x<-fitdist(Measurement, "norm")
x7a<-denscomp(x, main="Age 7 (Outer Ring)", xlab="Scale increments (mm)", addlegend=T,
              legendtext=c("mature (n=377)"),lwd=2,xlim = c(0,0.7))
dev.off()
citation("mixtools")


#ANOVA tests with different cohorts and years
Data <- read.csv("scale_increments.csv", sep=",", header = TRUE, check.names = FALSE)  #Read in csv file, subset the spawning stock
Dataset <- melt(Data, id=c("Year", "AGE", "SEX","LENGTH"), na.rm=TRUE)
#create dataset for 2011 with the last ring
#AGE 3 DATASET
Age3LR<-subset(Dataset, Dataset$AGE==3 & Dataset$variable=='CSW3') 
Age4LR<-subset(Dataset, Dataset$AGE==4 & Dataset$variable=='CSW4') 
Age5LR<-subset(Dataset, Dataset$AGE==5 & Dataset$variable=='CSW5') 
Age6LR<-subset(Dataset, Dataset$AGE==6 & Dataset$variable=='CSW6') 
Age7LR<-subset(Dataset, Dataset$AGE==7 & Dataset$variable=='CSW7') 
Age8LR<-subset(Dataset, Dataset$AGE==8 & Dataset$variable=='CSW8') 
Age3LR$Year <- factor(Age3LR$Year)
Age4LR$Year <- factor(Age4LR$Year)
Age5LR$Year <- factor(Age5LR$Year)
Age6LR$Year <- factor(Age6LR$Year)
Age7LR$Year <- factor(Age7LR$Year)
Age8LR$Year <- factor(Age8LR$Year)
fit1 <- aov(value~Year, data=Age3LR)
TukeyHSD(fit1)
summaryBy(value ~ Year, data = Age3LR,
          FUN = function(x) { c(mean = mean(x), sd = sd(x), n=length(x), se=sd(x)/sqrt(length(x))) } )

fit2 <- aov(value~Year, data=Age4LR)
TukeyHSD(fit2)
summaryBy(value ~ Year, data = Age4LR,
          FUN = function(x) { c(mean = mean(x), sd = sd(x), n=length(x), se=sd(x)/sqrt(length(x))) })

fit3 <- aov(value~Year, data=Age5LR)
TukeyHSD(fit3)
summaryBy(value ~ Year, data = Age5LR,
          FUN = function(x) { c(mean = mean(x), sd = sd(x), n=length(x), se=sd(x)/sqrt(length(x))) })

fit4 <- aov(value~Year, data=Age6LR)
TukeyHSD(fit4)
summaryBy(value ~ Year, data = Age6LR,
          FUN = function(x) { c(mean = mean(x), sd = sd(x), n=length(x), se=sd(x)/sqrt(length(x))) })

fit5 <- aov(value~Year, data=Age7LR)
TukeyHSD(fit5)
summaryBy(value ~ Year, data = Age7LR,
          FUN = function(x) { c(mean = mean(x), sd = sd(x), n=length(x), se=sd(x)/sqrt(length(x))) })

fit6 <- aov(value~Year, data=Age8LR)
TukeyHSD(fit6)
summaryBy(value ~ Year, data = Age8LR,
          FUN = function(x) { c(mean = mean(x), sd = sd(x), n=length(x), se=sd(x)/sqrt(length(x))) })
C<-subset(FIGDATA, select=c(Year, tot_obs_egg,tot_est_egg)) 

Dataset <- read.csv("scale_incrementsCI.csv", sep=",", header = TRUE, check.names = FALSE)  #Read in csv file, subset the spawning stock
Dataset$Year <- factor(Dataset$Year)
Dataset$Age <- factor(Dataset$Age)
prepanel.ci <- function(x, y, ly, uy, subscripts, ...)
{
  x <- as.numeric(x)
  ly <- as.numeric(ly[subscripts])
  uy <- as.numeric(uy[subscripts])
  list(ylim = range(y, uy, ly, finite = TRUE)
}

panel.ci <- function(x, y, ly, uy, subscripts, pch = c(15,22), col.line =
                       c(1,2), ...)
{
  x <- as.numeric(x)
  y <- as.numeric(y)
  ly <- as.numeric(ly[subscripts])
  uy <- as.numeric(uy[subscripts])
  panel.arrows(x, ly, x, uy, col=col.line,
               length = 0.25, unit = "native",
               angle = 90, code = 3)
  panel.xyplot(x, y, pch = pch, col.line=col.line,...)
  panel.abline(h=mean(y), lwd=2, lty=1, col="lightgreen" )
}

library(lattice)  
tiff(file="Mean_CI.tiff", width=7,height=4,units="in", res=600)                             
xyplot(Mean~ Year|Age,  data=Dataset,                                    
       groups=Age, col=1, pch=c(16), ylim=c(0,1.2),font.lab=3, font.axis=3,                                                    
       ylab=list("Mean increment length (mm)",font=6, cex=1.2), xlab=list("Year",font=6, cex=1.2),
       cex=0.9, layout=c(3,2), scales=list(x=list(rot=90)),             
       ly = Dataset$Lower_CI,                                                                          
       uy = Dataset$Upper_CI,                                                                          
       panel = panel.superpose,                                                                 
       panel.groups = panel.ci,                                                                 
       type="p",  
       
       par.settings = list(superpose.symbol = list(pch=c(16),
                                                   col=1, cex=1)))
update(trellis.last.object(),strip=strip.custom(strip.levels=TRUE,bg=8))                         
update(trellis.last.object(),scales=list(cex=0.75, font=6))
dev.off()





########################################################################################################
help(hist)
boxplot(Age_3$value~Age_3$variable, Data=Age_3, ylab="Scale increment length (mm)", xlab="Scale increment")
boxplot(Age_4$value~Age_4$variable, Data=Age_4, ylab="Scale increment length (mm)", xlab="Scale increment")

#Power Test values (pg. 208 Cohen 1977; Chapter 6.5; example 6.8])
cohen.ES(test = c("p", "t", "r", "anov", "chisq", "f2"),
         size = c("small", "medium", "large"))
cohen.ES(test="p", size="small")
pwr.p.test(h=0.321,n=100,sig.level=0.05,alternative="two.sided")
pwr.p.test(h=0.2,power=0.99,sig.level=0.05,alternative="two.sided")
pwr.p.test(h=0.20,power=0.95,sig.level=0.05,alternative="two.sided")
pwr.p.test(h=0.2,power=0.90,sig.level=0.05,alternative="two.sided")
pwr.p.test(h=0.2,power=0.85,sig.level=0.05,alternative="two.sided")
pwr.p.test(h=0.2,power=0.80,sig.level=0.05,alternative="two.sided")
pwr.p.test(h=0.2,power=0.75,sig.level=0.05,alternative="two.sided")
pwr.p.test(h=0.2,power=0.70,sig.level=0.05,alternative="two.sided")
pwr.p.test(h=0.2,power=0.65,sig.level=0.05,alternative="two.sided")
pwr.p.test(h=0.2,power=0.60,sig.level=0.05,alternative="two.sided")
pwr.p.test(h=0.2,power=0.55,sig.level=0.05,alternative="two.sided")
pwr.p.test(h=0.2,power=0.50,sig.level=0.05,alternative="two.sided")

pwr.p.test(h=0.3,power=0.99,sig.level=0.05,alternative="two.sided")
pwr.p.test(h=0.3,power=0.95,sig.level=0.05,alternative="two.sided")
pwr.p.test(h=0.3,power=0.90,sig.level=0.05,alternative="two.sided")
pwr.p.test(h=0.3,power=0.85,sig.level=0.05,alternative="two.sided")
pwr.p.test(h=0.3,power=0.80,sig.level=0.05,alternative="two.sided")
pwr.p.test(h=0.3,power=0.75,sig.level=0.05,alternative="two.sided")
pwr.p.test(h=0.3,power=0.70,sig.level=0.05,alternative="two.sided")
pwr.p.test(h=0.3,power=0.65,sig.level=0.05,alternative="two.sided")
pwr.p.test(h=0.3,power=0.60,sig.level=0.05,alternative="two.sided")
pwr.p.test(h=0.3,power=0.55,sig.level=0.05,alternative="two.sided")
pwr.p.test(h=0.3,power=0.50,sig.level=0.05,alternative="two.sided")

pwr.p.test(h=0.142,n=100,sig.level=0.20,alternative="two.sided")
p.p.two <- pwr.p.test(h=0.3, power = 0.99, alternative = "two.sided")
x<-plot(p.p.two, xlab="Sample Size", xlim=c(0,500), ylim=c(0,100) )
citation(package = "pwr")
readCitationFile(file, meta = NULL)