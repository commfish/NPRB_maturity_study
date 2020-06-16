#CSW1 is from focus to annulus1, CSW2 is from annulus1 to annulus 2 , etc.  
#The last scale measurement for each specimen is the plus growth (last annulus to edge of scale).
#The number of increment measurements for herring is equal to the number of annuli and the plus growth, which adds up to the age of the fish in this dataset

#libraries---
library(plyr)
library(extrafont)
library(ggplot2)
library(readr)
library(reshape2)
library(grid)
library(dplyr)
library(tidyr)
source("code/functions.R")
loadfonts(device="win")
windowsFonts(Times=windowsFont("TT Times New Roman"))
theme_set(theme_bw(base_size=12, base_family='Times New Roman')+ 
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank())+
            theme(axis.text.x = element_text(angle = 90, hjust = 1)))

#SITKA SCALE MEASUREMENTS (DETLEF-MTA Lab)
# data ----
data <- read_csv("data/scale_measurements.csv") 
names(data) <- c('year', 'system_code', 'species_code', 'location',
                 'image_name', 'age', 'sex', 'length','CSWa','CSWb','CSWc',
                 'CSWd','CSWe','CSWf','CSWg','CSWh','CSWi','CSWj',
                 'CSWk','CSWl','CSWm','CSWn')
data %>% 
  mutate(year=factor(year),
         age=factor(age),
         CSWa=as.numeric(CSWa),
         CSWb=as.numeric(CSWb),
         CSWc=as.numeric(CSWc),
         CSWd=as.numeric(CSWd),
         CSWe=as.numeric(CSWe),
         CSWf=as.numeric(CSWf),
         CSWg=as.numeric(CSWg),
         CSWh=as.numeric(CSWh),
         CSWi=as.numeric(CSWi),
         CSWj=as.numeric(CSWj),
         CSWk=as.numeric(CSWk),
         CSWl=as.numeric(CSWl),
         CSWm=as.numeric(CSWm),
         CSWn=as.numeric(CSWn))->data

#total is the sum of the total scale length from focus to outer edge
data %>%
  rowwise() %>% 
  mutate(total=sum(CSWa,CSWb,CSWc,CSWd,CSWe,CSWf,CSWg,CSWh,CSWi,CSWj, 
                  CSWk, CSWl,CSWm, CSWn, na.rm=TRUE))->data
#proportions from focus to each increment/total increment length
data %>%
  mutate(prop1=CSWa/total,#focus to annulus one
         prop2=(CSWa+CSWb)/total, #focus to annulus two
         prop3=(CSWa+CSWb+CSWc)/total,
         prop4=(CSWa+CSWb+CSWc+CSWd)/total,
         prop5=(CSWa+CSWb+CSWc+CSWd+CSWe)/total,
         prop6=(CSWa+CSWb+CSWc+CSWd+CSWe+CSWf)/total,
         prop7=(CSWa+CSWb+CSWc+CSWd+CSWe+CSWf+CSWg)/total,
         prop8=(CSWa+CSWb+CSWc+CSWd+CSWe+CSWf+CSWg+CSWh)/total,
         prop9=(CSWa+CSWb+CSWc+CSWd+CSWe+CSWf+CSWg+CSWh+CSWi)/total,
         prop10=(CSWa+CSWb+CSWc+CSWd+CSWe+CSWf+CSWg+CSWh+CSWi+CSWj)/total,
         prop11=(CSWa+CSWb+CSWc+CSWd+CSWe+CSWf+CSWg+CSWh+CSWi+CSWj+CSWk)/total,
         prop12=(CSWa+CSWb+CSWc+CSWd+CSWe+CSWf+CSWg+CSWh+CSWi+CSWj+CSWk+CSWl)/total,
         prop13=(CSWa+CSWb+CSWc+CSWd+CSWe+CSWf+CSWg+CSWh+CSWi+CSWj+CSWk+CSWl+CSWm)/total,
         prop14=(CSWa+CSWb+CSWc+CSWd+CSWe+CSWf+CSWg+CSWh+CSWi+CSWj+CSWk+CSWl+CSWm+CSWn)/total)-> dataset 

# analysis ---- 
summary_table<- table(dataset$year,dataset$age) #calculate sample sizes
three<-dataset[dataset$age=="3",] #create dataset for age-3
four<-dataset[dataset$age=="4",] #create dataset for age-4
five<-dataset[dataset$age=="5",]  #create dataset for age-5
six<-dataset[dataset$age=="6",]  #create dataset for age-6

#age-3 dataset (complete cases only)
drop <- c("system_code", "species_code", "location",
          "image_name", "length",
          "CSWd","CSWe","CSWf","CSWg","CSWh","CSWi","CSWj",
          "CSWk","CSWl","CSWm","CSWn","prop3", "prop4", "prop5","prop6",
          "prop7", "prop8","prop9", "prop10", "prop11", "prop12","prop13","prop14")
three <- three[,!(names(three) %in% drop)] #drop variables
complete.cases(three) #only keep complete cases (delete cases where one scale measurement is missing)
three<-three[complete.cases(three),]
three["CSWd"] <-NA 
three["prop3"] <-NA 

#age-4 dataset (compelte cases only)
four<-dataset[dataset$age=="4",]  
drop <- c("system_code", "species_code", "location",
          "image_name", "length",
          "CSWe","CSWf","CSWg","CSWh","CSWi","CSWj",
          "CSWk","CSWl","CSWm","CSWn","prop4", "prop5","prop6",
          "prop7", "prop8","prop9", "prop10", "prop11", "prop12","prop13","prop14")
four <- four[,!(names(four) %in% drop)]
complete.cases(four)
four<-four[complete.cases(four),]

#age-5 dataset (complete cases only)
five<-dataset[dataset$age=="5",]  
drop <- c("system_code", "species_code", "location",
          "image_name", "length",
          "CSWf","CSWg","CSWh","CSWi","CSWj",
          "CSWk","CSWl","CSWm","CSWn","prop5","prop6",
          "prop7", "prop8","prop9", "prop10", "prop11", "prop12","prop13","prop14")
five <- five[,!(names(five) %in% drop)]
complete.cases(five)
five<-five[complete.cases(five),]
five["CSWf"] <-NA 
five["prop5"] <-NA 

#age-6 dataset (complete cases only)
six<-dataset[dataset$age=="6",]  
drop <- c("system_code", "species_code", "location",
          "image_name", "length",
          "CSWg","CSWh","CSWi","CSWj",
          "CSWk","CSWl","CSWm","CSWn","prop6",
          "prop7", "prop8","prop9", "prop10", "prop11", "prop12","prop13","prop14")
six <- six[,!(names(six) %in% drop)]
complete.cases(six)
six<-six[complete.cases(six),]

#merge age-3 and age-4 datasets
match.1 <- merge(three, four, by=c("year","age", "sex", "CSWa", "CSWb","CSWc","CSWd",
                                   "total", "prop1", "prop2", "prop3"),all=TRUE) 

#create NA variables so can match with age-5 and age-6 datasets
match.1["CSWe"] <-NA 
match.1["CSWf"] <-NA 
match.1["prop4"] <-NA 
match.1["prop5"] <-NA 

match.2 <- merge(five, six, by=c("year","age", "sex", "CSWa", "CSWb","CSWc","CSWd","CSWe","CSWf",
                                 "total", "prop1", "prop2", "prop3", "prop4", "prop5"),all=TRUE) 
match.3 <- merge(match.1, match.2, by=c("year","age", "sex", "CSWa", "CSWb","CSWc","CSWd","CSWe", "CSWf",
                                 "total", "prop1", "prop2", "prop3", "prop4", "prop5"),all=TRUE) 

summary_table<- table(match.3$year,match.3$age) #calculate sample sizes

#reconfigure data so that focus to last annulus correponds to correct prop. for each age
match.3["last.annulus"] <- ifelse (match.3$age=="3",match.3$prop2,
                                ifelse (match.3$age=="4",match.3$prop3,												
                                        ifelse (match.3$age=="5",match.3$prop4,match.3$prop5)))	                         

drop <- c("CSWa","CSWb","CSWc","CSWd",
          "CSWe","CSWf","total","prop3","prop4",
          "prop5")
match.3 <- match.3[,!(names(match.3) %in% drop)]

match.3 <- match.3[order(match.3$year, match.3$age, match.3$prop1),] #order dataset

#summarize dataset by props
match.3 %>% group_by(age, year) %>% 
  summarise(n = length(prop1), 
            prop1.mean = mean(prop1), prop1.var = var(prop1), prop1.se = sqrt(prop1.var), 
            prop1.max = max(prop1),
            prop1.min = min(prop1),
            prop1.mean.asintrans = asintrans(mean(prop1)),
            prop1.max.asintrans = asintrans(max(prop1)),
            prop1.min.asintrans  = asintrans(min(prop1)),
            prop2.mean = mean(prop2), prop2.var = var(prop2), prop2.se = sqrt(prop2.var),
            prop2.max = max(prop2),
            prop2.min = min(prop2),
            prop2.mean.asintrans = asintrans(mean(prop2)),
            prop2.max.asintrans = asintrans(max(prop2)),
            prop2.min.asintrans  = asintrans(min(prop2)),
            last.annulus.mean = mean(last.annulus), last.annulus.var = var(last.annulus), last.annulus.se = sqrt(last.annulus.var),
            last.annulus.max = max(last.annulus),
            last.annulus.min = min(last.annulus),
            last.annulus.mean.asintrans = asintrans(mean(last.annulus)),
            last.annulus.max.asintrans = asintrans(max(last.annulus)),
            last.annulus.min.asintrans  = asintrans(min(last.annulus)))->dataset

#calculate f
dataset %>% 
  mutate(prop1.h1 = prop1.mean.asintrans-prop1.min.asintrans,
            prop1.h2 = prop1.max.asintrans-prop1.mean.asintrans,
            prop2.h1 = prop2.mean.asintrans-prop2.min.asintrans,
            prop2.h2 = prop2.max.asintrans-prop2.mean.asintrans,
            last.annulus.h1 = last.annulus.mean.asintrans-last.annulus.min.asintrans,
            last.annulus.h2 = last.annulus.max.asintrans-last.annulus.mean.asintrans)->dataset
dataset<-as.data.frame(dataset)
write.csv(dataset, "output/summary.csv")
sample.size<-dataset[dataset$n>=100,] #dataset where sample size is >=100
write.csv(sample.size, "output/summary_samplesize.csv")

#figure of focus to first annulus/total proportion
mean_plot1<-ggplot(data=dataset,aes(x=age, y=prop1.mean, fill=age))+
  facet_wrap(~year,ncol=3,as.table=TRUE)+
  geom_bar(stat="identity")+
  geom_errorbar(aes(ymin=prop1.mean-prop1.se, ymax=prop1.mean+prop1.se, width=.2))+
  theme_bw()+
  xlab ("Age")+ ylab("Focus to first annulus proportion")+
  theme(text=element_text(family="Times New Roman", face="bold", size=12))

png(file='figs/mean_prop1.png', res=200, width=9, height=6, units ="in")
pushViewport(viewport(layout=grid.layout(1,1)))
vplayout<-function(x,y) viewport (layout.pos.row=x, layout.pos.col=y)
print(mean_plot1,vp=vplayout(1,1:1))
dev.off()

#figure of focus to last annulus/total proportion
mean_plot2<-ggplot(data=dataset,aes(x=age, y=last.mean, fill=age))+
  facet_wrap(~year,ncol=3,as.table=TRUE)+
  geom_bar(stat="identity")+
  geom_errorbar(aes(ymin=last.mean-last.se, ymax=last.mean+last.se, width=.2))+
  theme_bw()+
  xlab ("Age")+ ylab("Focus to last annulus proportion")+
  theme(text=element_text(family="Times New Roman", face="bold", size=12))
png(file='figs/mean_prop2.png', res=200, width=9, height=6, units ="in")
pushViewport(viewport(layout=grid.layout(1,1)))
vplayout<-function(x,y) viewport (layout.pos.row=x, layout.pos.col=y)
print(mean_plot2,vp=vplayout(1,1:1))
dev.off()


##MADE UP DATA
# data ----
data <- read_csv("data/madeup_data.csv") 
names(data) <- c('year', 'region', 'CSW1', 'CSW2','CSW3', 'CSW4', 'CSW5', 'CSW6')
data %>% 
  mutate(year=factor(year),
         age=factor(age),
         CSW1=as.numeric(CSW1),
         CSW2=as.numeric(CSW2),
         CSW3=as.numeric(CSW3),
         CSW4=as.numeric(CSW4),
         CSW5=as.numeric(CSW5),
         CSW6=as.numeric(CSW6), 
         total=CSW1,CSW2+CSW3+CSW4+CSW5+CSW6,
         prop=CSW1/total)-> dataset 
drop <- c("CSW1","CSW2", "CSW3", "CSW4", "CSW5", "CSW6", "total")
df <- dataset[,!(names(dataset) %in% drop)] #drop the above variables from the dataset
#prop is a variable for the increment length from the focus to the first annulus
#divided by the total scale length
#subset data by fish number
one<-df[df$fish=="1",]      
two<-df[df$fish=="2",] 
three<-df[df$fish=="3",]      
four<-df[df$fish=="4",]  
five<-df[df$fish=="5",]      
six<-df[df$fish=="6",] 
seven<-df[df$fish=="7",]      
eight<-df[df$fish=="8",]  
nine<-df[df$fish=="9",]      
ten<-df[df$fish=="10",] 
eleven<-df[df$fish=="11",]      
twelve<-df[df$fish=="12",]  

# analysis ---- 
#----ANALYSIS BY FISH
#create a figure of mean proportions by fish and region
melted <- melt(df, id.vars=c("fish", "region"))
a<-ddply(melted, c("fish", "region", "variable"), summarise,
         mean = as.numeric(mean(value)), sd = sd(value),
         se = sd(value)/sqrt(length(value)))
write.csv(a, file = "output/output1.csv",row.names=FALSE)
mean_plot<-ggplot(data=a,aes(x=region, y=mean,fill=region))+
  facet_wrap(~fish,ncol=4,as.table=TRUE)+
  geom_bar(stat="identity")+
  theme_bw()+
  xlab ("Region")+ ylab("Proportion")+
  theme(text=element_text(family="Times New Roman", face="bold", size=12))
png(file='figs/mean_prop.png', res=200, width=9, height=6, units ="in")
pushViewport(viewport(layout=grid.layout(1,1)))
vplayout<-function(x,y) viewport (layout.pos.row=x, layout.pos.col=y)
print(mean_plot,vp=vplayout(1,1:1))
dev.off()

#run anova for each fish to determine if the factor region is significantly different (p<0.05)
increment<-lm(prop~region, data=one)
x1<-anova(increment)
x1["fish"] <- 'fish1' 

increment<-lm(prop~region, data=two)
x2<-anova(increment)
x2["fish"] <- 'fish2' 

increment<-lm(prop~region, data=three)
x3<-anova(increment)
x3["fish"] <- 'fish3' 

increment<-lm(prop~region, data=four)
x4<-anova(increment)
x4["fish"] <- 'fish4' 

increment<-lm(prop~region, data=five)
x5<-anova(increment)
x5["fish"] <- 'fish5' 

increment<-lm(prop~region, data=six)
x6<-anova(increment)
x6["fish"] <- 'fish6'

increment<-lm(prop~region, data=seven)
x7<-anova(increment)
x7["fish"] <- 'fish7' 

increment<-lm(prop~region, data=eight)
x8<-anova(increment)
x8["fish"] <- 'fish8' 

increment<-lm(prop~region, data=nine)
x9<-anova(increment)
x9["fish"] <- 'fish9'

increment<-lm(prop~region, data=ten)
x10<-anova(increment)
x10["fish"] <- 'fish10' 

increment<-lm(prop~region, data=eleven)
x11<-anova(increment)
x11["fish"] <- 'fish11' 

increment<-lm(prop~region, data=twelve)
x12<-anova(increment)
x12["fish"] <- 'fish12'

x<-rbind(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12)
write.csv(x, "output/anova.csv") 

#funtion to run one fish's data for TukeyHSD to check multiple runs below
p_value_prop(three)

#run TukeyHSD on each anova (by fish) that had a significant region factor (p<0.05) and combine p-value output
tukey_full <- 
  list(TukeyHSD(aov(lm(prop~region, one))), 
       TukeyHSD(aov(lm(prop~region, seven))))
table_full <- NULL #object to store the overall table
for(i in 1:length(tukey_full))#loop all elements in the list
{
  tukey <- tukey_full[[i]]
  factor_table <- unlist(lapply(tukey, function(x) nrow(x)))
  factor_table <- rep(names(factor_table), factor_table)
  tukey_bound <- NULL
  for (j in 1:length(tukey))
  {
    tukey_bound <- rbind(tukey_bound, tukey[[j]])
  }
  pairs <- rownames(tukey_bound)
  rownames(tukey_bound) <- NULL
  tukey_bound <- as.data.frame(tukey_bound)
  tukey_bound$parameter <- factor_table  # add a column indicating the parameter that is compaired
  tukey_bound$pairs <- pairs# add a column with the pairs of levels that are compaired
  # add column with the model number, which corresponds to the index of the list of TukeyHSD test objects
  tukey_bound$model <- i
  table_full <- rbind(table_full, tukey_bound)
}
names(table_full) <- gsub(" ", "_", names(table_full)) #replace spaces in variables with a "-"
table_full["p_adj_0.05"] <- ifelse(table_full$p_adj<=0.05,1,0) #create variable of p_value<=0.05
table_full["fish"] <- table_full$model #change model to fish
table_full["fish"] <- ifelse(table_full$fish==2,7,1)
table_full["region_pairs"] <- table_full$pairs #change model to fish
drop <- c("model", "pairs", "X1")
table_full <- table_full[,!(names(table_full) %in% drop)] #drop the above variables from the dataset
write.csv(table_full, "output/Tukey.csv") #ouput table of all ANOVAs

#create frequency table
drop <- c("diff", "lwr", "upr","parameter", "p_adj", "fish")
table_full1 <- table_full[,!(names(table_full) %in% drop)] #drop the above variables from the dataset
melted <- melt(table_full1, id.vars=c("region_pairs"))
a<-ddply(melted, c("region_pairs"), summarise,
         count = sum(value),
         n = 12) #n is the number of fish
a["frequency"] <- as.numeric(a$count/a$n) #% of p_adj that are <=0.05 (includes all fish)
write.csv(a, file = "output/frequency.csv",row.names=FALSE)
data <- read_csv("output/frequency.csv")
a<-data[data$frequency!=0,]
drop <- c("count", "n")
a <- a[,!(names(a) %in% drop)] 
write.csv(a, file = "output/frequency_no_zeros.csv",row.names=FALSE)

#----ANALYSIS WITH ALL FISH
#run anova for each fish to determine if the factor region is significantly different (p<0.05)
df["col"] <- "fish"
df$fish <-
  paste(df$col, df$fish, sep = "")
drop <- c("col")
df <- df[,!(names(df) %in% drop)]``
increment<-lm(prop~fish+region, data=df)
x<-anova(increment)
write.csv(x, "output/anova_allfish.csv") 
increment.aov<-aov(prop~fish+region, data=df)
x<-TukeyHSD(increment.aov, "fish")
x<-as.data.frame(x$fish)
y<-TukeyHSD(increment.aov, "region")
y<-as.data.frame(y$region)
z<-rbind(x,y)
write.csv(z, "output/Tukey_allfish.csv") 


