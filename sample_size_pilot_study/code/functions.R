#to get rid of scientific notation
value_formatter <- function (val) {
  format(val, scientific =FALSE)
}
#funtion to run one fish's data
p_value_prop <- function(data){
  increment<-lm(prop~region, data=data)
  increment.aov<-aov(prop~region, data=data)
  model.tables(increment.aov, type="means")
  tk<-TukeyHSD(increment.aov, "region", odered=TRUE)
  result<-as.data.frame(tk$region)
  write.csv(result, file = "output/output2.csv",row.names=FALSE)} 

#funtion to arcsince transform a proportions
asintrans <- function(p) { 2*asin(sqrt(p)) }