setwd("/home/hyunseok/Desktop/YS/YS")
foo = read.csv("result.csv", header=F)
foo[is.na(foo)]<-0
t = 1:nrow(foo)
plot(foo[,1]~t,type='l',ylim=c(0,1))
for(i in 2:ncol(foo)){
lines(foo[,i],col=i)
}
