# This model seeks for heterozygosity advantage in shifting seasons
# through simulating population frequency each time step. 
# depending on the expression level difference due to mutation allele, the 
# survival probability changes along the normal curve.
#set.seed(1001)

x = c(499/500,1/500) #frequency
freq = c(x[2])
t=1000 # generation amount
s = 1 # season either 1 or 2.
switch = 5 #generation amount after a season change.
u1 = 0.333; u2 = 0.666; sig = 0.1665 #means of normal curves for different seasons.
a1 = runif(1,0,0.5) #expression level of a wildtype
d = ((-2*a1+u1+u2)/3 + (-2*a1+u1+u2))/2
#a1 = 0.33; d = 0.17;
a = c(a1, a1+d, a1+2*d)
w1 = dnorm(a,u1,sig)/dnorm(u1,u1,sig)
w2 = dnorm(a,u2,sig)/dnorm(u2,u2,sig)
w = w1

for(i in 1:t){ #loooping time step
  if(i%%switch == 1 && i>1){ # change season after specified generation time
    if(s==1){
      s=2
      w = w2
    } else {
      s=1
      w = w1
    }
  }
  wbar = w[3]*x[2]^2 + w[2]*2*x[1]*x[2] + w[1]*x[1]^2
  x[2] = (w[3]*x[2]^2 + w[2]*x[1]*x[2])/wbar
  x[1] = 1 - x[2]
  freq = c(freq,x[2])
}
par(mfrow=c(2,1))
par(mar=c(2,2,2,2))
title = sprintf('d= %f, a1= %f',d,a1)
plot(freq, type='l',main = title)
plot(w1,ylim=c(0,1),col=1)
points(w2,col=2)
legend("topright",legend=c('season1 w','season2 w'),pch=1,col=c("black",'red'), cex=0.8)
print(w1*w2)
