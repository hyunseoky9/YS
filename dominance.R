# parameters
d=0.25
n = 100
adv = 0
u1 = 0.25; u2 = 0.75
sig = 0.125
ratio = c()
#distribution graph
#x = seq(0,1,0.01)
#y1 = dnorm(x,u1,sig)
#y2 = dnorm(x,u2,sig)
#plot(y1~x,type='l')
#lines(y2~x)
ratios = list()
list_num = 1
par(mfrow=c(3,3))
for (d in seq(0.05,0.45,0.05)){
  for (j in c(1:100)){
    for (i in c(1:n)){
      a11 = runif(1, min = 0, max = (1-2*d))
      a = c(a11, a11+d, a11+(2*d))
      w1 = dnorm(a, u1, sig)
      w2 = dnorm(a, u2, sig)
      ww = w1*w2
      if (ww[1] < ww[2] & ww[2] > ww[3]){
        adv = adv + 1
      }
    }
  ratio[j] = adv/n
  adv = 0
  }
ratios[[list_num]] = ratio
list_num = list_num +1
hist(ratio, main=NULL)
title(sprintf('dist. of htzg. adv freq. d=%f',d))
}


