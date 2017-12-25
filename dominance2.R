# parameters
n = 100
adv = 0
u1 = 0.25; u2 = 0.75
sig = 0.125
ratio = c()
#distribution graph
x = seq(0,5,0.01)
#y1 = dnorm(x,u1,sig)
#y2 = dnorm(x,u2,sig)
y3 = dexp(x, rate=5)
plot(y3~x,type='l')
#lines(y2~x)
ratios = list()
list_num = 1
for (j in c(1:100)){
  for (i in c(1:n)){
    d = rexp(1, rate=5)
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
hist(ratio, main=NULL, xlim=c(0,1))
title(sprintf('dist. of htzg. adv freq.'))
