rm(list=ls()) # removes all variables before starting
season = 1
change = 5
x = 0.5 #A allele freq
xs = c()
xs[1] = x
t = 200
val1= 0.8
val2 = 0.4
val3 = 0.3

vaa = val1
vAa = val2
vAA = val3
for(i in 2:t){
  if(i%%change == 1){
    if(season ==1 ){
      season = 2
      vaa = 0.3
      vAa = 0.4
      vAA = 0.8
    } else {
      season = 1
      vaa = val1
      vAa = val2
      vAA = val3
    }
  }
  vbar = vAA*x^2 + vAa*2*x*(1-x) + vaa*(1-x)^2
  xs[i] = (vAA*x^2 + vAa*x*(1-x))/vbar
  x = xs[i]
}
#print(xs)
plot(1:t,xs,type='l',ylab='frequency of A',xlab='time')
m = mean(xs[500:length(xs)])
sprintf("mean xs is %f",m)
#lines(1:t,rep(m,t),col=2)