rm(list=ls())
season = 1
change = 5
x = 0.3 #A allele freq
xs = c()
xs[1] = x
t = 500
vaa = 0.84
vAa = 0.8
vAA = 0.9
for(i in 2:t){
  if(i%%change == 1){
    if(season ==1 ){
      season = 2
      vaa = 0.84
      vAa = 0.8
      vAA = 0.9
    } else {
      season = 1
      vaa = 0.84
      vAa = 0.8
      vAA = 0.9
    }
  }
  vbar = vAA*x^2 + vAa*2*x*(1-x) + vaa*(1-x)^2
  xs[i] = (vAA*x^2 + vAa*x*(1-x))/vbar
  x = xs[i]
}
print(xs)
plot(1:500,xs)
sprintf("mean xs is %f",mean(xs))