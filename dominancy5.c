/* This model seeks for heterozygosity advantage in shifting seasons
through simulating population frequency each time step. 
depending on the expression level difference due to mutation allele, the 
survival probability changes along the normal curve.
In dominancy5, we have multiple mutations rising in a rate of mutation in the loci.*/

# include <stdio.h>

float ran1(long *seed);
int genotype2index(l)

int main() {
  int survivors[] = {};
  float means[] = {};
  float  coex_ratio = {}
  int topdog = {};
  int max_alleles = 0;
  int n;
  for(n=0; n<10; n++){
    float mu = 0.001f;
    float x[] = {499/500, 1/500}; //frequency
    float *freq[] = malloc(sizeof(x)); //track record of frequency over time
    memcpy(freq,x,sizeof(x));
    int t = 5000; //generation amount
    int swit = 5; //generation amount after a season change.
    float u1 = 0.333f;
    float u2 = 0.666f;
    float sig = 0.1665; //means of normal curves for different seasons.
    int seed = -1;
    float a1 = ran1(&seed);
    printf("success!!");
  }
}
/*
for
(n in 1:10){
  mu = 0.001
  x = c(499/500,1/500) #frequency
  freq = c(x)
  t = 5000 # generation amount
  switch = 5 #generation amount after a season change.
  u1 = 0.333; u2 = 0.666; sig = 0.1665 #means of normal curves for different seasons.
  #a1 = 0.07450218; d = 0.4570124
  a1 = runif(1,0,0.5) #expression level of a wildtype
  d = runif(1, (-2*a1+u1+u2)/3, -2*a1+u1+u2)
  #d = ((-2*a1+u1+u2)/3 + (-2*a1+u1+u2))/2
  exp_d = c(0,d)
  a = c(a1, a1+d, a1+2*d)
  alleles_list = matrix(c(1,1,1,2,2,2),ncol=3)
  w1 = dnorm(a,u1,sig)/dnorm(u1,u1,sig)
  w2 = dnorm(a,u2,sig)/dnorm(u2,u2,sig)
  w = w1
  new_mut = 100 #time step when new mutation can come in

  for(i in 1:new_mut){ #loooping time step
    if(i%%switch == 1 && i>1){ # change season after specified generation time
      if(all(w==w1)){
        w = w2
      } else {
        w = w1
      }
    }
    wbar = w[1]*x[1]^2 + w[2]*2*x[1]*x[2] + w[3]*x[2]^2
    new_x = rep(0,length(x))
    new_x[2] = (w[3]*x[2]^2 + w[2]*x[1]*x[2])/wbar
    new_x[1] = (w[2]*x[1]*x[2]+w[1]*x[1]^2)/wbar
    x = new_x
    freq = cbind(freq,x)
  }

  #time steps when more mutations can occur
  for(i in (new_mut+1):t){ #loooping time step
    if(i%%switch == 1){ # change season after specified generation time
      if(all(w==w1)){
        w = w2
      } else {
        w = w1
      }
    }
    arise = runif(1,0,1)
    if (arise < mu){ # algorithm for new mutation arising
      d_new = runif(1,0,0.5) # new mutation expression dist.
      exp_d = c(exp_d,d_new)
      a = c()
      alleles_list = c() #keeps record of which index is what genotype
      #allele combination orders stored
      for(i in 1:length(exp_d)){ #refill all the expression levels
        for(j in i:length(exp_d)){
          a = c(a, a1+exp_d[i]+exp_d[j])
          alleles_list = cbind(alleles_list,c(i,j))
        }
      }
      if(all(w == w1)){
        w1 = dnorm(a,u1,sig)/dnorm(u1,u1,sig) #remake all the w1&w2
        w2 = dnorm(a,u2,sig)/dnorm(u2,u2,sig)
        w = w1 
      } else {
        w1 = dnorm(a,u1,sig)/dnorm(u1,u1,sig)
        w2 = dnorm(a,u2,sig)/dnorm(u2,u2,sig)
        w = w2 
      }
      arise_from = ceiling(runif(1,0,length(x)))
      x = c(x,1/500)
      x[arise_from] = x[arise_from] - 1/500
      freq = rbind(freq,rep(0,dim(freq)[2]))
      freq = cbind(freq,x) 
    }
    x_square = c()
    for(j in 1:ncol(alleles_list)){ #get all the factors when x is squared
      if(alleles_list[1,j]==alleles_list[2,j]){
        x_square = c(x_square, x[alleles_list[1,j]]*x[alleles_list[2,j]])
      } else{
        x_square = c(x_square, 2*x[alleles_list[1,j]]*x[alleles_list[2,j]])
      }
    }
    wbar = sum(x_square*w)
    new_x = rep(0,length(x))
    for(j in 1:length(x)){ #update each frequency in x
      factor_sum = 0
      for(k in 1:ncol(alleles_list)){ #for every genotype
        if(all(alleles_list[,k]==j)){
          factor_sum = factor_sum + x_square[k]*w[k]
        } else if(any(alleles_list[,k]==j)){
          factor_sum = factor_sum + 0.5*x_square[k]*w[k]
        }
      }
      new_x[j] = factor_sum/wbar
    }
    x = new_x
    freq = cbind(freq,x)
  }
  avg = c() #average frequency for each allele in the last 50 gen.
  for( row in 1:nrow(freq)){
    avg = c(avg,mean(freq[row,(ncol(freq)-50):ncol(freq)]))
  }
  extinct = which(avg < 0.01)
  if (length(extinct)>0){
    survivors[[n]] = seq(1,length(x))[-extinct] # surviving alleles in each run
  } else if (length(extinct)==0){
    survivors[[n]] = seq(1,length(x))
    coex_ratio[[length(coex_ratio)+1]] = cbind(w1,w2,w1*w2) # coex_ratio stores w1,w2, and w1w2 when all alleles coexist. 
  }
  means[[n]] = w1*w2 #ww for each genotype
  if(ncol(alleles_list)<max_alleles){ # update max_alleles
    max_alleles = ncol(alleles_list)
  }
  allele_label = c() # make allele labels
  for( i in 1:ncol(alleles_list)){
    label = sprintf('A%s%s',alleles_list[1,i],alleles_list[2,i])
    allele_label = c(allele_label, label)
  }
  par(mfrow=c(2,2))
  par(mar=c(2,2,2,2))
  title = sprintf('bl=1,r=2,g=3')
  time_line = 100:5000
  plot(time_line,freq[1,time_line], ylim=c(0,1), type='l',main = title) #plot of frequency
  for(i in 2:length(x)){ # graph all the frequencies
    lines(time_line,freq[i,time_line],type='l',col=i)
  }
  ww = w1*w2
  plot(ww, axes=FALSE, ylim=c(0,max(ww)), col=1) #plot of w1*w2
  axis(2)
  axis(1, at=seq_along(ww),labels=allele_label)
  title('geometric mean of two seasons (w1*w2)')
  plot(w1,main='w1')
  plot(w2,main='w2')
  
}
*/









#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran1(long *idum)
{
  int j;
  long k;
  static long iy=0;
  static long iv[NTAB];
  float temp;

  if (*idum <= 0 || !iy) {
    if (-(*idum) < 1) *idum=1;
    else *idum = -(*idum);
    for (j=NTAB+7;j>=0;j--) {
      k=(*idum)/IQ;
      *idum=IA*(*idum-k*IQ)-IR*k;
      if (*idum < 0) *idum += IM;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ;
  *idum=IA*(*idum-k*IQ)-IR*k;
  if (*idum < 0) *idum += IM;
  j=iy/NDIV;
  iy=iv[j];
  iv[j] = *idum;
  if ((temp=AM*iy) > RNMX) return RNMX;
  else return temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX