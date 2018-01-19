# This model seeks for heterozygosity advantage in shifting seasons
# through simulating population frequency each time step. 
# depending on the expression level difference due to mutation allele, the 
# survival probability changes along the normal curve.
# In dominancy4, a new, third mutation is introduced to the population
# at time stip 100 with a random change in expression level. As a spinoff
# it calculates each allele's fitness based on genotype's geometric mean.
set.seed(1009)
rm(list=ls())
# in stat in each round, who survived, ww order.
hypright = 0 # number of times my hypothesis is right.
hypwrong = c() #trial index of a instance where hypo is wrong.
coex_find = 0 # coexistence find mode. put 1 to turn on.
coex_freq = c() # records frequency history if coexisting at the end
survivors = list() # list of survivor allele each trial
means = c() # matrix of geometric mean of w's. Each column is a trial
coex_ratio = list() # list that stores w1,w2,ww when all 3 alleles coexist
coex_index = c() #array of trial index where all 3 alleles survive
allele_fitnesses = c() # matrix of allele fitnesses each column is a trial
initfreq = c() # matrix of initial frequency when 3rd allele come in
n = 1
for(trial in 1:n){
  x = c(499/500,1/500) #frequency
  freq = c()
  t = 2000 # generation amount
  s = 1 # season either 1 or 2.
  switch = 5 #generation amount after a season change.
  u1 = 0.4; u2 = 0.6; sig = 0.1665; #means of normal curves for different seasons.
  #a1 = 0.07450218; d = 0.4570124
  a1 = runif(1,0, 0.5) #expression level of a wildtype
  d = runif(1, 0, 0.5)
  #d = ((-2*a1+u1+u2)/3 + (-2*a1+u1+u2))/2
  exp_d = c(a1,d)
  a = c(2*a1, a1+d, 2*d)
  w1 = dnorm(a,u1,sig)/dnorm(u1,u1,sig)
  w2 = dnorm(a,u2,sig)/dnorm(u2,u2,sig)
  w = w1
  new_mut = 1 #time step when new mutation comes in


  #for(i in 1:new_mut){ #loooping time step
  #  if(i%%switch == 1 && i>1){ # change season after specified generation time
  #    if(s==1){
  #      s=2 
  #      w = w2
  #    } else {
  #      s=1
  #      w = w1
  #    }
  #  }
  #  wbar = w[1]*x[1]^2 + w[2]*2*x[1]*x[2] + w[3]*x[2]^2
  # new_x = rep(0,length(x))
  # new_x[2] = (w[3]*x[2]^2 + w[2]*x[1]*x[2])/wbar
  #  new_x[1] = (w[2]*x[1]*x[2]+w[1]*x[1]^2)/wbar
  #  x = new_x
  #  freq = cbind(freq,x)
  #}
  
  #NEW MUTATION
  d2 = runif(1,0,0.5) # new mutation expression dist.
  exp_d = c(exp_d,d2)
  a = c()
  alleles_list = c() #keeps record of which index is what genotype
  #allele combination orders stored
  for(i in 1:length(exp_d)){ #refill all the expression levels
    for(j in i:length(exp_d)){
      a = c(a, exp_d[i]+exp_d[j])
      alleles_list = cbind(alleles_list,c(i,j))
    }
  }
  w1 = dnorm(a,u1,sig)/dnorm(u1,u1,sig) #remake all the w1&w2
  w2 = dnorm(a,u2,sig)/dnorm(u2,u2,sig)
  
  # PERTURBATION EXPERMIENT
  w2[4] = w2[4]  + 0.2
  
  w2[2] = w2[2] + 0.4 #coexistence
  #w2[1] = w2[1] + 0.2
  #w1[1] = w1[1] + 0.2
  
  #w1[1]=w1[1] +0.2
  
  
  w = w1
  
  
  # same frequency for all
  x = c(0.3333334,0.3333333,0.3333333)
  #put all the past freq for third allele (which is 0)
  freq = cbind(freq,x)
  initfreq = cbind(initfreq,x)
  #time steps after 3rd mutation
  for(i in 1:t){ #loooping time step
    if(i%%switch == 1){ # change season after specified generation time
      if(s==1){
        s = 2
        w = w2
      } else {
        s = 1
        w = w1
      }
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
    allele_label = c() # make allele labels
  for( i in 1:ncol(alleles_list)){
    label = sprintf('A%s%s',alleles_list[1,i],alleles_list[2,i])
    allele_label = c(allele_label, label)
  }
  ww = w1*w2 # geometric mean
  wamean = (w1+w2)/2 # arithmetic mean
  
  # calculate allele fitness.
  allele_fitness = rep(0,length(x))
  for(j in 1:ncol(alleles_list)){
    if(alleles_list[1,j] == alleles_list[2,j]){
      allele_fitness[alleles_list[1,j]] = allele_fitness[alleles_list[1,j]] + ww[j]
    } else {
      allele_fitness[alleles_list[1,j]] = allele_fitness[alleles_list[1,j]] + ww[j]/2
      allele_fitness[alleles_list[2,j]] = allele_fitness[alleles_list[2,j]] + ww[j]/2
    }
  }
  allele_fitnesses = cbind(allele_fitnesses,allele_fitness)
  
  avg = c()
  for( row in 1:nrow(freq)){
    avg = c(avg,mean(freq[row,(ncol(freq)-50):ncol(freq)]))
  }
  extinct = which(avg < 0.01)
  if (length(extinct)>0){
    survivors[[trial]] = c(1,2,3)[-extinct]
  } else if (length(extinct)==0){
    survivors[[trial]] = c(1,2,3)
    coex_ratio[[length(coex_ratio)+1]] = cbind(w1,w2,w1*w2)
    coex_index = c(coex_index,trial)
  }
  means = cbind(means,w1*w2)
  
  #result output for a single instance
  if(coex_find){
    if(length(extinct)==0){
      par(mfrow=c(2,3))
      par(mar=c(2,2,2,2))
      title = sprintf('bl=1,r=2,g=3 trial#=%d',trial)
      time_line = 500:t
      plot(time_line,freq[1,time_line], ylim=c(0,1), type='l',main = title) #plot of frequency
      lines(time_line,freq[2,time_line],type='l',col=2)
      lines(time_line,freq[3,time_line],type='l',col=3)
      plot(ww, axes=FALSE, ylim=c(0,max(ww)), col=1) #plot of w1*w2
      axis(2)
      axis(1, at=seq_along(ww),labels=allele_label)
      title('geometric mean')
      plot(w1,main='w1') # plot w1
      plot(w2,main='w2') # plot w2
      plot(wamean, axes=FALSE, ylim=c(0,max(wamean)),col=1) # plot mean of w1 and w2
      axis(2)
      axis(1, at=seq_along(wamean),labels= allele_label)
      title('arithmetic mean')
      plot(seq(0,1,0.05),dnorm(seq(0,1,0.05),u1,sig)/dnorm(u2,u2,sig),type='l',col=1)
      lines(seq(0,1,0.05),dnorm(seq(0,1,0.05),u2,sig)/dnorm(u2,u2,sig),type='l',col=2)
      points(a,w1,col=1)
      points(a,w2,col=2)
      text(a,w2, labels = allele_label, pos = 4)
      coex_freq = cbind(coex_freq, freq)
    }
  } else {
    par(mfrow=c(2,3))
    par(mar=c(2,2,2,2))
    title = sprintf('bl=1,r=2,g=3 trial#=%d',trial)
    time_line = 500:t
    plot(time_line,freq[1,time_line], ylim=c(0,1), type='l',main = title) #plot of frequency
    lines(time_line,freq[2,time_line],type='l',col=2)
    lines(time_line,freq[3,time_line],type='l',col=3)
    plot(ww, axes=FALSE, ylim=c(0,max(ww)), col=1) #plot of w1*w2
    axis(2)
    axis(1, at=seq_along(ww),labels=allele_label)
    title('geometric mean')
    plot(w1,main='w1') # plot w1
    plot(w2,main='w2') # plot w2
    plot(wamean, axes=FALSE, ylim=c(0,max(wamean)),col=1) # plot mean of w1 and w2
    axis(2)
    axis(1, at=seq_along(wamean),labels= allele_label)
    title('arithmetic mean')
    plot(seq(0,1,0.05),dnorm(seq(0,1,0.05),u1,sig)/dnorm(u2,u2,sig),type='l',col=1)
    lines(seq(0,1,0.05),dnorm(seq(0,1,0.05),u2,sig)/dnorm(u2,u2,sig),type='l',col=2)
    points(a,w1,col=1)
    points(a,w2,col=2)
    text(a,w2, labels = allele_label, pos = 4)
  }
  
  checker = hypright
  top2 = order(ww, decreasing=TRUE)[1:2]
  tophomo = order(ww[c(1,4,6)])[3]
  if(tophomo == 3){
    tophomo_index = 6
    tophet_index = c(3,5)
  } else if(tophomo ==2){
    tophomo_index = 4
    tophet_index = c(2,5)
  } else {
    tophomo_index = 1
    tophet_index = c(2,3)
  }
  if(ww[tophet_index[1]]>ww[tophet_index[2]]){
    tophet_index = tophet_index[1]
  } else {
    tophet_index = tophet_index[2]
  }
  if(ww[tophet_index] > ww[tophomo_index]){
    hypsurvive = alleles_list[,tophet_index]
  } else {
    hypsurvive = tophomo
  }
  if(all(hypsurvive == survivors[[trial]])){
    hypright = hypright + 1
  }
  
  if(checker == hypright){
    hypwrong = c(hypwrong,trial)
  }
}

topdog = c()
for(col in 1:ncol(means)){
  top = alleles_list[,which(means[,col]==max(means[,col]))]
  topdog = cbind(topdog, top)
}
colnames(topdog)<-c()
#result output for a stat of multiple instances
counts = rep(0,3)
for(i in 1:length(survivors)){
  for(j in 1:length(survivors[[i]])){
   counts[survivors[[i]][j]] = counts[survivors[[i]][j]] + 1
  }
}
rownames(means) = allele_label
ratios = counts/sum(counts) # ratios of alleles surviving in all trials
names(ratios)<-c(1,2,3)
#print(ratios)
print(survivors)
colnames(allele_fitnesses)<-c()
print(allele_fitnesses)
print(means)
