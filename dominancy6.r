#  This model seeks for heterozygosity advantage in shifting seasons
# through simulating population frequency each time step. 
# depending on the expression level difference due to mutation allele, the 
# survival probability changes along the normal curve.
# In dominancy6, we have multiple mutations rising in a rate of
# mutation in the loci. In addition the frequency now has stochastic
# element through WF model simulation.

rm(list=ls())

genotype2index <- function(genotype, new_al_num) {
# converts the genotype to the index in new population array 
  if(min(genotype)==1){
    index = max(genotype)
  } else {
    index = sum(seq(new_al_num,1)[1:(min(genotype)-1)]) +
    abs(genotype[2]-genotype[1]) + 1
  }
  return(index)
}
# in stat in each round, who survived, ww order.
survivors = list()
means = list()
coex_ratio = list()
topdog = c()
max_alleles = 0

for(n in 1:10){
  N = 100
  mu = 0.001
  x = c(1) #frequency of alleles
  pop = matrix(c(N),ncol=1,nrow=1) #array of population of genotypes
  freq = matrix(x,ncol=1,nrow=1)
  t = 5000 # generation amount
  switch = 5 #generation amount after a season change.
  u1 = 0.333; u2 = 0.666; sig = 0.1665 #means of normal curves for different seasons.
  #a1 = 0.07450218; d = 0.4570124
  a1 = runif(1,0,0.5) #expression level of a wildtype
  all_exp = c(a1) #expression level of each allele
  a = c(2*a1) #expression level of all genotypes
  genotypes_list = matrix(c(1,1),ncol=1)
  w1 = dnorm(a,u1,sig)/dnorm(u1,u1,sig)
  w2 = dnorm(a,u2,sig)/dnorm(u2,u2,sig)
  w = w1

  for(time in 1:t){ #loooping time step
    if(time%%switch == 1 && time>1){ # change season after specified generation time
      if(all(w==w1)){
        w = w2
      } else {
        w = w1
      }
    }
    arise = runif(1)
    #assumption: mutation arises only once per generation
    #protocol for when mutation arises
    if (arise < mu){ # algorithm for new mutation arising
      positive_pop = which(pop>0) #select genotypes that have counts
      select_mutant = positive_pop[ceiling(runif(1)*length(positive_pop))] #select the genotype that will mutate out of the genotypes with positive counts
      mutant_arisen_genotype = genotypes_list[,select_mutant]
      arise_from = mutant_arisen_genotype[ceiling(runif(1)*length(mutant_arisen_genotype))] #allele number that the mutation arose from
      other_allele = mutant_arisen_genotype[which(mutant_arisen_genotype != arise_from)]
      mutant_genotype = sort(c(length(x)+1,other_allele))
      d_new = runif(1,-0.5,0.5) # new mutation expression dist.
      new_exp = all_exp[arise_from]+d_new
      if (new_exp<0){
        new_exp = 0
      }
      all_exp = c(all_exp,new_exp)
      newpop = rep(0, sum(seq(1,length(all_exp)))) #new population with added allele
      j=1
      tempop = pop
      for(i in length(x):1){ #put old pop into new pop count
        newpop[j:(j+i-1)] = tempop[1:i]
        tempop = tempop[-(1:i)]
        j = j + (i-1) + 2
      }
      #find index of the mutant genotypes in the newpop
      from_index = genotype2index(mutant_arisen_genotype,length(x)+1)
      to_index = genotype2index(mutant_genotype,length(x)+1)
      new_genotypes_list = c() #keeps record of which index is what genotype
      #allele combination orders stored
      a = c() #reset genotype expression level
      #refill all the expression levels for genotypes and 
      for(i in 1:length(all_exp)){
        for(j in i:length(all_exp)){
          a = c(a,all_exp[i]+all_exp[j])
          new_genotypes_list = cbind(new_genotypes_list,c(i,j))
        }
      }
      # add mutation to the newpop
      newpop[from_index] = newpop[from_index] - 1
      newpop[to_index] = newpop[to_index] + 1
      genotypes_list = new_genotypes_list
      if(all(w == w1)){
        w1 = dnorm(a,u1,sig)/dnorm(u1,u1,sig) #remake all the w1&w2
        w2 = dnorm(a,u2,sig)/dnorm(u2,u2,sig)
        w = w1 
      } else {
        w1 = dnorm(a,u1,sig)/dnorm(u1,u1,sig)
        w2 = dnorm(a,u2,sig)/dnorm(u2,u2,sig)
        w = w2 
      }
      pop = newpop
      x[arise_from] = x[arise_from] - 1/(2*N)
      x = c(x,1/(2*N))
      freq = rbind(freq,rep(0,dim(freq)[2]))
    }
    x_square = c()
    for(j in 1:ncol(genotypes_list)){ #get all the factors when x is squared
      if(genotypes_list[1,j] == genotypes_list[2,j]){
        x_square = c(x_square, x[genotypes_list[1,j]]*x[genotypes_list[2,j]])
      } else {
        x_square = c(x_square, 2*x[genotypes_list[1,j]]*x[genotypes_list[2,j]])
      }
    }
    wbar = sum(x_square*w)
    pop = rmultinom(1, size=N, x_square/wbar*w) #stochastic reproduction
    for(j in 1:length(x)){ #calculate allele frequencies for the next gen.
      factor_sum = 0
      for(i in 1:ncol(genotypes_list)){
        if(any(genotypes_list[,i] == j)){
          if(genotypes_list[1,i] == genotypes_list[2,i]){
            factor_sum = factor_sum + pop[i]
          } else {
            factor_sum = factor_sum + pop[i]/2
          }
        }
      }
      x[j] = factor_sum/N
    }
    freq = cbind(freq,x) # update freq
  }
}
  extinct = which(freq[,ncol(freq)] == 0)
  if (length(extinct)>0){
    survivors[[n]] = seq(1,length(x))[-extinct] # surviving alleles in each run in 'n'.
    ww = w1[-extinct]*w2[-extinct]
  } else if (length(extinct)==0){
    survivors[[n]] = seq(1,length(x))
    coex_ratio[[length(coex_ratio)+1]] = cbind(w1,w2,w1*w2) # coex_ratio stores w1,w2, and w1w2 when all alleles coexist. 
    ww = w1*w2
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
  time_line = 1:5000
  plot(time_line,freq[1,time_line], ylim=c(0,1), type='l',main = title) #plot of frequency
  for(i in 2:length(x)){ # graph all the frequencies
    lines(time_line,freq[i,time_line],type='l',col=i)
  }
  plot(ww, axes=FALSE, ylim=c(0,max(ww)), col=1) #plot of w1*w2 only the surviving alleles in the end
  axis(2)
  axis(1, at=seq_along(ww),labels=allele_label)
  title('geometric mean of two seasons (w1*w2)')
  plot(w1,main='w1') #plot of w1 and w2 (only the surviving ones).
  plot(w2,main='w2')
}









genotype2index <- function(genotype, new_al_num) {
# converts the genotype to the index in new population array 
  if(min(genotype)==1){
    index = max(genotype)
  } else {
    index = sum(seq(new_al_num,1)[1:(min(genotype)-1)]) +
    abs(genotype[2]-genotype[1]) + 1
  }
  return(index)
}