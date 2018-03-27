# Influenza Competition model 1.2.2
# Same concept as comp model 1 but more efficient. 
# Some key differences:
# - IT HAS BACK MUTATION (NO BACKMUTATION IN 1.2)
# - No sequence information using array. Only keeps track of mutant allele amount

import numpy as np
import datetime
import timeit
from progress.bar import Bar
import sys
import os

# Parameters
# L = sequence length for a virus
# N = Population size (N0 = initial population)
# K = Carrying capacity
# s = fitness decrease from deleterious mutation
# mu = mutation rate
# gen_num = generation amount 
# r = reassortment rate
# rep = repetition amount
# N1r = ratio of 1segment virus

#default parameters
rep = 100
L = 300
s = 0.05
N0 = 1000
K = 1000
mu = 0.0005
gen_num = 500
cost = 0.00
r = 0.5
N1r = 0.5

class Virus1():
    """
    This class produces objects which are single agents of a influenza virus with 1 segment.
    k = number of deleterious mutation in entire genome
    s = fitness decrease from deleterious mutation
    L = sequence length for a virus
    cost = cost of having multisegments
    w = fitness
    """
    def __init__(self,k):
        self.k = k
        self.w = (1 - s)**self.k
    
    def mutate(self,mu):
        """
        Mutation in sequence before reproduction
        mu = mutation rate
        """
        mut_num = int(np.random.binomial(L, mu)) # number of mutation
        self.k += mut_num
                
                
class Virus2():
    """
    This class produces objects which are single agents of a influenza virus with 2 segments.
    k1 = number of deleterious mutation in segment1
    k2 = number of deleterious mutation in segment2
    k = number of deleterious mutation in entire genome
    cost = cost of having multisegments
    w = fitness
    progeny_n = number of progenies a virus agent will have during reproduction. Default is 0.
    """
    def __init__(self,k1, k2):
        self.k1 = k1
        self.k2 = k2
        self.k = self.k1 + self.k2
        self.progeny_n = 0
        self.w = (1 - s)**self.k - cost
    
    def mutate(self,mu):
        """
        Mutation in sequence before reproduction
        mu = mutation rate
        """
        mut_num = int(np.random.binomial(L, mu)) # number of mutation
        split_pt = int(np.ceil(np.random.uniform(0,1)*mut_num)) # how much of the mutation segment1 is getting
        self.k1 += split_pt
        self.k2 += mut_num - split_pt

def reproduce(viruses1, viruses2):
    """
    input
        viruses1 = array of 1segmented virus agents.
        viruses2 = array of multi-segmented virus agents.
    output
        next_gen1 = array of 1segemented virus of next generation.
        next_gen2 = array of multi-segemented virus of next generation.
    """
    next_gen1 = []
    next_gen2 = []
    
    for i in range(len(viruses1)): # reproduction of 1segs
        progeny_n = np.random.poisson(viruses1[i].w*(2/(1+N/K))) # number of progeny for that virus
        for j in range(progeny_n):
            next_gen1.append(Virus1(viruses1[i].k))
            
    remaining = [] # virus parents with remaining progenies to reproduce (for multi-segs)
    for i in range(len(viruses2)): # reproduction of multi-segs
        viruses2[i].progeny_n = np.random.poisson(viruses2[i].w*(2/(1+N/K)))
        if viruses2[i].progeny_n > 0:
            remaining.append(i)
    while len(remaining) >= 2: # while there's at least 2 remaining viruses to reproduce together
        samp = sorted(np.random.choice(remaining, 2, replace=False)) # pick 2 from the remaining
        offspring1, offspring2 = reassort(viruses2[samp[0]], viruses2[samp[1]])
        next_gen2.append(offspring1)
        next_gen2.append(offspring2)
        viruses2[samp[0]].progeny_n -= 1
        viruses2[samp[1]].progeny_n -= 1
        if viruses2[samp[0]].progeny_n == 0: # remove index from remaining if progeny_n is 0.
            del remaining[np.where(np.array(remaining)==samp[0])[0][0]]
        if viruses2[samp[1]].progeny_n == 0:
            del remaining[np.where(np.array(remaining)==samp[1])[0][0]]
    # when there's only 1 index reamining, go through regular reproduction (no reassortment)
    for i in remaining:
        for j in range(viruses2[i].progeny_n):
            next_gen2.append(Virus2(viruses2[i].k1, viruses2[i].k2))
    return next_gen1, next_gen2

def reassort(v1, v2):
    prob = np.random.uniform(0,1)
    if prob < r: # reassortment
        offspring1 = Virus2(v1.k1, v2.k2) # v1 gives first segment
        offspring2 = Virus2(v2.k1, v1.k2) # v2 gives first segment
        return offspring1, offspring2
    else: # no reassortment
        offspring1 = Virus2(v1.k1,v1.k2)
        offspring2 = Virus2(v2.k1,v2.k2)
        return offspring1, offspring2



if len(sys.argv) > 1: # input parameters from command line
    try:
        rep = int(sys.argv[1])
        L = int(sys.argv[2])
        s = float(sys.argv[3])
        N0 = int(sys.argv[4])
        K = int(sys.argv[5])
        mu = float(sys.argv[6])
        gen_num = int(sys.argv[7])
        cost = float(sys.argv[8])
        r = float(sys.argv[9])
        N1r = float(sys.argv[10])
    except:
        pass


start = timeit.default_timer() # timer start
bar = Bar('Processing', max=rep) # progress bar start
# write out data with file name indicating time it started collecting
now = datetime.datetime.now()
destination = input('Enter file destination?')
if destination not in os.listdir('./data'):
    os.command('mkdir ' + destination)
params = '%d,%d,%.2f,%d,%d,%.5f,%d,%.2f,%.2f,%.2f'%(rep,L,s,N0,K,mu,gen_num,cost,r,N1r)
file_name = './data/'+destination+'/c1.2.2s_%s(0).csv'%(params)
while file_name[7::] in os.listdir('./data'):
    lastnum = int(file_name[-6])
    file_name = file_name[0:-6] + str(lastnum+1) + file_name[-5::]
fh = open(file_name,'w')
detail = input('Do you want to record the population of each time step? 0=no 1=yes (Will only record last step if typed "no")')
if detail: 
    detail2 = input('Do you want to record each individual\'s k? 0=no 1=yes (Will only record means if "no"')
if detail:
    fh.write('rep,t,pop1,pop2,k1,k2\n')
    for repe in range(rep):
        # Initialize the virus population. Each subpopulation gets half the population.
        N = N0
        viruses1 = []
        viruses2 = []
        for i in range(int(N*(1-N1r))):
            viruses2.append(Virus2(0,0))
        for i in range(int(N*N1r)):
            viruses1.append(Virus1(0))
        fh.write('%d,%d,%d,%d,%.2f,%.2f\n'%(repe+1,0,len(viruses1),len(viruses2),0,0))

        # run the simulation
        for gen in range(gen_num):
            if len(viruses1) == 0 or len(viruses2) == 0: # terminate if either subpop -> 0
                break
            for i in range(len(viruses1)):
                viruses1[i].mutate(mu)
            for j in range(len(viruses2)):
                viruses2[j].mutate(mu)
            viruses1, viruses2 = reproduce(viruses1,viruses2)

            # get kmean for each subpop
            if not detail2:
                ks = [] # k's for each virus in a subpop
                if len(viruses1)>0:
                    for i in range(len(viruses1)):
                        ks.append(viruses1[i].k)
                    k_means1 = np.mean(np.array(ks))
                else:
                    k_means1 = -1
                ks = []
                if len(viruses2)>0:
                    for j in range(len(viruses2)):
                        ks.append(viruses2[j].k)
                    k_means2 = np.mean(np.array(ks))
                else:
                    k_means2 = -1

                fh.write('%d,%d,%d,%d,%.2f,%.2f\n'%(repe+1,gen+1,len(viruses1),len(viruses2),k_means1, k_means2))
                N = len(viruses1) + len(viruses2)
            else:
                ks1 = [] # k's for each virus in a subpop
                if len(viruses1)>0:
                    for i in range(len(viruses1)):
                        ks1.append(viruses1[i].k)
                else:
                    ks1.append('NA')
                ks2 = []
                if len(viruses2)>0:
                    for j in range(len(viruses2)):
                        ks.append(viruses2[j].k)
                else:
                    ks2.append('NA')
                fh.write('%d,%d,%d,%d,%s,%s\n'%(repe+1,gen+1,len(viruses1),len(viruses2),str(ks1),str(ks2)))
                N = len(viruses1) + len(viruses2)
        bar.next()

else:
    fh.write('pop1,pop2,k1,k2\n')
    for repe in range(rep):
        # Initialize the virus population. Each subpopulation gets half the population.
        N = N0
        viruses1 = []
        viruses2 = []
        for i in range(int(N*(1-N1r))):
            viruses2.append(Virus2(0,0))
        for i in range(int(N*N1r)):
            viruses1.append(Virus1(0))

        # run the simulation
        for gen in range(gen_num):
            if len(viruses1) == 0 or len(viruses2) == 0: # terminate if either subpop -> 0
                break
            for i in range(len(viruses1)):
                viruses1[i].mutate(mu)
            for j in range(len(viruses2)):
                viruses2[j].mutate(mu)
            viruses1, viruses2 = reproduce(viruses1,viruses2)
            N = len(viruses1) + len(viruses2)

        # get kmean for each subpop
        ks = [] # k's for each virus in a subpop
        if len(viruses1)>0:
            for i in range(len(viruses1)):
                ks.append(viruses1[i].k)
            k_means1 = np.mean(np.array(ks))
        else:
            k_means1 = -1
        ks = []
        if len(viruses2)>0:
            for j in range(len(viruses2)):
                ks.append(viruses2[j].k)
            k_means2 = np.mean(np.array(ks))
        else:
            k_means2 = -1

        fh.write('%d,%d,%.2f,%.2f\n'%(len(viruses1),len(viruses2),
                                     k_means1, k_means2)) 
        bar.next()

fh.close()
stop = timeit.default_timer()
print('\nthe simulation ran for %.2f min'%((stop - start)/60))
bar.finish()
