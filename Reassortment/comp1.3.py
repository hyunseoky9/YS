# Influenza Competition model 1.2.2
# Same concept as comp model 1 but more efficient. 
# Some key differences:
# - back mutation, tracking info, and program cut strategy parameterized.
# - No sequence information using array. Only keeps track of mutant allele amount

import numpy as np
import datetime
import timeit
from progress.bar import Bar
import sys
import os

# Parameters
# N0 = Population size
# K = Carrying capacity
# L = sequence length
# s = fitness decrease from deleterious mutation
# mu = mutation rate per site
# gen_num = generation amount 
# r = reassortment rate
# rep = repetition amount
# N1r = ratio of 1segment virus

#default parameters
back = 0
timestep = 0
krecord = 0
rep = 1
L = 300
s = 0.05
N0 = 1000
K = 1000
mu = 0.0005
gen_num = 10
cost = 0.00
r = 0.5
N1r = 0.5

if len(sys.argv) > 1: # input parameters from command line
    try:
        back = int(sys.argv[2])
        timestep = int(sys.argv[3])
        krecord = int(sys.argv[4])
        rep = int(sys.argv[5])
        L = int(sys.argv[6])
        s = float(sys.argv[7])
        N0 = int(sys.argv[8])
        K = int(sys.argv[9])
        mu = float(sys.argv[10])
        gen_num = int(sys.argv[11])
        cost = float(sys.argv[12])
        r = float(sys.argv[13])
        N1r = float(sys.argv[14])
    except:
        pass

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
    	self.id = 1
        self.k = k
    
    def mutate(self,mu):
        """
        Mutation in sequence before reproduction
        mu = mutation rate
        """
        mut_num = int(np.random.binomial(L, mu)) # number of mutation
        if not back:
            self.k += mut_num
        else:
            for i in range(mut_num):
                if np.random.uniform(0,L) < self.k:
                    self.k -= 1 # back mutation
                else:
                    self.k += 1 # normal mutation
                    
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
    	self.id = 2
        self.k1 = k1
        self.k2 = k2
        self.k = self.k1 + self.k2
        self.progeny_n = 0
    
    def mutate(self,mu):
        """
        Mutation in sequence before reproduction
        mu = mutation rate
        """
        mut_num = int(np.random.binomial(L, mu)) # number of mutation
        if not back:
            split_pt = int(np.ceil(np.random.uniform(0,1)*mut_num)) # how much of the mutation segment1 is getting
            self.k1 += split_pt
            self.k2 += mut_num - split_pt
        else:
            for i in range(mut_num):
                p = np.random.uniform(0,1)
                if np.random.uniform(0,L) < self.k: # back mutation
                    if self.k1 == 0:
                        self.k2 -= 1
                        self.k -= 1
                    elif self.k2 == 0:
                        self.k1 -= 1
                        self.k -= 1
                    else: 
                        if p < 0.5: # seg1
                            self.k1 -= 1
                            self.k -= 1
                        else: # seg2
                            self.k2 -= 1
                            self.k -= 1
                else: # normal mutation
                    if p < 0.5: # seg1
                        self.k1 += 1
                        self.k += 1
                    else: # seg2
                        self.k2 += 1
                        self.k += 1

def step(pop):
	while len(next_gen) < N0:
		next_gen = []
		popsize = [0,0]	
		sample = np.random.choice(pop, 2, replace=False)
		if sample[0].id == 1 or sample[1].id == 1:
			if np.random.uniform(0,1) < 0.5: # sample 0 or 1
				if np.random.uniform(0,1) < (1-s)**sample[0].k: # progeny live or die
					if sample[0].id == 1:
						next_gen.append(Virus1(sample[0].k))
						popsize[0] += 1
					else:
						next_gen.append(Virus2(sample[0].k1,sample[0].k2))
						popsize[1] += 1
			else:
				if np.random.uniform(0,1) > (1-s)**sample[1].k:
					if sample[1].id == 1:
						next_gen.append(Virus1(sample[1].k))
						popsize[0] += 1
					else:
						next_gen.append(Virus2(sample[1].k1,sample[1].k2))
						popsize[1] += 1
		elif sample[0].id == 2 and sample[1].id == 2:
			if np.random.uniform(0,1) < 0.5:
				if np.random.uniform(0,1) > (1-s)**(sample[0].k1+sample[1].k2):
					next_gen.append(Virus2(sample[0].k1,sample[1].k2))
					popsize[1] += 1
			else:
				if np.random.uniform(0,1) > (1-s)**(sample[0].k2+sample[1].k1):
					next_gen.append(Virus2(sample[0].k2,sample[1].k1))
					popsize[1] += 1

	if np.sum(popsize) != N0:
		raise ValueError('popsize doesn\'t add up to N0!!')
	return next_gen, popsize


start = timeit.default_timer() # timer start
bar = Bar('Processing', max=rep) # progress bar start
# write out data with file name indicating time it started collecting
now = datetime.datetime.now()
destination = 'test'
if len(sys.argv) > 1:
    try:
        destination = sys.argv[1] 
    except:
        pass
if destination not in os.listdir('./data'):
    os.system('mkdir ./data/' + destination)
params = '%d,%d,%d,%.2f,%d,%d,%.5f,%d,%.2f,%.2f,%.2f'%(back,rep,L,s,N0,K,mu,gen_num,cost,r,N1r)
tail = 'c1.3s_%s(0).csv'%(params)
while tail in os.listdir('./data/'+destination):
    lastnum = int(tail[-6])
    tail = tail[0:-6] + str(lastnum+1) + tail[-5::] 
file_name = './data/' + destination + '/' + tail
fh = open(file_name,'w')

if timestep:
    fh.write('rep,t,pop1,pop2,k1,k2\n')
    for repe in range(rep):
		# initiate
		for i in N0*N1r:
			pop.append(Virus1(0))
		for i in N0*(1-N1r):
			pop.append(Virus2(0,0))
		popsize = [N0*N1r,N0*(1-N1r)]

		# run through generation
		for gen in gen_num:
			for i in range(len(pop)):
				pop[i].mutate(mu)
			pop, popsize = step(pop)
            # recording k
            ks1 = [] # k's for each virus in a subpop
            ks2 = []
            for i in range(len(pop)):
            	if pop[i].id == 1:
            		ks1.append(pop[i].k)
            	else:
            		ks2.append(pop[i].k)
            if krecord == 1:                	
                if len(ks1) == 0:
                	ks1.append('NA')
                if len(ks2) == 0:
                	ks2.append('NA')
                ks1 = str(ks1).replace(', ','.')[1:-1]
                ks2 = str(ks2).replace(', ','.')[1:-1]
                fh.write('%d,%d,%d,%d,%s,%s\n'%(repe+1,gen+1,popsize[0],popsize[1],ks1,ks2))
            elif krecord == 0:
                if len(ks1) == 0:
                	ks1 = -1
                else:
                	ks1 = np.mean(np.array(ks1))
                if len(ks2) == 0:
                	ks2 = -1
                else krecord == 0:
                	ks2 = np.mean(np.array(ks2))
                fh.write('%d,%d,%d,%d,%.2f,%.2f\n'%(repe+1,gen+1,popsize[0],popsize[1],ks1,ks2))
            elif krecord == 2:
                if len(ks1) == 0:
                	ks1 = -1
                else:
                	ks1 = np.min(np.array(ks1))
                if len(ks2) == 0:
                	ks2 = -1
                else krecord == 0:
                	ks2 = np.min(np.array(ks2))
                fh.write('%d,%d,%d,%d,%d,%d\n'%(repe+1,gen+1,popsize[0],popsize[1],ks1, ks2))
        bar.next()

else:
    fh.write('pop1,pop2,k1,k2\n')
    for repe in range(rep):
		# initiate
		for i in N0*N1r:
			pop.append(Virus1(0))
		for i in N0*(1-N1r):
			pop.append(Virus2(0,0))
		popsize = [N0*N1r,N0*(1-N1r)]

		# run through generation
		for gen in gen_num:
			for i in range(len(pop)):
				pop[i].mutate(mu)
			pop, popsize = step(pop)

        # recording k
        ks1 = [] # k's for each virus in a subpop
        ks2 = []
        for i in range(len(pop)):
        	if pop[i].id == 1:
        		ks1.append(pop[i].k)
        	else:
        		ks2.append(pop[i].k)
        if krecord == 1:                	
            if len(ks1) == 0:
            	ks1.append('NA')
            if len(ks2) == 0:
            	ks2.append('NA')
            ks1 = str(ks1).replace(', ','.')[1:-1]
            ks2 = str(ks2).replace(', ','.')[1:-1]
            fh.write('%d,%d,%s,%s\n'%(popsize[0],popsize[1],ks1,ks2))
        elif krecord == 0:
            if len(ks1) == 0:
            	ks1 = -1
            else:
            	ks1 = np.mean(np.array(ks1))
            if len(ks2) == 0:
            	ks2 = -1
            else krecord == 0:
            	ks2 = np.mean(np.array(ks2))
            fh.write('%d,%d,%.2f,%.2f\n'%(popsize[0],popsize[1],ks1,ks2))
        elif krecord == 2:
            if len(ks1) == 0:
            	ks1 = -1
            else:
            	ks1 = np.min(np.array(ks1))
            if len(ks2) == 0:
            	ks2 = -1
            else krecord == 0:
            	ks2 = np.min(np.array(ks2))
            fh.write('%d,%d,%d,%d\n'%(popsize[0],popsize[1],ks1, ks2))
        bar.next()

fh.close()
stop = timeit.default_timer()
print('\nthe simulation ran for %.2f min'%((stop - start)/60))
bar.finish()
