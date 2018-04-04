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
indivk = 0
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
        indivk = int(sys.argv[4])
        untilext = int(sys.argv[5])

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
        self.k = k
        self.w = (1 - s)**self.k
    
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
        if not back:
            split_pt = int(np.ceil(np.random.uniform(0,1)*mut_num)) # how much of the mutation segment1 is getting
            self.k1 += split_pt
            self.k2 += mut_num - split_pt
        else:
            for i in range(mut_num):
                p = np.random.uniform(0,1)
                if np.random.uniform(0,L) < self.k: # back mutation
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

