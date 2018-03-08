from os import system as cd
import sys
import numpy as np
rep = 100
L = 300
s = 0.05
N0 = 1000
K = 1000
mu = list(np.linspace(0.0000,0.0009,10))
gen_num = 500
cost = 0
r = 0.5

for i in range(len(mu)):
	params = '%d %d %f %d %d %.5f %d %f %f'%(rep,L,s,N0,K,mu[i],gen_num,cost,r)
	file2run = sys.argv[1]
	cd('python %s %s'%(file2run, params))
