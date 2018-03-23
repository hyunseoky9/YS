from os import system as cd
import sys
import numpy as np
rep = 100
L = 300
s = 0.05
N0 = 1000
K = 1000
mu = list(np.linspace(0.0001, 0.0009, 9))
del mu[4]
gen_num = 500
cost = 0
r = 0.5
N1r = np.linspace(0.1,0.9,9)
for j in range(len(mu)):
	for i in range(len(N1r)):
		params = '%d %d %f %d %d %.5f %d %f %f %f'%(rep,L,s,N0,K,mu[j],gen_num,cost,r,N1r[i])
		file2run = sys.argv[1]
		cd('python %s %s'%(file2run, params))
