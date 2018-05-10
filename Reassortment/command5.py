# was used for making drift_test data.

from os import system as cd
import sys
import numpy as np
seed = np.random.randint(-9223372036854775808,-1)
back = 0
timestep = 1
krecord = 0 #  2= smallest k of the subpop. 1= all indiv's k. 0=mean k.
untilext = 0
rep = 50
L = 300
s = [0.9,0.1,0.01,0.001]
N0 = [50,100]
K = [50,100]
mu = [0.001/300,0.01/300,0.1/300,1/300]
gen_num = 500
cost = 0
r = 0.5
N1r = [0,1]
destination = 'felsenstein_test_1.3c'
file2run = sys.argv[1]
version = file2run[4:7]

count = 0
for i in range(len(mu)):
	for k in range(len(s)):
		for j in range(len(N0)):
			for l in range(len(N1r)):
				if version == '1.2':
					params = '%s %d %d %d %d %d %d %f %d %d %.5f %d %f %f %f'%(destination,back,timestep,krecord,untilext,
																				rep,L,s[k],N0[j],K[j],mu[i],gen_num,cost,r,N1r[l])
				elif version == '1.3':
					params = '%s %d %d %d %d %d %f %d %d %.5f %d %f %f %f'%(destination,back,timestep,krecord,
																			rep,L,s[k],N0[j],K[j],mu[i],gen_num,cost,r,N1r[l])
				if file2run[-1] == "y": # file is .py
					cd('python %s %s'%(file2run, params))
				else: # file is .c
					cd('gcc -Wall %s -o cfile -lm'%(file2run))
					cd('./cfile %s %d'%(params, seed))
				count += 1
				print("%d/%d DONE\n"%(count, len(N0)*len(N1r)*len(s)*len(mu)))
