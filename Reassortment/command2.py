# was used for making back_test2 data.
from os import system as cd
import sys
import numpy as np
seed = np.random.randint(-9223372036854775808,-1)
back = 0
timestep = 1
krecord = 0 #  2= smallest k of the subpop. 1= all indiv's k. 0=mean k.
untilext = 0
rep = 100
L = 300
s = 0.05 # [0,0.05,0.1,0.15,0.2,0.25,0.3]
N0 = 1000
K = 1000
mu = 0.0007 #[0.0002,0.0003,0.0004,0.0005,0.0006,0.0007,0.0008,0.0009,0.0010,0.0011,0.0012]
gen_num = 500
cost = 0.0 #[0.0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16]
r = 0.5
N1r = [0,1] #[0.6,0.7,0.8,0.9]
destination = 'test'
file2run = sys.argv[1]
version = file2run[4:7]

count = 0
#for j in range(len(cost)):
#	for k in range(len(N1r)):
#		for i in range(len(s)):
for i in range(len(N1r)):
	if version == '1.2':
		params = '%s %d %d %d %d %d %d %f %d %d %.5f %d %f %f %f'%(destination,back,timestep,krecord,untilext,
																	rep,L,s,N0,K,mu,gen_num,cost,r,N1r[i])
	elif version == '1.3':
		params = '%s %d %d %d %d %d %f %d %d %.5f %d %f %f %f'%(destination,backtimestep,krecord,
																rep,L,s,N0,K,mu,gen_num,cost,r,N1r[i])
	if file2run[-1] == "y": # file is .py
		cd('python %s %s'%(file2run, params))
	else: # file is .c
		cd('gcc -Wall %s -o cfile -lm'%(file2run))
		cd('./cfile %s %d'%(params, seed))
	count += 1
	print("%d/%d DONE\n"%(count, len(N1r)))
