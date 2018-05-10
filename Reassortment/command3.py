# was used for making ratio_test, cost_test, and N0ratio_test, r_test data.

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
s = 0.05
N0 = 1000
K = 1000
mu = 0.00080 #[0.0002,0.0003,0.0004,0.0005,0.0006,0.0007,0.00080,0.0009,0.0010,0.0011,0.0012]
gen_num = 500
cost = [0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95]
r = 0.5
N1r = 0.5 #[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
destination = 'cost_test'
file2run = sys.argv[1]
version = file2run[4:7]
count = 0

for i in range(len(cost)):
	if version == '1.2':
		params = '%s %d %d %d %d %d %d %f %d %d %.5f %d %f %f %f'%(destination,back,timestep,krecord,untilext,
																	rep,L,s,N0,K,mu,gen_num,cost[i],r,N1r)
	elif version == '1.3':
		params = '%s %d %d %d %d %d %f %d %d %.5f %d %f %f %f'%(destination,back,timestep,krecord,
																rep,L,s,N0,K,mu,gen_num,cost[i],r,N1r)
	if file2run[-1] == "y": # file is .py
		cd('python %s %s'%(file2run, params))
	else: # file is .c
		cd('gcc -Wall %s -o cfile -lm'%(file2run))
		cd('./cfile %s %d'%(params, seed))
	count += 1
	print("%d/%d DONE\n"%(count, len(cost)))
count = 0
cost = 0.0
N1r = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
N0 = [1000,10000,100000,200000]
K = [1000,10000,100000,200000]
destination = 'N0ratio_test'

for i in range(len(N0)):
	for j in range(len(N1r)):
		if version == '1.2':
			params = '%s %d %d %d %d %d %d %f %d %d %.5f %d %f %f %f'%(destination,back,timestep,krecord,untilext,
																		rep,L,s,N0[i],K[i],mu,gen_num,cost,r,N1r[j])
		elif version == '1.3':
			params = '%s %d %d %d %d %d %f %d %d %.5f %d %f %f %f'%(destination,back,timestep,krecord,
																	rep,L,s,N0[i],K[i],mu,gen_num,cost,r,N1r[j])
		if file2run[-1] == "y": # file is .py
			cd('python %s %s'%(file2run, params))
		else: # file is .c
			cd('gcc -Wall %s -o cfile -lm'%(file2run))
			cd('./cfile %s %d'%(params, seed))
		count += 1
		print("%d/%d DONE\n"%(count, len(N0)*len(N1r)))
count = 0
N1r = 0.5
N0 = 1000
K = 1000
r = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
destination = 'r_test'

for i in range(len(r)):
	if version == '1.2':
		params = '%s %d %d %d %d %d %d %f %d %d %.5f %d %f %f %f'%(destination,back,timestep,krecord,untilext,
																	rep,L,s,N0,K,mu,gen_num,cost,r[i],N1r)
	elif version == '1.3':
		params = '%s %d %d %d %d %d %f %d %d %.5f %d %f %f %f'%(destination,back,timestep,krecord,
																rep,L,s,N0,K,mu,gen_num,cost,r[i],N1r)
	if file2run[-1] == "y": # file is .py
		cd('python %s %s'%(file2run, params))
	else: # file is .c
		cd('gcc -Wall %s -o cfile -lm'%(file2run))
		cd('./cfile %s %d'%(params, seed))
	count += 1
	print("%d/%d DONE\n"%(count, len(r)))

count = 0
r = 0.5
s = [0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95]
destination = 's_test'
for i in range(len(s)):
	if version == '1.2':
		params = '%s %d %d %d %d %d %d %f %d %d %.5f %d %f %f %f'%(destination,back,timestep,krecord,untilext,
																	rep,L,s[i],N0,K,mu,gen_num,cost,r,N1r)
	elif version == '1.3':
		params = '%s %d %d %d %d %d %f %d %d %.5f %d %f %f %f'%(destination,back,timestep,krecord,
																rep,L,s[i],N0,K,mu,gen_num,cost,r,N1r)
	if file2run[-1] == "y": # file is .py
		cd('python %s %s'%(file2run, params))
	else: # file is .c
		cd('gcc -Wall %s -o cfile -lm'%(file2run))
		cd('./cfile %s %d'%(params, seed))
	count += 1
	print("%d/%d DONE\n"%(count, len(s)))

count = 0
s = 0.05
cost = [0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95]
N1r = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
destination = 'cost_test'
for i in range(len(cost)):
	for j in range(len(N1r)):
		if version == '1.2':
			params = '%s %d %d %d %d %d %d %f %d %d %.5f %d %f %f %f'%(destination,back,timestep,krecord,untilext,
																		rep,L,s,N0,K,mu,gen_num,cost[i],r,N1r[j])
		elif version == '1.3':
			params = '%s %d %d %d %d %d %f %d %d %.5f %d %f %f %f'%(destination,back,timestep,krecord,
																	rep,L,s,N0,K,mu,gen_num,cost[i],r,N1r[j])
		if file2run[-1] == "y": # file is .py
			cd('python %s %s'%(file2run, params))
		else: # file is .c
			cd('gcc -Wall %s -o cfile -lm'%(file2run))
			cd('./cfile %s %d'%(params, seed))
			count += 1
		print("%d/%d DONE\n"%(count, len(cost)*len(N1r)))


