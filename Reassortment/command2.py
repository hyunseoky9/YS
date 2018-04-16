# was used for making back_test2 data.

from os import system as cd
import sys
import numpy as np
back = [0,1]
timestep = 1
krecord = 0 #  2= smallest k of the subpop. 1= all indiv's k. 0=mean k.
untilext = 0
rep = 100
L = 300
s = 0.05
N0 = 1000
K = 1000
mu = [0.0002,0.0003,0.0004,0.0005,0.0006,0.0007,0.0008,0.0009,0.0010,0.0011,0.0012]
gen_num = 500
cost = 0
r = 0.5
N1r = [0,1]
destination = 'back_test2_1.3c'
file2run = sys.argv[1]
version = file2run[4:7]

count = 0
for j in range(len(back)):
	for k in range(len(N1r)):
		for i in range(len(mu)):

			if version == '1.2':
				params = '%s %d %d %d %d %d %d %f %d %d %.5f %d %f %f %f'%(destination,back[j],timestep,krecord,untilext,
																			rep,L,s,N0,K,mu[i],gen_num,cost,r,N1r[k])
			elif version == '1.3':
				params = '%s %d %d %d %d %d %f %d %d %.5f %d %f %f %f'%(destination,back[j],timestep,krecord,
																		rep,L,s,N0,K,mu[i],gen_num,cost,r,N1r[k])
			if file2run[-1] == "y": # file is .py
				cd('python %s %s'%(file2run, params))
			else: # file is .c
				cd('gcc -Wall %s -o cfile -lm'%(file2run))
				cd('./cfile %s'%(params))
			count += 1
			print("%d/%d DONE"%(count, len(mu)*len(back)*len(N1r)))

#elif version == '1.3':
#	params = '%s %d %d %d %d %f %d %d %.5f %d %f %f %f'%(destination,back,timestep,krecord,
#															rep,L,s,N0,K,mu[i],gen_num,cost,r,N1r)

#for j in range(len(mu)):
#	for i in range(len(N1r)):
#		params = '%d %d %f %d %d %.5f %d %f %f %f'%(rep,L,s,N0,K,mu[j],gen_num,cost,r,N1r[i])
#		file2run = sys.argv[1]
#		cd('python %s %s'%(file2run, params))
