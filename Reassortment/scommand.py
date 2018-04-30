# was used for making step data

from os import system as cd
import sys
import numpy as np
seed = np.random.randint(-9223372036854775808,-1)
krecord = 0 #  2= smallest k of the subpop. 1= all indiv's k. 0=mean k.
rep = 10000
s = 0.05
N0 = 1000
initkmu = [2,3,4,5,6,7,8,9,10,11,12];
cost = 0.0
r = 0.5
N1orN2 = [1,2]
destination = 'step'
file2run = sys.argv[1]
version = file2run[4:7]
count = 0
for i in range(len(initkmu)):
	for j in range(len(N1orN2)):
		params = '%s %d %d %.2f %d %.3f %.2f %.2f %d'%(destination,krecord,rep,s,N0,initkmu[i],cost,r,N1orN2[j])
		cd('gcc -Wall %s -o cfile -lm'%(file2run))
		cd('./cfile %s %d'%(params, seed))
		count += 1
		print("%d/%d DONE\n"%(count, len(initkmu)*len(N1orN2)))
