# was used to make ultimate_test and ultimate_test2 data.
# parameter adding needs to be fixed for comp1.2!!
from os import system as cd
import sys
import numpy as np
import time
from datetime import date
td = date.today() #simulation start date
start = time.time()
seed = -1 #np.random.randint(-9223372036854775808,-1)
timestep = 0
krecord = 0
untilext = 0
N0 = 1000
r = 0.5
host_num = 1
kmax = 20
rep = 10000
gen_num = 500
u = 0.24
s = 0.05
c = 0.00
K = 1000
destination = 'test'
file2run = sys.argv[1]
version = file2run[4:7]

if version == '1.1':
	params = '%s %d %d %d %d %f %d %d %.5f %d %f %f %d %d %d'%(destination,timestep,krecord,untilext,rep,s,N0,K,u,gen_num,c,r,seed,host_num,kmax)
if file2run[-1] == "y": # file is .py
	cd('python %s %s'%(file2run, params))
else: # file is .c
	cd('gcc -Wall %s -o cfile -lm'%(file2run))
	cd('./cfile %s'%(params))

end = time.time()
print('it took following minutes: ', (end - start)/60)

