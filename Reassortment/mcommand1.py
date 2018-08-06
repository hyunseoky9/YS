# was used to make ultimate_test and ultimate_test2 data.
# parameter adding needs to be fixed for comp1.2!!
from os import system as cd
import sys
import numpy as np
import time
from datetime import date

def recursion(list_index,data,data2,z,params):
	if len(list_index) == 0:
		string = ''
		for j in range(len(data2)):
			string = '%s %s'%(string,data2[j])
		params.append(string)
	else:
		if z < len(list_index) - 1:
			for i in range(len(data[list_index[z]])):
				data2[list_index[z]] = data[list_index[z]][i]
				params = recursion(list_index,data,data2,z+1,params)
		else:
			for i in range(len(data2[list_index[z]])):
				string = ''
				for j in range(len(data2)):
					if type(data2[j]) == type('string'):
						string = '%s %s'%(string,data2[j])
					if type(data2[j]) == type([]):
						string = '%s %s'%(string,data2[j][i])
				params.append(string)
				#print(string)
	return params

td = date.today() #simulation start date
start = time.time()
handle = open('paraminfo.csv')

code_match = 2 # select correct code to run
line_num = 0
list_index = []
dic = {}
for line in handle:
	if line_num:
		data = line.strip().split(",")
		data2 = line.strip().split(",")
		if data[1] == str(code_match):
			for i in range(len(data)):
				if len(data[i].split('*')) > 1:
					list_index.append(i - 2) # considering we're going to take out 'code' and 'description' from data array 
					data[i] = data[i].split('*')
					data2[i] = data2[i].split('*')
			break
	else:
		data = line.strip().split(",")
		for i in range(len(data)):
			dic[data[i]] = i - 2 # considering we're going to take out 'code' and 'description' from data array 
	line_num += 1
handle.close()
data = data[2::]
data2 = data2[2::]

#change seed
if data[dic['seed']] == 'random':
	seed = str(np.random.randint(-9223372036854775808,-1))
	data[dic['seed']] = seed
	data2[dic['seed']] = seed

# check if init info and host_num match
hostcount = 0
for i in range(len(data[dic['pop2init']])):
	if data[dic['pop2init']][i] == '~':
		hostcount += 1
if hostcount != int(data[dic['host_num']]):
	raise ValueError("pop2init info and host number doesn't match")

# set pop1init length and pop2 init length
data[dic['pop2initlen']] = str(len(data[dic['pop2init']]))
data[dic['pop1initlen']] = str(len(data[dic['pop1init']]))
data2[dic['pop2initlen']] = str(len(data2[dic['pop2init']]))
data2[dic['pop1initlen']] = str(len(data2[dic['pop1init']]))

# fill in all the parameter vectors to try
params = []
params = recursion(list_index,data,data2,0,params)
file2run = sys.argv[1]
version = file2run[4:7]

for i in range(len(params)):
	#print('./cfile%s'%(params[i]))
	cd('gcc -Wall %s -o cfile -lm'%(file2run))
	cd('./cfile%s'%(params[i]))
	print('%d/%d DONE'%((i+1),len(params)))

end = time.time()
print('it took following minutes: ', (end - start)/60)
