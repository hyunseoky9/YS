How to navigate through this project:

1. running simulations
2. data storage
3. data analysis


1. running simulations 
	There's 2 types of models: comp and meta
	1) comp model
		a) explanation
		comp model is a single population model where you can simulate 
		1segment and 2 segment to compete with each other. They don't have to
		compete if you set the parameter 'N1r' to 1 or 0.
		comp model has 3 different versions.
		ver 1 is a prototype, and practically useless. I first made the model with python so it has .py and .c (comp1.1)
		ver 2 has flexible population size with carrying capacity K (comp1.2). The latest one 1.2.3 is only logically sound and work fine with commandx.py 
		ver 3 is a WF model with constant population (comp1.3). 1.3 will only work perfectly. 1.3.2 is garbage now and 1.3.3 is iffy.

		comp2.3 and comp3.3 adds 3 segments and 8 segments respectively but needs a lot more work to run fine.

		b) how to run
		same as running meta model.
		To run simply type: 
		```python command.py comp(version).c (code)``` 

	2) meta model
		a) explanation
		meta model is a multi metapopulation model where you can simulate 1segments and 2 segments growing from multiple hosts.
		Each step, the viruses can migrate to the migration pool and transmit to other hosts with some transmission rate. 
		it's an expansion of the comp model.
		ver 1 and 2 allow viruses to catch multiple deleterious mutations, but ver 2 uses matrix computation. (meta1.1 & meta1.1.2)
		ver 3 allows a virus to only accumulates 1 mutation per time in order to save computation time. (meta1.1.3)

		** 1.1.2 and 1.1 doesn't yet allow multihost system as it was still under development. It cannot use some parameters: pop2init, pop2initlen, pop1init, pop1initlen, tr, mig.

		b) how to run
		meta model has more easy and convenient method of running the simulation.
		Instead of making multiple commands for different parameter sets in command*.py, here we use a csv file called 'paraminfo.csv' to store all the 
		different parameter combinations. Each row contains a combination of parameters that is used for a test that is described on the first column, 'description'. 
		To run the set of parameters on a row. You specify the code on a 'code' column for the row on shell command.
		When making some parameters, theres are some rules you need to follow:
			-set pop2initlen and pop1initlen to 0. Its just used for parsing in the c file.
			-pop2init and pop1init has to have a number between 0 and 1 followed by '~'sign. Each (number)~ pair is whatever segment's initial frequency at a host. 
			-Number of (number)~ pair in pop1init and pop2init has to match the host_num parameter.
			-if you want random seed, type 'random'. If you need a consistent random number, put any negative integer below -9223372036854775808.
		To run, simply type:
		```python command.py meta(ver).c (code)```

2. data storage
	All data is stored in ./data.
	Each simulation data is stored in a file. Each file has exp.txt that explains what the simulation was about.
	./data/test is a scrap and can be erased at your convenience.
3. data analysis.
	In ./data, there's 2 .ipynb jupyter notebook file named data_viewer and data_viewer2. 
	data_viewer has analysis of data from comp model and data_viewer2 has that of meta model.
	data_viewer model is more messy and you'll likely have to finegle around to get it working.
	data_viewer2 is more refined and you have to specify some parameters in the cell and the wanted result will show up.
	Follow the instructions in the notebook.






















