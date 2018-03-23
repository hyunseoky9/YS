Todo:
check if 1.2 follow f0 eq when there's no mutation for 1 segment and 2 segment
- find a parameter setting where 2seg is winning
	- try mutation from 0 to 0.0009 with 100reps

- change starting ratio, when the params are set for 2seg winning.







problem with parallel:
somehow when running the generations, it's not running the same way as it was sequential....
Gotta check where its wrong starting from the top of the script.

-Let's see if it's different with 1.2 ver

File nomenclature:
c1.2s_(params) -> c = comp, 1.2 = ver, s = stat
param order: rep,L,s,N0,K,mu,gen_num,cost,r

Questions:
- should mutation reflect fitness of the current generation or only of the next generation?

Comments/findings:

- as ratio of 1seg:2seg gets to around 7:3, 1segments start to win with mutation rate 0.005. 