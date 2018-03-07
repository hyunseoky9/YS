to make it efficient:
np.where rid -> use index array? only use k to do all the mutation work. Pick a random number and if it's less or equal to k back mutation, otherwise increase mutation

keep k_means during the generation simulation instead of going through all the k... (only for comp.py, rest is ok i think?)

terminate when either 1 subpop goes to 0

try graphing k_means

problem with parallel:
somehow when running the generations, it's not running the same way as it was sequential....
Gotta check where its wrong starting from the top of the script.

things to check:
np.random.poisson()


checked:
k_means in the beginnning
progeny_n nothing wrong
w_mean
