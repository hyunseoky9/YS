import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt

n = 100
adv = 0
u1 = 0.25
u2 = 0.75
d=0.25
sig = 0.125
ratios = []
for i in range(100):
	for j in range(n):
		a11 = np.random.random_sample()*2*d
		a = np.array([a11, (a11+d), (a11+2*d)])
		w1 = st.norm.pdf(a,u1,sig)
		w2 = st.norm.pdf(a,u2,sig)
		ww = w1*w2
		if ( ww[0] < ww[1] and ww[1]> ww[2]):
			adv += 1
	ratios.append(adv/n)
	adv = 0

ratios = np.array(ratios)
plt.hist(ratios)
plt.show()