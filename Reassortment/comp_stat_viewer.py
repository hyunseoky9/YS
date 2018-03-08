import numpy as np
import pandas
import matplotlib.pyplot as plt
import sys

file_name = sys.argv[1]
data = pandas.read_csv(file_name)

print('1seg won:',len(np.where(data['pop1']>0)[0]))
print('2seg won:',len(np.where(data['pop2']>0)[0]))

#plt.plot(data['pop1'], label='1seg')
#plt.plot(data['pop2'], label='2seg')
#plt.title('population ')
#plt.legend()
#plt.show()

#plt.plot(data['k1'], label='1seg')
#plt.plot(data['k2'], label='2seg')
#plt.legend()
#plt.title('mean_k')
#plt.show()
