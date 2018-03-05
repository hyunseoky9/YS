import numpy as np
import pandas
import matplotlib.pyplot as plt
import sys

file_name = sys.argv[1]
data = pandas.read_csv(file_name)
plt.plot(data['pop1'], label='1seg')
plt.plot(data['pop2'], label='2seg')
plt.title('population ')
plt.legend()
plt.show()

plt.plot(data['k1'], label='1seg')
plt.plot(data['k2'], label='2seg')
plt.legend()
plt.title('mean_k')
plt.show()
