import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as const

arr = np.logspace(0, 1, 1000)
plt.plot(arr)
plt.savefig('build/test.pdf')
plt.savefig('build/test.jpg')
plt.close()
