import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as const

arr1 = np.logspace(0, 1, 1000)
arr2 = np.logspace(0.2, 1.2, 1000)
arr3 = np.exp(arr1) / 1000
arr4 = np.sqrt(arr2)

plt.plot(arr1, label='log')
plt.plot(arr2, label='log')
plt.plot(arr3, label='exp')
plt.plot(arr4, label='sqrt')

plt.xlabel(r'$t$ / s')
plt.ylabel(r'$E_\nu^2\dot{\phi}_\nu$ / G$\kern0.5pt$e$\kern-0.5pt$V$\kern1pt$s$^{-1}$cm$^{-2}$')

plt.legend()

plt.savefig('build/test.pdf')
plt.savefig('build/test.jpg')
plt.close()
