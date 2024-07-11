import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

from functional import *

def integral(i, a, b, c, d):
	return quad(i, a, b, args = (c, d))

def test_charmed_hadron_production():
	s = 13e3**2
	M = 0.938
	E = (s - 2 * M**2) / (2 * M)

	x = np.logspace(-10, 0, 100000)
	dx = x[1:] - x[:-1]

	had = ['d0', 'd+', 'd+s']
	val = [2.072, 0.834, 0.353]

	f = 0.0

	for h in zip(had, val):
		dsig = charmed_hadron_differential_production(x, E, h[0])[1:]
		sig = np.sum(dx * dsig)
		r = 2 * sig / h[1]
		f += r / len(had)
		print(f'{h[0]} : {r:.6f}')

	print(f'\navg : {f:.6f}')


#	x = np.linspace(8e-8, 8e-6, 100)

#	return quad(charmed_hadron_differential_production, 0, 1, args = (E, 'd0'))

	

#	return quad(charm_quark_differential_production, 1e-8, 1, args = (E))
#	
#	plt.plot(x, y)
#	plt.show()
#
#	return 0
#
#	x = np.linspace(0, 1, 100)
#	dx = x[1] - x[0]
#	dsig = charmed_hadron_differential_production(x, E, 'd0')[1:]
#	sig = np.sum(dx * dsig)
#	return sig


#print(
test_charmed_hadron_production()
#)
