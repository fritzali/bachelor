import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

from functional import *

def integral(i, a, b, c, d):
	return quad(i, a, b, args = (c, d))

def test_charmed_hadron_production():
	return quad(charmed_hadron_production, 0, 1, args = (1e8, 'd0'))
#	s = 13e3**2
#	M = 0.938
#	m = 1.86
#	E = (s - m**2 - M**2 / (2 * M))
#	
#	x = np.linspace(8e-8, 8e-6, 100)
#
#	return quad(charm_quark_differential_production, 1e-8, 1, args = (E))
#	y = np.array(np.vectorize(integral)(charmed_hadron_differential_production, x, 1, 1e4, 'd0'))[0, :]
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


print(test_charmed_hadron_production())
