import numpy as np
import matplotlib.pyplot as plt

from functional import *

def test_charmed_hadron_cross_section():
	s = 13e3**2
	M = 0.938
	E = (s - 2 * M**2) / (2 * M)

	x = np.logspace(-10, 0, 100000)
	dx = x[1:] - x[:-1]

	had = ['D0', 'D+', 'D+s']
	val = [2.072, 0.834, 0.353]

	f = 0.0

	for h in zip(had, val):
		dsig = charmed_hadron_differential_production(x, E, h[0])[1:]
		sig = np.sum(dx * dsig)
		r = 2 * sig / h[1]
		f += r / len(had)
		print(f'{h[0]} : {r:.6f}')

	print(f'\navg : {f:.6f}')

test_charmed_hadron_cross_section()
