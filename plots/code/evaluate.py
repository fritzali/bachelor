import numpy as np
import matplotlib.pyplot as plt
import time

from functional import *

def test_charmed_hadron_cross_section():
	start = time.perf_counter()

	s = 13e3**2
	M = 0.938
	E = (s - 2 * M**2) / (2 * M)

	x = np.logspace(-10, 0, 100000)
	dx = x[1:] - x[:-1]

	had = ['D0', 'D+', 'D+s']
	val = [2.072, 0.834, 0.353]

	f = np.empty(0)

	for h in zip(had, val):
		dsig = charmed_hadron_differential_production(x, E, h[0])[1:]
		sig = np.sum(dx * dsig)
		r = 2 * sig / h[1]
		f = np.append(f, r)

	print(f'\nenergy: {E / 1e6:.0f} PeV\n')

	print('hadron:\t\tratio:\t\tvalue:')
	for h in zip(had, f, val):
		print(f'{h[0]}\t\t{h[1]:.6f}\t{h[2]:1.3f} mb')

	print(f'\naverage: {np.mean(f):.6f} +- {np.std(f):.6f}')

	end = time.perf_counter()

	print(f'\nelapsed: {end - start:.3f} s\n')

test_charmed_hadron_cross_section()
