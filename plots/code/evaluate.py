import numpy as np
import matplotlib.pyplot as plt
import time

from code.functional import *


def sample_inelastic_hadron_scattering():
	s = np.logspace(1, 10, 1000)

	p = inelastic_hadron_proton_scattering(s, 'p')
	pi = inelastic_hadron_proton_scattering(s, 'pi')
	K = inelastic_hadron_proton_scattering(s, 'k')

	with open('code/tabulate/other/sample_inelastic_scattering.txt', 'w') as f:
		f.write(f'# Hadron Inelastic Scattering Cross Section\n')
		f.write(f'# s / GeV**2 # sig (proton) / mb # sig (pion) / mb # sig (kaon) / mb\n')
		for row in zip(s, p, pi, K):
			f.write(r'{0}   {1}   {2}   {3}'.format(*row))
			f.write(f'\n')


def sample_charmed_hadron_cross_section():
	x = np.logspace(-7, 0, 1000)

	sig1 = charmed_hadron_differential_production(x, 1e12, 'd0')
	sig2 = charmed_hadron_differential_production(x, 1e10, 'd0')
	sig3 = charmed_hadron_differential_production(x, 1e8, 'd0')

	with open('code/tabulate/other/sample_charm_hadron.txt', 'w') as f:
		f.write(f'# `D0` Sample Charmed Hadron Cross Section\n')
		f.write(f'# x # sig (1e12 GeV) / mb # sig (1e10 GeV) / mb # sig (1e8 GeV) / mb\n')
		for row in zip(x, sig1, sig2, sig3):
			f.write(r'{0}   {1}   {2}   {3}'.format(*row))
			f.write(f'\n')


def test_charmed_hadron_cross_section():
	start = time.perf_counter()

	s = 13e3**2
	M = 0.938
	E = (s - 2 * M**2) / (2 * M)

	x = np.logspace(-10, 0, 10000)

	had = ['D0', 'D+', 'D+s']
	val = [2.072, 0.834, 0.353]

	f = np.empty(0)

	for h in zip(had, val):
		y = charmed_hadron_differential_production(x, E, h[0])
		sig = np.trapezoid(y, x)
		r = 2 * sig / h[1]
		f = np.append(f, r)

	print(f'\nenergy: {E / 1e6:.0f} PeV\n')

	print('hadron:\t\tratio:\t\tvalue:')
	for h in zip(had, f, val):
		print(f'{h[0]}\t\t{h[1]:.6f}\t{h[2]:1.3f} mb')

	print(f'\naverage: {np.mean(f):.6f} +- {np.std(f):.6f}')

	end = time.perf_counter()

	print(f'\nelapsed: {end - start:.3f} s\n')


sample_inelastic_hadron_scattering()
sample_charmed_hadron_cross_section()
test_charmed_hadron_cross_section()
