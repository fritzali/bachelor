'''Parametrizations of hadronic cross sections as described in the thesis document.'''

import numpy as np
from warnings import warn


def total_hadron_proton_scattering(s, h):
	'''Return the total hadron-proton scattering cross section.

	Parameters
	----------
	s : float
		The squared center of mass energy in GeV
	h : {'p', 'pi', 'k'}
		The incident hadron on the proton target

	Returns
	-------
	float
		The total hadron-proton scattering cross section for `h` at `s` in mb
	'''
	match h.lower():
		case 'p':
			P  = 34.41
			R1 = 13.07
			R2 =  7.39
			sh = 15.98
		case 'pi':
			P  = 18.75
			R1 =  9.56
			R2 =  1.767
			sh = 10.23
		case 'k':
			P  = 16.36
			R1 =  4.29
			R2 =  3.408
			sh = 12.62
		case _:
			raise ValueError(f'`{h.lower()}` is an invalid hadron identifyer, use `p`, `pi` or `k` instead')
	H  = 0.272
	n1 = 0.447
	n2 = 0.5486
	return H * np.log(s / sh)**2 + P + R1 * (sh / s)**n1 + R2 * (sh / s)**n2


def hadron_elastic_total_ratio(s):
	'''Return the universal ratio of elastic to total hadron-proton cross section.

	Parameters
	----------
	s : float
		The squared center of mass energy in GeV

	Returns
	-------
	float
		The universal ratio of elastic to total hadron-proton cross section
	'''
	A  = 1/2
	g1 = 0.466
	g2 = 0.0259
	g3 = 0.00177
	return A * np.tanh(g1 - g2 * np.log(s) + g3 * np.log(s)**2)


def inelastic_hadron_proton_scattering(s, h):
	'''Return the inelastic hadron-proton scattering cross section.

	Parameters
	----------
	s : float
		The squared center of mass energy in GeV
	h : {'p', 'pi', 'K'}
		The incident hadron on the proton target

	Returns
	-------
	float
		The inelastic hadron-proton scattering cross section for `h` at `s` in mb
	'''
	return total_hadron_proton_scattering(s, h) * (1 - hadron_elastic_total_ratio(s))


def charm_quark_differential_production(x, E):
	'''Return the charm quark differential cross section for production in proton-proton collisions.

	Parameters
	----------
	x : float
		The energy ratio Ec / Ep of charm quark to incident proton in proton target rest coordinates
	E : float
		The projectile energy Ep in GeV from proton target rest coordinates

	Returns
	-------
	float
		The charm quark differential cross section for production in proton-proton collisions in mb
	'''
	if E < 1e4 or E > 1e11:
		warn(f'{E} is outside of bounds {1e4} to {1e11}')
	if E >= 1e8:
		a1 = 0.403
		a2 = 2.002
		b1 = 0.237
		b2 = 0.023
		n1 = 7.639
		n2 = 0.102
	else:
		a1 = 0.826
		a2 = 8.411
		b1 = 0.197
		b2 = 0.016
		n1 = 1.061
		n2 = 0.107
	a = a1 * np.log(E) - a2
	b = b1 - b2 * np.log(E) - 1
	n = n1 + n2 * np.log(E)
	m = 1.2
	return a * x**b * (1 - x**m)**n / 14.5
