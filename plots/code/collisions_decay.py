'''
Parametrization of decay and collision as described in the thesis document.

	Functions
	---------
	hadron_proton_cooling_factor
	proton_proton_optical_depth

'''
import numpy as np
import scipy.constants as con

import cross_sections as cr


def hadron_proton_cooling_factor(E, n, h, d = None):
	'''
	Returns the cooling factor for hadrons scattered by protons.

		Parameters
		----------
		E : float
			The energy Eh in GeV as viewed from resting proton coordinates
		n : float
			The nucleon number density per cm cubed
		h : {'pi', 'k', 'd0', 'd+', 'd+s', 'lam+c'}
			The incident hadron type
		d : float, optional
			The target field size in cm, assumed to be infinite if `None`

		Returns
		-------
		float
			The cooling factor for hadrons scattered by protons
	'''
	M = 0.938
	match h.lower():
		case 'pi':
			m   = 0.140
			tau = 26.03e-9
			s   = 2 * E * M + M**2 + m**2
			sig = cr.inelastic_hadron_proton_scattering(s, 'pi') * 1e-24
		case 'k':
			m   = 0.494
			tau = 12.38e-9
			s   = 2 * E * M + M**2 + m**2
			sig = cr.inelastic_hadron_proton_scattering(s, 'k') * 1e-24
		case 'd0':
			m   = 1.86
			tau = 410e-15
			s   = 2 * E * M + M**2 + m**2
			sig = cr.inelastic_hadron_proton_scattering(s, 'k') * 1e-24
		case 'd+':
			m   = 1.87
			tau = 1033e-15
			s   = 2 * E * M + M**2 + m**2
			sig = cr.inelastic_hadron_proton_scattering(s, 'k') * 1e-24
		case 'd+s':
			m   = 1.97
			tau = 501e-15
			s   = 2 * E * M + M**2 + m**2
			sig = cr.inelastic_hadron_proton_scattering(s, 'k') * 1e-24
		case 'lam+c':
			m   = 2.29
			tau = 203e-15
			s   = 2 * E * M + M**2 + m**2
			sig = cr.inelastic_hadron_proton_scattering(s, 'k') * 1e-24
		case _:
			raise ValueError(f'`{h.lower()}` is an invalid hadron identifyer, use `pi`, `k`, `d0`, `d+`, `d+s` or `lam+c` instead')
	kap = 0.8
	c = con.c * 100
	t_cool = 1 / (kap * sig * n * c)
	t_dec = tau * E / m
	d_dec = t_dec * c
	if d is not None and d < d_dec:
		d_free = t_cool * c
		return 1 - np.exp(- d_free / d)
	return 1 - np.exp(- t_cool / t_dec)


def proton_proton_optical_depth(E, n, d):
	'''
	Returns the effective optical depth for protons hitting protons.

		Parameters
		----------
		E : float
			The energy Ep in GeV as viewed from resting proton coordinates
		n : float
			The nucleon number density per cm cubed
		d : float
			The target field size in cm

		Returns
		-------
		float
			The effective optical depth for protons hitting protons
	'''
	M = 0.938
	s = 2 * (E * M + M**2)
	kap = 0.5
	sig = cr.inelastic_hadron_proton_scattering(s, 'p') * 1e-24
	d_free = 1 / (kap * sig * n)
	return d / d_free
