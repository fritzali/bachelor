'''
Parametrization of decay and collision as described in the thesis document.

	Functions
	---------
	hadron_proton_cooling_factor
		Returns the cooling factor for hadrons scattered by protons

	proton_proton_optical_depth
		Returns the effective optical depth for protons hitting protons

'''
import numpy as np

import code.parametrizations.cross_sections as cr


def hadron_proton_cooling_factor(E, n, h, d = None):
	'''
	Returns the cooling factor for hadrons scattered by protons.

		Parameters
		----------
		E : float
			The energy Eh as viewed from target rest coordinates in GeV
		n : float
			The nucleon number density in 1 / cm**3
		h : {'pi', 'k', 'd0', 'd+', 'd+s', 'lam+c'}
			The incident hadron type
		d : float, optional
			The target field size, assumed to be infinite if `None`, in cm

		Returns
		-------
		float
			The dimensionless cooling factor for hadrons scattered by protons
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
			tau = 0.410e-12
			s   = 2 * E * M + M**2 + m**2
			sig = cr.inelastic_hadron_proton_scattering(s, 'k') * 1e-24
		case 'd+':
			m   = 1.87
			tau = 1.033e-12
			s   = 2 * E * M + M**2 + m**2
			sig = cr.inelastic_hadron_proton_scattering(s, 'k') * 1e-24
		case 'd+s':
			m   = 1.97
			tau = 0.501e-12
			s   = 2 * E * M + M**2 + m**2
			sig = cr.inelastic_hadron_proton_scattering(s, 'k') * 1e-24
		case 'lam+c':
			m   = 2.29
			tau = 0.203e-12
			s   = 2 * E * M + M**2 + m**2
			sig = cr.inelastic_hadron_proton_scattering(s, 'k') * 1e-24
		case _:
			raise ValueError(f'`{h.lower()}` is an invalid hadron identifyer, use `pi`, `k`, `d0`, `d+`, `d+s` or `lam+c` instead')
	kap = 0.8
	c = 29979245800
	tcool = 1 / (kap * sig * n * c)
	tdec = tau * E / m
	ddec = tdec * c
	if d is not None and d < ddec:
		dfree = tcool * c
		return 1 - np.exp(- dfree / d)
	return 1 - np.exp(- tcool / tdec)


def proton_proton_optical_depth(E, n, d):
	'''
	Returns the effective optical depth for protons hitting protons.

		Parameters
		----------
		E : float
			The energy Ep as viewed from target rest coordinates in GeV
		n : float
			The nucleon number density in 1 / cm**3
		d : float
			The target field size in cm

		Returns
		-------
		float
			The dimensionless effective optical depth for protons hitting protons
	'''
	M = 0.938
	s = 2 * (E * M + M**2)
	kap = 0.5
	sig = cr.inelastic_hadron_proton_scattering(s, 'p') * 1e-24
	dfree = 1 / (kap * sig * n)
	return d / dfree
