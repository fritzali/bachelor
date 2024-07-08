'''Parametrization of decay and collision as described in the thesis document.'''

import numpy as np
import scipy.constants as con

import cross_sections as cr


def hadron_proton_cooling_factor(E, h, d):
	'''Return the cooling factor for hadrons scattered by protons

	Parameters
	----------
	E : float
		The energy Eh as viewed in resting proton coordinates
	h : {'pi', 'k', 'd0', 'd+', 'd+s', 'lam+c'}
		The incident hadron type

	Returns
	-------
	float
		The charmed hadron `h` from charm quark `c` fragmentation function
	'''
	match h.lower():
		case 'pi':
			m   = 
			tau = 
			s   = 
			sig = cr.inelastic_hadron_proton_scattering(s, 'pi')
		case 'k':
			m   = 
			tau = 
			s   = 
			sig = cr.inelastic_hadron_proton_scattering(s, 'k')
		case 'd0':
			m   = 
			tau = 
			s   = 
			sig = cr.inelastic_hadron_proton_scattering(s, 'k')
		case 'd+':
			m   = 
			tau = 
			s   = 
			sig = cr.inelastic_hadron_proton_scattering(s, 'k')
		case 'd+s':
			m   = 
			tau = 
			s   = 
			sig = cr.inelastic_hadron_proton_scattering(s, 'k')
		case 'lam+c':
			m   = 
			tau = 
			s   = 
			sig = cr.inelastic_hadron_proton_scattering(s, 'k')
		case _:
			raise ValueError(f'`{h.lower()}` is an invalid hadron identifyer, use `pi`, `k`, `d0`, `d+`, `d+s` or `lam+c` instead')
	return N * z * (1 - z)**2 / ((1 - z)**2 + eps * z)**2


def proton_proton_optical_depth(E):
