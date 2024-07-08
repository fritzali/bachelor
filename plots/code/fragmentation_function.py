'''
Parametrization of a fragmentation function as described in the thesis document.

	Functions
	---------
	charmed_hadron_fragmentation_function
		Returns the charmed hadrons from charm quarks fragmentation function

'''
import numpy as np


def charmed_hadron_fragmentation_function(z, h):
	'''
	Returns the charmed hadron from charm quarks fragmentation function.

		Parameters
		----------
		z : float
			The energy ratio Eh / Ec of resulting charmed hadron to produced charm quark in target rest coordinates
		h : {'d0', 'd+', 'd+s', 'lam+c'}
			The type of hadronic final state observed

		Returns
		-------
		float
			The dimensionless charmed hadron from charm quarks fragmentation function
	'''
	match h.lower():
		case 'd0':
			N   = 0.577
			eps = 0.101
		case 'd+':
			N   = 0.238
			eps = 0.104
		case 'd+s':
			N   = 0.0327
			eps = 0.0322
		case 'lam+c':
			N   = 0.0067
			eps = 0.00418
		case _:
			raise ValueError(f'`{h.lower()}` is an invalid charmed hadron identifyer, use `d0`, `d+`, `d+s` or `lam+c` instead')
	return N * z * (1 - z)**2 / ((1 - z)**2 + eps * z)**2
