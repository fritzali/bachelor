'''Parametrizations of spectral distributions as described in the thesis document.'''

import numpy as np
from warnings import warn


def charmed_hadron_neutrinos(Enu, Eh, h):
	'''Return the charmed hadron to neutrino singular decay spectrum.

	Parameters
	----------
	Enu : float
		The energy of produced neutrinos in a proton rest frame view
	Eh : float
		The energy of decayed charmed hadrons in a proton rest frame view
	h : {'d0', 'd+', 'd+s', 'lam+c'}
		The type of hadronic initial state observed

	Returns
	-------
	float
		The charmed hadron to neutrino singular decay spectrum
	'''
	match h.lower():
		case 'd0':
			l = 0.67**2 / 1.86**2
		case 'd+':
			l = 0.63**2 / 1.87**2
		case 'd+s':
			l = 0.84**2 / 1.97**2
		case 'lam+c':
			l = 1.27**2 / 2.29**2
		case _:
			raise ValueError(f'`{h.lower()}` is an invalid charmed hadron identifyer, use `d0`, `d+`, `d+s` or `lam+c` instead')
	y = Enu / Eh
	if y < 0 or y > 1 - l:
		warn(f'{y} is outside of bounds {0.} to {1 - l}')
	a = 1 - l
	b = 1 - 2 * l
	D = 1 - 8 * l - 12 * l**2 * np.log(l) + 8 * l**3 - l**4
	F = (6 * b * a**2 - 4 * a**3 - 12 * l**3 * a + 12 * l**2 * y - 6 * b * y**2 + 4 * y**3 + 12 * l**2 * np.log((1 - y) / l)) / D
	return F / Eh







import matplotlib.pyplot as plt

x = np.linspace(0, 1, 1000)
plt.plot(x, np.vectorize(charmed_hadron_neutrinos)(x, 1, 'd0'))
#plt.plot(x, charmed_hadron_neutrinos(x, 1, 'd+'))
#plt.plot(x, charmed_hadron_neutrinos(x, 1, 'd+s'))
#plt.plot(x, charmed_hadron_neutrinos(x, 1, 'lam+c'))
plt.show()
