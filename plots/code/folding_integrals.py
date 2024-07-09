'''
Calculation of folding integrals to convert between spectra.

	Functions
	---------

'''
import numpy as np

import distribution_spectra as ds

def hadron_spectrum(Eh, h, lim = 1e15, N = 1000, NN = 1000):
	Ep = np.linspace(Eh, lim, N)
	dEp = Ep[1] - Ep[0]
	x = Eh / Ep
	match h.lower():
		case 'pi':
			prod = np.vectorize(ds.meson_production)(x, Ep, h)
		case 'k':
			prod = np.vectorize(ds.meson_production)(x, Ep, h)
		case 'd0':
			prod = np.vectorize(ds.charmed_hadron_production)(x, Ep, h, NN)
		case 'd+':
			prod = np.vectorize(ds.charmed_hadron_production)(x, Ep, h, NN)
		case 'd+s':
			prod = np.vectorize(ds.charmed_hadron_production)(x, Ep, h, NN)
		case 'lam+c':
			prod = np.vectorize(ds.charmed_hadron_production)(x, Ep, h, NN)
		case _:
			raise ValueError(f'`{h.lower()}` is an invalid hadron identifyer, use `pi`, `k`, `d0`, `d+`, `d+s` or `lam+c` instead')
	dspec = (prod / Ep**2)[1:]
	return np.sum(dEp * dspec)
