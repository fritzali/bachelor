'''
Parametrizations of spectral distributions as described in the thesis document.

	Functions
	---------
	meson_production
		Returns the proton-proton to pion or kaon singular production spectrum

	meson_decay_neutrinos
		Returns the pion or kaon to neutrino singular decay spectrum

	charmed_hadron_decay_neutrinos
		Returns the charmed hadron to neutrino singular decay spectrum

'''
import numpy as np
from warnings import warn


def meson_production(x, E, h):
	'''
	Returns the proton-proton to pion or kaon singular production spectrum.

		Parameters
		----------
		x : float
			The energy ratio Eh / Ep of produced meson to incident proton in proton target rest coordinates
		E : float
			The projectile energy Ep from proton target rest coordinates in GeV
		h : {'pi', 'k'}
			The type of charged meson produced

		Returns
		-------
		float
			The proton-proton to pion or kaon singular production spectrum in 1 / GeV
	'''
	match h.lower():
		case 'pi':
			m = 0.140
			f = 1
		case 'k':
			m = 0.494
			f = 0.12
		case _:
			raise ValueError(f'`{h.lower()}` is an invalid hadron identifyer, use `pi` or `k` instead')
	B0 = 0.25
	a0 = 0.98
	r0 = 2.6
	c1 = 1.515
	c2 = 0.206
	c3 = 0.075
	C = c1 - c2 * np.log(E) + c3 * np.log(E)**2
	B = B0 + C
	a = a0 / np.sqrt(C)
	r = r0 / np.sqrt(C)
	u = (1 - m / (x * E))**(1/2)
	v = 1 - x**a
	w = 1 + r * x**a * v
	F = 4 * a * B * x**(a - 1) * (v / w)**4 * (1 / v + r * (1 - 2 * x**a) / w) * u
	return f * F / E


def meson_decay_neutrinos(Enu, Eh, h):
	'''
	Returns the pion or kaon to neutrino singular decay spectrum.

		Parameters
		----------
		Enu : float
			The energy of produced neutrinos from a proton rest frame view in GeV
		Eh : float
			The energy of decayed mesons from a proton rest frame view in GeV
		h : {'pi', 'k'}
			The type of meson initital state observed

		Returns
		-------
		float
			The pion or kaon to neutrino singular decay spectrum in 1 / GeV
	'''
	match h.lower():
		case 'pi':
			m = 0.140
			f = 0.9999
		case 'k':
			m = 0.494
			f = 0.6356
		case _:
			raise ValueError(f'`{h.lower()}` is an invalid hadron identifyer, use `pi` or `k` instead')
	l = 0.106**2 / m**2
	y = Enu / Eh
	if y > 1 - l:
		warn(f'{y} exceeds bound {1 - l}')
	return f / (Eh * (1 - l))


def charmed_hadron_decay_neutrinos(Enu, Eh, h):
	'''
	Returns the charmed hadron to neutrino singular decay spectrum.

		Parameters
		----------
		Enu : float
			The energy of produced neutrinos from a proton rest frame view in GeV
		Eh : float
			The energy of decayed charmed hadrons from a proton rest frame view in GeV
		h : {'d0', 'd+', 'd+s', 'lam+c'}
			The type of hadronic initial state observed

		Returns
		-------
		float
			The charmed hadron to neutrino singular decay spectrum in 1 / GeV
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
	if y > 1 - l:
		warn(f'{y} exceeds bound {0.} to {1 - l}')
	a = 1 - l
	b = 1 - 2 * l
	D = 1 - 8 * l - 12 * l**2 * np.log(l) + 8 * l**3 - l**4
	F = (6 * b * a**2 - 4 * a**3 - 12 * l**3 * a + 12 * l**2 * y - 6 * b * y**2 + 4 * y**3 + 12 * l**2 * np.log((1 - y) / l)) / D
	return F / Eh
