'''
Parametrizations of spectral distributions as described in the thesis document.

	Functions
	---------
	meson_production
		Returns the proton-proton to pion or kaon singular production spectrum

	meson_decay_neutrinos
		Returns the pion or kaon to neutrino singular decay spectrum

	charmed_hadron_production
		Returns the proton-proton to charmed hadron singular production spectrum

	charmed_hadron_decay_neutrinos
		Returns the charmed hadron to neutrino singular decay spectrum

'''
import numpy as np
from warnings import warn

import cross_sections as cr
import fragmentation_function as ff


def meson_production(x, E, h):
	'''
	Returns the proton-proton to pion or kaon singular production spectrum.

		Parameters
		----------
		x : float
			The energy ratio Eh / Ep of produced meson to incident proton in target rest coordinates
		E : float
			The projectile energy Ep as viewed from target rest coordinates in GeV
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
	if x < 0 or x > 1:
		warn(f'`{x}` is outside of bounds {0.0} and {1.0}')
		return 0.0
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
			The energy of produced neutrinos as viewed from target rest coordinates in GeV
		Eh : float
			The energy of decayed mesons as viewed from target rest coordinates in GeV
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
		return 0.0
	return f / (Eh * (1 - l))


def charmed_hadron_production(x, E, h, N = 100):
	'''
	Returns the proton-proton to charmed hadron singular production spectrum.

		Parameters
		----------
		x : float
			The energy ratio Eh / Ep of produced hadron to incident proton in target rest coordinates
		E : float
			The projectile energy Ep as viewed from target rest coordinates in GeV
		h : {'d0', 'd+', 'd+s', 'lam+c'}
			The type of hadron produced
		N : int, optional
			The number of steps for integration accuracy

		Returns
		-------
		float
			The proton-proton to charmed hadron singular production spectrum in 1 / GeV
	'''
	if x < 0 or x > 1:
		warn(f'`{x}` is outside of bounds {0.0} to {1.0}')
		return 0.0
	z = np.linspace(x, 1, N)
	dz = z[1] - z[0]
	dsig = (cr.charm_quark_differential_production(x / z, E) * ff.charmed_hadron_fragmentation_function(z, h) / z)[1:]
	sig = np.sum(dz * dsig)
	M = 0.938
	s = 2 * (E * M + M**2)
	return sig / (E * cr.inelastic_hadron_proton_scattering(s, 'p'))


def charmed_hadron_decay_neutrinos(Enu, Eh, h):
	'''
	Returns the charmed hadron to neutrino singular decay spectrum.

		Parameters
		----------
		Enu : float
			The energy of produced neutrinos as viewed from target rest coordinates in GeV
		Eh : float
			The energy of decayed charmed hadrons as viewed from target rest coordinates in GeV
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
			f = 0.067
		case 'd+':
			l = 0.63**2 / 1.87**2
			f = 0.176
		case 'd+s':
			l = 0.84**2 / 1.97**2
			f = 0.065
		case 'lam+c':
			l = 1.27**2 / 2.29**2
			f = 0.045
		case _:
			raise ValueError(f'`{h.lower()}` is an invalid charmed hadron identifyer, use `d0`, `d+`, `d+s` or `lam+c` instead')
	y = Enu / Eh
	if y > 1 - l:
		warn(f'{y} exceeds bound {1 - l}')
		return 0.0
	a = 1 - l
	b = 1 - 2 * l
	D = 1 - 8 * l - 12 * l**2 * np.log(l) + 8 * l**3 - l**4
	F = (6 * b * a**2 - 4 * a**3 - 12 * l**2 * a + 12 * l**2 * y - 6 * b * y**2 + 4 * y**3 + 12 * l**2 * np.log((1 - y) / l)) / D
	return f * F / Eh
