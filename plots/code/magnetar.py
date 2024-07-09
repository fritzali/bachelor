'''
Object oriented implementation of magnetar class as described in the thesis document.

	Classes
	-------
	magnetar
		Collects parameters and methods associated with the magnetar model

'''
import numpy as np
from warnings import warn

from parametrizations import *


class magnetar:
	'''
	Collects parameters and methods associated with the magnetar model.

	Attributes
	----------
	R : float
		The stellar radius in cm
	B : float
		The polar magnetic field strength in G
	o : float
		The initial angular frequency in rad / s
	chi : float
		The relative dipole tilt to rotational axis angle in rad
	I : float
		The moment of inertia in g * cm**2
	mu : float
		The magnetic moment in erg / G
	tsd : float
		The spindown time in s
	lum : float
		The initital luminosity in erg / s
	c : float
		The vacuum speed of light in cm / s
	e : float
		The elementary charge in esu

	Methods
	-------
	'''

	def __init__(self, R = 1e6, B = 1e15, o = 1e4, chi = 95e-2, I = 1e45, m = 'force free'):
		'''
		Constructs all attributes of the magnetar class.

		Parameters
		----------
		R : float, optional
			The stellar radius in cm
		B : float, optional
			The polar magnetic field strength in G
		o : float, optional
			The initial angular frequency in rad / s
		chi : float, optional
			The relative dipole tilt to rotational axis angle in rad
		I : float, optional
			The moment of inertia in g * cm**2
		m : {'force free', 'vacuum'}, optional
			The global magnetosphere model
		'''
		self.R = R
		self.B = B
		self.o = o
		self.chi = chi
		self.I = I
		mu = B * R**3 / 2
		c = 2.99792458e10
		e = 4.80320471e-10
		self.mu = mu
		self.c = c
		self.e = e
		match m.lower():
			case 'force free':
				K = mu**2 * (1 + np.sin(chi)**2) / c**3
			case 'vacuum':
				K = 2 * mu**2 * np.sin(chi)**2 / (3 * c**3)
			case _:
				raise ValueError(f'`{m.lower()}` is not a valid magnetosphere model, use `force free` or `vacuum` instead')
		self.tsd = I / (2 * K * o**2)
		self.lum = K * o**4

	def __str__(self):
		'''Defines string output for printing the magnetar object.'''
		str1 = f'Magnetar:\n    R = {self.R:.3} cm\n    B = {self.B:.3} G\n    o = {self.o:.3} rad / s\n    '
		str2 = f'chi = {self.chi:.3} rad\n    I = {self.I:.3} g * cm**2\n    mu = {self.mu:.3} erg / G\n    '
		str3 = f'tsd = {self.tsd:.3} s\n    lum = {self.lum:.3} erg / s\n    E = {self.E(0):.3} GeV\n    '
		str4 = f'spec = {self.proton_spectrum_prefactor(0):.3}\n    n = {self.number_density(self.tsd):.3} 1 / cm**3'
		return str1 + str2 + str3 + str4

	def L(self, t):
		'''
		Returns the spindown luminosity.

		Parameters
		----------
		t : float
			The time passed from magnetar formation in s

		Returns
		-------
			The spindown luminosity in erg / s
		'''
		return self.lum / (1 + t / self.tsd)**2

	def E(self, t, f = 1e-1):
		'''
		Returns the monotonic proton energy.

		Parameters
		----------
		t : float
			The time passed from magnetar formation in s
		f : float, optional
			The efficiency fraction of potential drop acceleration

		Returns
		-------
			The monotonic proton energy in GeV
		'''
		return f * self.e * self.B * self.R**3 * self.o**2 / (2 * self.c**2 * (1 + t / self.tsd)) * 624.150907

	def proton_spectrum_prefactor(self, t):
		'''
		Returns the prefactor of a delta functional proton spectrum.

		Parameters
		----------
		t : float
			The time passed from magnetar formation in s

		Returns
		-------
			The dimensionless prefactor of a delta functional proton spectrum
		'''
		return self.B * self.R**3 * self.o**2 / (self.c * self.e * (1 + t / self.tsd))

	def ejecta_radius(self, t, b = 1e-1):
		'''
		Returns the supernova ejecta radius.

		Parameters
		----------
		t : float
			The time passed from magnetar formation in s
		b : float, optional
			The relativistic velocity fraction

		Returns
		-------
			The supernova ejecta radius in cm
		'''
		return b * self.c * t

	def number_density(self, t, b = 1e-1, M = 1e1):
		'''
		Returns the supernova shell nucleon density.

		Parameters
		----------
		t : float
			The time passed from magnetar formation in s
		b : float, optional
			The relativistic velocity fraction
		M : float, optional
			The total ejecta mass in solar masses

		Returns
		-------
			The supernova shell nucleon density in 1 / cm**3
		'''
		r = self.ejecta_radius(t, b)
		return 3 * M * 1.9884e30/ (4 * np.pi * r**3 * 1.672621926e-27)

	def _hadron_spectrum(self, t, E, h, f = 1e-1, b = 1e-1, M = 1e1, D = False, N = 100):
		'''
		Returns the hadron production spectrum from injection of protons.

		Parameters
		----------
		t : float
			The time passed from magnetar formation in s
		E : float
			The projectile energy Eh as viewed from target rest coordinates in GeV
		h : {'pi', 'k', 'd0', 'd+', 'd+s', 'lam+c'}
			The type of hadronic particle
		f : float, optional
			The efficiency fraction of potential drop acceleration
		b : float, optional
			The relativistic velocity fraction
		M : float, optional
			The total ejecta mass in solar masses
		D : bool, optional
			The option to consider ejecta size for cooling, assumed to be infinite if `False`
		N : int, optional
			The number of steps for integration accuracy

		Returns
		-------
			The hadron production spectrum from injection of protons in 1 / (GeV s)
		'''
		Ep = self.E(t, f)
		x = E / Ep
		if x < 0 or x > 1:
			warn(f'{x} is out of bounds {0.0} to {1.0}')
			return 0.0
		n = self.number_density(t, b, M)
		d = self.ejecta_radius(t, b)
		if D is True:
			cf = hadron_proton_cooling_factor(E, n, h, d)
		else:
			cf = hadron_proton_cooling_factor(E, n, h)
		od = proton_proton_optical_depth(Ep, n, d)
		sig = self.proton_spectrum_prefactor(t)
		match h.lower():
			case 'pi':
				prod = meson_production(x, Ep, h)
			case 'k':
				prod = meson_production(x, Ep, h)
			case 'd0':
				prod = charmed_hadron_production(x, Ep, h, N)
			case 'd+':
				prod = charmed_hadron_production(x, Ep, h, N)
			case 'd+s':
				prod = charmed_hadron_production(x, Ep, h, N)
			case 'lam+c':
				prod = charmed_hadron_production(x, Ep, h, N)
			case _:
				raise ValueError(f'`{h.lower()}` is an invalid hadron identifyer, use `pi`, `k`, `d0`, `d+`, `d+s` or `lam+c` instead')
		return cf * od * sig * prod

	def hadron_spectrum(self, t, E, h, f = 1e-1, b = 1e-1, M = 1e1, D = False, N = 100):
		'''
		Returns the hadron spectrum from injection of protons.

		Parameters
		----------
		t : float
			The time passed from magnetar formation in s
		E : float
			The energy Eh as viewed from target rest coordinates in GeV
		h : {'pi', 'k', 'd0', 'd+', 'd+s', 'lam+c'}
			The type of hadronic particle
		f : float, optional
			The efficiency fraction of potential drop acceleration
		b : float, optional
			The relativistic velocity fraction
		M : float, optional
			The total ejecta mass in solar masses
		D : bool, optional
			The option to consider ejecta size for cooling, assumed to be infinite if `False`
		N : int, optional
			The number of steps for integration accuracy

		Returns
		-------
			The hadron production spectrum from injection of protons in 1 / (GeV s)
		'''
		return np.vectorize(self._hadron_spectrum)(t, E, h, f, b, M, D, N)

	def _neutrino_spectrum(self, t, E, h, f = 1e-1, b = 1e-1, M = 1e1, D = False, N = 100):
		'''
		Returns the neutrino spectrum from decay of hadrons.

		Parameters
		----------
		t : float
			The time passed from magnetar formation in s
		E : float
			The energy Enu as viewed from target rest coordinates in GeV
		h : {'pi', 'k', 'd0', 'd+', 'd+s', 'lam+c'}
			The type of hadronic particle
		f : float, optional
			The efficiency fraction of potential drop acceleration
		b : float, optional
			The relativistic velocity fraction
		M : float, optional
			The total ejecta mass in solar masses
		D : bool, optional
			The option to consider ejecta size for cooling, assumed to be infinite if `False`
		N : int, optional
			The number of steps for integration accuracy

		Returns
		-------
			The neutrino spectrum from decay of hadrons in 1 / (GeV s)
		'''
		match h.lower():
			case 'pi':
				l = 0.106**2 / 0.140**2
			case 'k':
				l = 0.106**2 / 0.494**2
			case 'd0':
				l = 0.67**2 / 1.86**2
			case 'd+':
				l = 0.63**2 / 1.87**2
			case 'd+s':
				l = 0.84**2 / 1.97**2
			case 'lam+c':
				l = 1.27**2 / 2.29**2
			case _:
				raise ValueError(f'`{h.lower()}` is an invalid hadron identifyer, use `pi`, `k`, `d0`, `d+`, `d+s` or `lam+c` instead')
		Ep = self.E(t, f)
		Eh = np.linspace(E / (1 - l), Ep, N)
		dEh = Eh[1] - Eh[0]
		dsig = (self.hadron_spectrum(t, Eh, h, f, b, M, D, N))[1:]
		match h.lower():
			case 'pi':
				ddec = (meson_decay_neutrinos(E, Eh, h))[1:]
			case 'k':
				ddec = (meson_decay_neutrinos(E, Eh, h))[1:]
			case 'd0':
				ddec = (charmed_hadron_decay_neutrinos(E, Eh, h))[1:]
			case 'd+':
				ddec = (charmed_hadron_decay_neutrinos(E, Eh, h))[1:]
			case 'd+s':
				ddec = (charmed_hadron_decay_neutrinos(E, Eh, h))[1:]
			case 'lam+c':
				ddec = (charmed_hadron_decay_neutrinos(E, Eh, h))[1:]
			case _:
				raise ValueError(f'`{h.lower()}` is an invalid hadron identifyer, use `pi`, `k`, `d0`, `d+`, `d+s` or `lam+c` instead')
		return np.sum(dEh * dsig * ddec)

	def neutrino_spectrum(self, t, E, h, f = 1e-1, b = 1e-1, M = 1e1, D = False, N = 100):
		'''
		Returns the neutrino spectrum from decay of hadrons.

		Parameters
		----------
		t : float
			The time passed from magnetar formation in s
		E : float
			The energy Enu as viewed from target rest coordinates in GeV
		h : {'pi', 'k', 'd0', 'd+', 'd+s', 'lam+c'}
			The type of hadronic particle
		f : float, optional
			The efficiency fraction of potential drop acceleration
		b : float, optional
			The relativistic velocity fraction
		M : float, optional
			The total ejecta mass in solar masses
		D : bool, optional
			The option to consider ejecta size for cooling, assumed to be infinite if `False`
		N : int, optional
			The number of steps for integration accuracy

		Returns
		-------
			The neutrino spectrum from decay of hadrons in 1 / (GeV s)
		'''
		return np.vectorize(self._neutrino_spectrum)(t, E, h, f, b, M, D, N)


mag = magnetar()
print(mag)

import matplotlib.pyplot as plt

t = np.logspace(1, 5, 10)

spec = mag._neutrino_spectrum(mag.tsd, 1e9, 'k')

print(spec)
#plt.plot(t, mag.neutrino_spectrum(t, 1e9, 'k'))
#plt.plot(t, mag.neutrino_spectrum(t, 1e9, 'd0', D = True))

#plt.xscale('log')
#plt.yscale('log')

#plt.show()






