'''
Object oriented implementation of magnetar class as described in the thesis document.

	Classes
	-------
	magnetar
		Collects parameters and methods associated with the magnetar model

'''
import numpy as np

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






mag = magnetar()
print(mag)
