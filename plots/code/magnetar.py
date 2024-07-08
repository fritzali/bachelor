'''
Object oriented implementation of magnetar class as described in the thesis document.

	Classes
	-------
	magnetar
		Collects parameters and methods associated with the magnetar model

'''
import numpy as np

import integration as cal


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
		str3 = f'tsd = {self.tsd:.3} s\n    lum = {self.lum:.3} erg / s\n    E = {self.E(0):.3} GeV'
		return str1 + str2 + str3

	def L(self, t):
		'''
		Returns the spindown luminosity.

		Parameters
		----------
		t : float
			The time passed from magnetar formation

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
			The time passed from magnetar formation
		f : float, optional
			The efficiency fraction of potential drop acceleration

		Returns
		-------
			The monotonic proton energy in GeV
		'''
		return f * self.e * self.B * self.R**3 * self.o**2 / (2 * self.c**2 * (1 + t / self.tsd)) * 624.150907

	def proton_spectrum_coefficient(self, t, f = 1e-1):
		'''
		Returns the prefactor of a delta functional proton spectrum.

		Parameters
		----------
		t : float
			The time passed from magnetar formation
		f : float, optional
			The efficiency fraction of potential drop acceleration

		Returns
		-------
			The dimensionless prefactor of a delta functional proton spectrum
		'''
		return f * self.e * self.B * self.R**3 * self.o**2 / (2 * self.c**2 * (1 + t / self.tsd)) * 624.150907
