'''
Object oriented implementation of magnetar class as described in the thesis document.

	Classes
	-------
	magnetar
		Collects parameters and methods associated with the magnetar model

'''
import numpy as np
from warnings import warn

import datetime
import time

from functional import *


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
		str1 = f'# Magnetar:\n#     R = {self.R:.3} cm\n#     B = {self.B:.3} G\n#     o = {self.o:.3} rad / s\n#     '
		str2 = f'chi = {self.chi:.3} rad\n#     I = {self.I:.3} g * cm**2\n#     mu = {self.mu:.3} erg / G\n#     '
		str3 = f'tsd = {self.tsd:.3} s\n#     lum = {self.lum:.3} erg / s\n#     E = {self.E(0):.3} GeV\n#     '
		str4 = f'spec = {self.proton_spectrum_prefactor(0):.3}\n#     n = {self.number_density(self.tsd):.3} 1 / cm**3'
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

	def _cooling_factor(self, t, E, h, b = 1e-1, M = 1e1, D = False):
		'''
		Returns the ejecta material cooling factor.

		Parameters
		----------
		t : float
			The time passed from magnetar formation in s
		E : float
			The projectile energy Eh as viewed from target rest coordinates in GeV
		h : {'pi', 'k', 'd0', 'd+', 'd+s', 'lam+c'}
			The type of hadronic particle
		b : float, optional
			The relativistic velocity fraction
		M : float, optional
			The total ejecta mass in solar masses
		D : bool, optional
			The option to consider ejecta size for cooling, assumed to be infinite if `False`

		Returns
		-------
			The unitless ejecta material cooling factor
		'''
		n = self.number_density(t, b, M)
		d = self.ejecta_radius(t, b)
		if D is True:
			return hadron_proton_cooling_factor(E, n, h, d)
		else:
			return hadron_proton_cooling_factor(E, n, h)

	def cooling_factor(self, t, E, h, b = 1e-1, M = 1e1, D = False):
		'''
		Returns the ejecta material cooling factor.

		Parameters
		----------
		t : float
			The time passed from magnetar formation in s
		E : float
			The projectile energy Eh as viewed from target rest coordinates in GeV
		h : {'pi', 'k', 'd0', 'd+', 'd+s', 'lam+c'}
			The type of hadronic particle
		b : float, optional
			The relativistic velocity fraction
		M : float, optional
			The total ejecta mass in solar masses
		D : bool, optional
			The option to consider ejecta size for cooling, assumed to be infinite if `False`

		Returns
		-------
			The unitless ejecta material cooling factor
		'''
		return np.vectorize(self._cooling_factor)(t, E, h, b, M, D)

	def _optical_depth(self, t, f = 1e-1, b = 1e-1, M = 1e1):
		'''
		Returns the ejecta material proton effective optical depth.

		Parameters
		----------
		t : float
			The time passed from magnetar formation in s
		f : float, optional
			The efficiency fraction of potential drop acceleration
		b : float, optional
			The relativistic velocity fraction
		M : float, optional
			The total ejecta mass in solar masses

		Returns
		-------
			The unitless ejecta material proton effective optical depth
		'''
		Ep = self.E(t, f)
		n = self.number_density(t, b, M)
		d = self.ejecta_radius(t, b)
		return proton_proton_optical_depth(Ep, n, d)

	def optical_depth(self, t, f = 1e-1, b = 1e-1, M = 1e1):
		'''
		Returns the ejecta material proton effective optical depth.

		Parameters
		----------
		t : float
			The time passed from magnetar formation in s
		f : float, optional
			The efficiency fraction of potential drop acceleration
		b : float, optional
			The relativistic velocity fraction
		M : float, optional
			The total ejecta mass in solar masses

		Returns
		-------
			The unitless ejecta material proton effective optical depth
		'''
		return np.vectorize(self._optical_depth)(t, f, b, M)

	def collision_factor(self, t, E, h, f = 1e-1, b = 1e-1, M = 1e1, D = False):
		'''
		Returns the ejecta material combined attenuation factor.

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

		Returns
		-------
			The unitless ejecta material combined attenuation factor
		'''
		cf = self.cooling_factor(t, E, h, b, M, D)
		od = self.optical_depth(t, f, b, M)
		return cf * od

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
		sig = self.proton_spectrum_prefactor(t)
		f = self.collision_factor(t, E, h, f, b, M, D)
		return prod * sig * f

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


def magnetar_hadron_spectrum(mag, Kt = 500, KE = 100, f = 1e-1, b = 1e-1, M = 1e1, D = False, N = 100):
	'''
	Prints calculated hadron spectra for all types to tabulated text files

	Parameters
	----------
	mag : magnetar
		The magnetar object of which respective methods are used
	Kt : int, optional
		The number of points in time
	KE : int, optional
		The number of energy values
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
		None
	'''
	start = time.perf_counter()
	t = np.logspace(1, 8, Kt)
	E = np.logspace(5, 12, KE)
	with open('code/tabulate/magnetar/hadrons/pi.txt', 'w') as f:
		spec = mag.hadron_spectrum(t[None, :], E[:, None], 'pi')
		f.write(f'# `pi` Hadron Spectrum / 1/(GeVs) - {datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")}\n')
		f.write(f'# Time / s (horizontal axis)\n')
		f.write(f'# Energy / GeV (vertical axis)\n')
		np.savetxt(f, spec)
	with open('code/tabulate/magnetar/hadrons/K.txt', 'w') as f:
		spec = mag.hadron_spectrum(t[None, :], E[:, None], 'k')
		f.write(f'# `K` Hadron Spectrum / 1/(GeVs) - {datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")}\n')
		f.write(f'# Time / s (horizontal axis)\n')
		f.write(f'# Energy / GeV (vertical axis)\n')
		np.savetxt(f, spec)
	with open('code/tabulate/magnetar/hadrons/D0.txt', 'w') as f:
		spec = mag.hadron_spectrum(t[None, :], E[:, None], 'd0')
		f.write(f'# `D0` Hadron Spectrum / 1/(GeVs) - {datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")}\n')
		f.write(f'# Time / s (horizontal axis)\n')
		f.write(f'# Energy / GeV (vertical axis)\n')
		np.savetxt(f, spec)
	with open('code/tabulate/magnetar/hadrons/Dplus.txt', 'w') as f:
		spec = mag.hadron_spectrum(t[None, :], E[:, None], 'd+')
		f.write(f'# `D+` Hadron Spectrum / 1/(GeVs) - {datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")}\n')
		f.write(f'# Time / s (horizontal axis)\n')
		f.write(f'# Energy / GeV (vertical axis)\n')
		np.savetxt(f, spec)
	with open('code/tabulate/magnetar/hadrons/DplusS.txt', 'w') as f:
		spec = mag.hadron_spectrum(t[None, :], E[:, None], 'd+s')
		f.write(f'# `D+s` Hadron Spectrum / 1/(GeVs) - {datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")}\n')
		f.write(f'# Time / s (horizontal axis)\n')
		f.write(f'# Energy / GeV (vertical axis)\n')
		np.savetxt(f, spec)
	with open('code/tabulate/magnetar/hadrons/LAMplusC.txt', 'w') as f:
		spec = mag.hadron_spectrum(t[None, :], E[:, None], 'lam+c')
		f.write(f'# `Lam+c` Hadron Spectrum / 1/(GeVs) - {datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")}\n')
		f.write(f'# Time / s (horizontal axis)\n')
		f.write(f'# Energy / GeV (vertical axis)\n')
		np.savetxt(f, spec)
	end = time.perf_counter()
	with open('code/tabulate/magnetar/hadrons/axes.txt', 'w') as f:
		f.write(f'# Hadron Spectrum / 1/(GeVs) - {datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")}\n')
		f.write(f'{mag}')
		f.write(f'\n# Time / s (horizontal axis)\n')
		np.savetxt(f, t, newline=' ')
		f.write(f'\n# Energy / GeV (vertical axis)\n')
		np.savetxt(f, E, newline=' ')
		f.write(f'\n# Elapsed Time / s\n')
		f.write(f'# {end - start}')


def magnetar_neutrino_spectrum(mag, K = 100):
	'''
	Prints calculated neutrino spectra for all types to tabulated text files

	Parameters
	----------
	mag : magnetar
		The magnetar object of which respective methods are used
	K : int, optional
		The number of energy values

	Returns
	-------
		None
	'''
	start = time.perf_counter()
	E = np.logspace(5, 12, K)
	t = np.genfromtxt('code/tabulate/magnetar/hadrons/axes.txt', skip_footer=1)
	x = np.genfromtxt('code/tabulate/magnetar/hadrons/axes.txt', skip_header=15)
	d = np.diag(np.insert((x[1:] - x[:-1]), 0, 0.0))
	with open('code/tabulate/magnetar/neutrinos/pi.txt', 'w') as f:
		f.write(f'# `pi` Neutrino Spectrum / 1/(GeVs) - {datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")}\n')
		f.write(f'# Time / s (horizontal axis)\n')
		f.write(f'# Energy / GeV (vertical axis)\n')
		had = np.genfromtxt('code/tabulate/magnetar/hadrons/pi.txt')
		dec = meson_decay_neutrinos(E[:, None], x[None, :], 'pi')
		spec = dec @ d @ had
		np.savetxt(f, spec)
	with open('code/tabulate/magnetar/neutrinos/K.txt', 'w') as f:
		f.write(f'# `K` Neutrino Spectrum / 1/(GeVs) - {datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")}\n')
		f.write(f'# Time / s (horizontal axis)\n')
		f.write(f'# Energy / GeV (vertical axis)\n')
		had = np.genfromtxt('code/tabulate/magnetar/hadrons/K.txt')
		dec = meson_decay_neutrinos(E[:, None], x[None, :], 'k')
		spec = dec @ d @ had
		np.savetxt(f, spec)
	with open('code/tabulate/magnetar/neutrinos/D0.txt', 'w') as f:
		f.write(f'# `D0` Neutrino Spectrum / 1/(GeVs) - {datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")}\n')
		f.write(f'# Time / s (horizontal axis)\n')
		f.write(f'# Energy / GeV (vertical axis)\n')
		had = np.genfromtxt('code/tabulate/magnetar/hadrons/D0.txt')
		dec = charmed_hadron_decay_neutrinos(E[:, None], x[None, :], 'd0')
		spec = dec @ d @ had
		np.savetxt(f, spec)
	with open('code/tabulate/magnetar/neutrinos/Dplus.txt', 'w') as f:
		f.write(f'# `D+` Neutrino Spectrum / 1/(GeVs) - {datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")}\n')
		f.write(f'# Time / s (horizontal axis)\n')
		f.write(f'# Energy / GeV (vertical axis)\n')
		had = np.genfromtxt('code/tabulate/magnetar/hadrons/Dplus.txt')
		dec = charmed_hadron_decay_neutrinos(E[:, None], x[None, :], 'd+')
		spec = dec @ d @ had
		np.savetxt(f, spec)
	with open('code/tabulate/magnetar/neutrinos/DplusS.txt', 'w') as f:
		f.write(f'# `D+s` Neutrino Spectrum / 1/(GeVs) - {datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")}\n')
		f.write(f'# Time / s (horizontal axis)\n')
		f.write(f'# Energy / GeV (vertical axis)\n')
		had = np.genfromtxt('code/tabulate/magnetar/hadrons/DplusS.txt')
		dec = charmed_hadron_decay_neutrinos(E[:, None], x[None, :], 'd+s')
		spec = dec @ d @ had
		np.savetxt(f, spec)
	with open('code/tabulate/magnetar/neutrinos/LAMplusC.txt', 'w') as f:
		f.write(f'# `Lam+c` Neutrino Spectrum / 1/(GeVs) - {datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")}\n')
		f.write(f'# Time / s (horizontal axis)\n')
		f.write(f'# Energy / GeV (vertical axis)\n')
		had = np.genfromtxt('code/tabulate/magnetar/hadrons/LAMplusC.txt')
		dec = charmed_hadron_decay_neutrinos(E[:, None], x[None, :], 'lam+c')
		spec = dec @ d @ had
		np.savetxt(f, spec)
	end = time.perf_counter()
	with open('code/tabulate/magnetar/neutrinos/axes.txt', 'w') as f:
		f.write(f'# Neutrino Spectrum / 1/(GeVs) - {datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")}\n')
		f.write(f'{mag}')
		f.write(f'\n# Time / s (horizontal axis)\n')
		np.savetxt(f, t, newline=' ')
		f.write(f'\n# Energy / GeV (vertical axis)\n')
		np.savetxt(f, E, newline=' ')
		f.write(f'\n# Elapsed Time / s\n')
		f.write(f'# {end - start}')


def magnetar_integrated_neutrino_spectrum(mag):
	'''
	Prints integrated neutrino spectra for all types to tabulated text files

	Parameters
	----------
	mag : magnetar
		The magnetar object of which respective methods are used
	ta : float
		The lower temproal bound of integration
	tb : float
		The upper temproal bound of integration

	Returns
	-------
		None
	'''
	t = np.genfromtxt('code/tabulate/magnetar/neutrinos/axes.txt', skip_footer=1)
	str1 = '(1e2 - 1e4)'
	str2 = '(1e4 - 1e5)'
	str3 = '(1e2 - 1e7)'
	con1 = (t > 1e2) & (t < 1e4)
	con2 = (t > 1e4) & (t < 1e5)
	con3 = (t > 1e2) & (t < 1e7)
	t1 = t[con1]
	t2 = t[con2]
	t3 = t[con3]
	E = np.genfromtxt('code/tabulate/magnetar/neutrinos/axes.txt', skip_header=15)
	with open('code/tabulate/magnetar/integrate/pi.txt', 'w') as f:
		f.write(f'# `pi` Integrated Spectrum / 1/GeV - {datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")}\n')
		f.write(f'{mag}\n')
		f.write(f'# Energy / GeV ')
		f.write(f'# Spectrum {str1} / 1/GeV ')
		f.write(f'# Spectrum {str2} / 1/GeV ')
		f.write(f'# Spectrum {str3} / 1/GeV\n')
		neu = np.genfromtxt('code/tabulate/magnetar/neutrinos/pi.txt')
		y1 = neu[:, con1]
		y2 = neu[:, con2]
		y3 = neu[:, con3]
		spec1 = np.trapezoid(y1, t1, axis=1)
		spec2 = np.trapezoid(y2, t2, axis=1)
		spec3 = np.trapezoid(y3, t3, axis=1)
		for row in zip(E, spec1, spec2, spec3):
			f.write(r'{0}   {1}   {2}   {3}'.format(*row))
			f.write(f'\n')
	with open('code/tabulate/magnetar/integrate/K.txt', 'w') as f:
		f.write(f'# `K` Integrated Spectrum / 1/GeV - {datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")}\n')
		f.write(f'{mag}\n')
		f.write(f'# Energy / GeV ')
		f.write(f'# Spectrum {str1} / 1/GeV ')
		f.write(f'# Spectrum {str2} / 1/GeV ')
		f.write(f'# Spectrum {str3} / 1/GeV\n')
		neu = np.genfromtxt('code/tabulate/magnetar/neutrinos/K.txt')
		y1 = neu[:, con1]
		y2 = neu[:, con2]
		y3 = neu[:, con3]
		spec1 = np.trapezoid(y1, t1, axis=1)
		spec2 = np.trapezoid(y2, t2, axis=1)
		spec3 = np.trapezoid(y3, t3, axis=1)
		for row in zip(E, spec1, spec2, spec3):
			f.write(r'{0}   {1}   {2}   {3}'.format(*row))
			f.write(f'\n')
	with open('code/tabulate/magnetar/integrate/D0.txt', 'w') as f:
		f.write(f'# `D0` Integrated Spectrum / 1/GeV - {datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")}\n')
		f.write(f'{mag}\n')
		f.write(f'# Energy / GeV ')
		f.write(f'# Spectrum {str1} / 1/GeV ')
		f.write(f'# Spectrum {str2} / 1/GeV ')
		f.write(f'# Spectrum {str3} / 1/GeV\n')
		neu = np.genfromtxt('code/tabulate/magnetar/neutrinos/D0.txt')
		y1 = neu[:, con1]
		y2 = neu[:, con2]
		y3 = neu[:, con3]
		spec1 = np.trapezoid(y1, t1, axis=1)
		spec2 = np.trapezoid(y2, t2, axis=1)
		spec3 = np.trapezoid(y3, t3, axis=1)
		for row in zip(E, spec1, spec2, spec3):
			f.write(r'{0}   {1}   {2}   {3}'.format(*row))
			f.write(f'\n')
	with open('code/tabulate/magnetar/integrate/Dplus.txt', 'w') as f:
		f.write(f'# `D+` Integrated Spectrum / 1/GeV - {datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")}\n')
		f.write(f'{mag}\n')
		f.write(f'# Energy / GeV ')
		f.write(f'# Spectrum {str1} / 1/GeV ')
		f.write(f'# Spectrum {str2} / 1/GeV ')
		f.write(f'# Spectrum {str3} / 1/GeV\n')
		neu = np.genfromtxt('code/tabulate/magnetar/neutrinos/Dplus.txt')
		y1 = neu[:, con1]
		y2 = neu[:, con2]
		y3 = neu[:, con3]
		spec1 = np.trapezoid(y1, t1, axis=1)
		spec2 = np.trapezoid(y2, t2, axis=1)
		spec3 = np.trapezoid(y3, t3, axis=1)
		for row in zip(E, spec1, spec2, spec3):
			f.write(r'{0}   {1}   {2}   {3}'.format(*row))
			f.write(f'\n')
	with open('code/tabulate/magnetar/integrate/DplusS.txt', 'w') as f:
		f.write(f'# `D+s` Integrated Spectrum / 1/GeV - {datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")}\n')
		f.write(f'{mag}\n')
		f.write(f'# Energy / GeV ')
		f.write(f'# Spectrum {str1} / 1/GeV ')
		f.write(f'# Spectrum {str2} / 1/GeV ')
		f.write(f'# Spectrum {str3} / 1/GeV\n')
		neu = np.genfromtxt('code/tabulate/magnetar/neutrinos/DplusS.txt')
		y1 = neu[:, con1]
		y2 = neu[:, con2]
		y3 = neu[:, con3]
		spec1 = np.trapezoid(y1, t1, axis=1)
		spec2 = np.trapezoid(y2, t2, axis=1)
		spec3 = np.trapezoid(y3, t3, axis=1)
		for row in zip(E, spec1, spec2, spec3):
			f.write(r'{0}   {1}   {2}   {3}'.format(*row))
			f.write(f'\n')
	with open('code/tabulate/magnetar/integrate/LAMplusC.txt', 'w') as f:
		f.write(f'# `Lam+c` Integrated Spectrum / 1/GeV - {datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")}\n')
		f.write(f'{mag}\n')
		f.write(f'# Energy / GeV ')
		f.write(f'# Spectrum {str1} / 1/GeV ')
		f.write(f'# Spectrum {str2} / 1/GeV ')
		f.write(f'# Spectrum {str3} / 1/GeV\n')
		neu = np.genfromtxt('code/tabulate/magnetar/neutrinos/LAMplusC.txt')
		y1 = neu[:, con1]
		y2 = neu[:, con2]
		y3 = neu[:, con3]
		spec1 = np.trapezoid(y1, t1, axis=1)
		spec2 = np.trapezoid(y2, t2, axis=1)
		spec3 = np.trapezoid(y3, t3, axis=1)
		for row in zip(E, spec1, spec2, spec3):
			f.write(r'{0}   {1}   {2}   {3}'.format(*row))
			f.write(f'\n')


mag = magnetar(B = 10**14.5)


#magnetar_hadron_spectrum(mag)
#magnetar_neutrino_spectrum(mag)
magnetar_integrated_neutrino_spectrum(mag)
