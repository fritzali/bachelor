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
		Eh = np.linspace(E / (0.9999 * (1 - l)), 0.9999 * Ep, 100 * N)
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

	def _integrated_neutrino_spectrum(self, t1, t2, E, h, f = 1e-1, b = 1e-1, M = 1e1, D = False, N = 100):
		'''
		Returns the integrated neutrino spectrum from decay of hadrons.

		Parameters
		----------
		t1 : float
			The lower integration limit as time passed from magnetar formation in s
		t2 : float
			The upper integration limit as time passed from magnetar formation in s
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
			The integrated neutrino spectrum from decay of hadrons in 1 / GeV
		'''
		t = np.linspace(t1, t2, 20)
		dt = t[1] - t[0]
		dspec = (self.neutrino_spectrum(t, E, h, f, b, M, D, N))[1:]
		return np.sum(dt * dspec)

	def integrated_neutrino_spectrum(self, t1, t2, E, h, f = 1e-1, b = 1e-1, M = 1e1, D = False, N = 100):
		'''
		Returns the integrated neutrino spectrum from decay of hadrons.

		Parameters
		----------
		t1 : float
			The lower integration limit as time passed from magnetar formation in s
		t2 : float
			The upper integration limit as time passed from magnetar formation in s
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
			The integrated neutrino spectrum from decay of hadrons in 1 / GeV
		'''
		return np.vectorize(self._integrated_neutrino_spectrum)(t1, t2, E, h, f, b, M, D, N)


mag = magnetar()

def magnetar_hadron_spectrum(mag, Kt = 1000, KE = 100, f = 1e-1, b = 1e-1, M = 1e1, D = False, N = 100):
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
	with open('code/tabulate/hadrons/pi.txt', 'w') as f:
		spec = mag.hadron_spectrum(t[None, :], E[:, None], 'pi')
		f.write(f'# `pi` Hadron Spectrum / 1/(GeVs) - {datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")}\n')
		f.write(f'# Time / s (horizontal axis)\n')
		f.write(f'# Energy / GeV (vertical axis)\n')
		np.savetxt(f, spec)
	end = time.perf_counter()
	with open('code/tabulate/hadrons/axes.txt', 'w') as f:
		f.write(f'# Hadron Spectrum / 1/(GeVs) - {datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")}\n')
		f.write(f'{mag}')
		f.write(f'\n# Time / s (horizontal axis)\n')
		np.savetxt(f, t, newline=' ')
		f.write(f'\n# Energy / GeV (vertical axis)\n')
		np.savetxt(f, E, newline=' ')
		f.write(f'\n# Elapsed Time / s\n')
		f.write(f'# {end - start}')


def magnetar_neutrino_spectrum(mag, K = 1000):
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
	t = np.genfromtxt('code/tabulate/hadrons/axes.txt', skip_footer=1)
	x = np.genfromtxt('code/tabulate/hadrons/axes.txt', skip_header=15)
	d = np.diag(np.insert((x[1:] - x[:-1]), 0, 0.0))
	with open('code/tabulate/neutrinos/pi.txt', 'w') as f:
		had = np.genfromtxt('code/tabulate/hadrons/pi.txt')
		dec = meson_decay_neutrinos(x[None, :], E[:, None], 'pi')
		spec = dec @ d @ had
		f.write(f'# `pi` Neutrino Spectrum / 1/(GeVs) - {datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")}\n')
		f.write(f'# Time / s (horizontal axis)\n')
		f.write(f'# Energy / GeV (vertical axis)\n')
		np.savetxt(f, spec)
	end = time.perf_counter()
	with open('code/tabulate/neutrinos/axes.txt', 'w') as f:
		f.write(f'# Neutrino Spectrum / 1/(GeVs) - {datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")}\n')
		f.write(f'{mag}')
		f.write(f'\n# Time / s (horizontal axis)\n')
		np.savetxt(f, t, newline=' ')
		f.write(f'\n# Energy / GeV (vertical axis)\n')
		np.savetxt(f, E, newline=' ')
		f.write(f'\n# Elapsed Time / s\n')
		f.write(f'# {end - start}')
	

magnetar_hadron_spectrum(mag)
magnetar_neutrino_spectrum(mag)



# t = np.logspace(1, 7, 50)
# 
# pi = mag.hadron_spectrum(t, 1e8, 'pi')
# k = mag.hadron_spectrum(t, 1e8, 'k')
# c = (mag.hadron_spectrum(t, 1e8, 'd0')
# 	+ mag.hadron_spectrum(t, 1e8, 'd+')
# 	+ mag.hadron_spectrum(t, 1e8, 'd+s')
# 	+ mag.hadron_spectrum(t, 1e8, 'lam+c'))
# 
# import matplotlib.pyplot as plt
# 
# plt.plot(t, pi)
# plt.plot(t, k)
# plt.plot(t, c)
# 
# plt.xscale('log')
# plt.yscale('log')
# 
# plt.xlim(5e1, 1e7)
# plt.ylim(1e23, 5e27)
# 
# plt.show()

# t_E_1e6 = np.logspace(3, 7, 100)

# pi_had_E_1e6 = mag.hadron_spectrum(t_E_1e6, 1e6, 'pi')
# k_had_E_1e6 = mag.hadron_spectrum(t_E_1e6, 1e6, 'k')
# c_had_E_1e6 = (mag.hadron_spectrum(t_E_1e6, 1e6, 'd0')
# 			+ mag.hadron_spectrum(t_E_1e6, 1e6, 'd+')
# 			+ mag.hadron_spectrum(t_E_1e6, 1e6, 'd+s')
# 			+ mag.hadron_spectrum(t_E_1e6, 1e6, 'lam+c'))
# with open('data/had_E_1e6_no_cf.txt', 'w') as f:
# 	f.write('# Eh = 1e6 GeV (no cooling)\n')
# 	f.write('# t / s # pi / 1/(GeV s) # k / 1/(GeV s) # c / 1/(GeV s)\n')
# 	for row in zip(t_E_1e6, pi_had_E_1e6, k_had_E_1e6, c_had_E_1e6):
# 		f.write(r'{0} {1} {2} {3}'.format(*row))
# 		f.write('\n')

# pi_neu_E_1e6 = mag.neutrino_spectrum(t_E_1e6, 1e6, 'pi')
# k_neu_E_1e6 = mag. neutrino_spectrum(t_E_1e6, 1e6, 'k')
# d0_neu_E_1e6 = mag.neutrino_spectrum(t_E_1e6, 1e6, 'd0')
# with open('data/neu_E_1e6_no_cf.txt', 'w') as f:
# 	f.write('# Enu = 1e6 GeV (no cooling)\n')
# 	f.write('# t / s # pi / 1/(GeV s) # k / 1/(GeV s) # d0 / 1/(GeV s)\n')
# 	for row in zip(t_E_1e6, pi_neu_E_1e6, k_neu_E_1e6, d0_neu_E_1e6):
# 		f.write(r'{0} {1} {2} {3}'.format(*row))
# 		f.write('\n')


# t_E_1e9 = np.logspace(2, 5, 100)
 
# pi_had_E_1e9 = mag.hadron_spectrum(t_E_1e9, 1e9, 'pi')
# k_had_E_1e9 = mag.hadron_spectrum(t_E_1e9, 1e9,  'k')
# c_had_E_1e9 = (mag.hadron_spectrum(t_E_1e9, 1e9, 'd0')
# 			+ mag.hadron_spectrum(t_E_1e9, 1e9, 'd+')
# 			+ mag.hadron_spectrum(t_E_1e9, 1e9, 'd+s')
# 			+ mag.hadron_spectrum(t_E_1e9, 1e9, 'lam+c'))
# with open('data/had_E_1e9_no_cf.txt', 'w') as f:
# 	f.write('# Eh = 1e9 GeV (no cooling)\n')
# 	f.write('# t / s # pi / 1/(GeV s) # k / 1/(GeV s) # c / 1/(GeV s)\n')
# 	for row in zip(t_E_1e9, pi_had_E_1e9, k_had_E_1e9, c_had_E_1e9):
# 		f.write(r'{0} {1} {2} {3}'.format(*row))
# 		f.write('\n')

# pi_neu_E_1e9 = mag.neutrino_spectrum(t_E_1e9, 1e9, 'pi')
# k_neu_E_1e9 = mag.neutrino_spectrum(t_E_1e9, 1e9, 'k')
# d0_neu_E_1e9 = mag.neutrino_spectrum(t_E_1e9, 1e9, 'd0')
# with open('data/neu_E_1e9_with_cf.txt', 'w') as f:
# 	f.write('# Enu = 1e9 GeV (with cooling)\n')
# 	f.write('# t / s # pi / 1/(GeV s) # k / 1/(GeV s) # d0 / 1/(GeV s)\n')
# 	for row in zip(t_E_1e9, pi_neu_E_1e9, k_neu_E_1e9, d0_neu_E_1e9):
# 		f.write(r'{0} {1} {2} {3}'.format(*row))
# 		f.write('\n')


# E = np.logspace(4, 11, 20)

# pi_neu = mag.integrated_neutrino_spectrum(1e2, 1e7, E, 'pi')
# k_neu = mag.integrated_neutrino_spectrum(1e2, 1e7, E, 'k')
# d0_neu = mag.integrated_neutrino_spectrum(1e2, 1e7, E, 'd0')
# with open('data/neu_with_cf.txt', 'w') as f:
# 	f.write('# t1, t2 = 1e2, 1e7 (with cooling)\n')
# 	f.write('# E / GeV # pi / 1/GeV # k / 1/GeV # d0 / 1/GeV\n')
# 	for row in zip(E, pi_neu, k_neu, d0_neu):
# 		f.write(r'{0} {1} {2} {3}'.format(*row))
# 		f.write('\n')
