'''
Functional implementation of nucleus module as described in the thesis document.

	Functions
	---------
		nucleus_hadron_spectrum
			Prints calculated hadron spectra for all types to tabulated text files

'''
import matplotlib.pyplot as plt
import numpy as np

import datetime

from code.functional import *

def nucleus_neutrino_spectrum(reg, n = 1e14, d =1e15):
	'''
	Prints calculated hadron spectra for all types to tabulated text files.

	Parameters
	----------
	reg : string
		The directory string to which files are saved
	n : float
		The ionized hydrogen number density
	d : float
		The distance or accretion disk height

	Returns
	-------
		None
	'''
	Ep = np.logspace(5, 12, 200)
	dp = np.diag(np.insert((Ep[1:] - Ep[:-1]), 0, 0.0))
	Sp = 1e30 / Ep**2
	od = proton_proton_optical_depth(Ep, n, d)
	Sp = od * Sp
	Eh = np.logspace(5, 12, 1000)
	dh = np.diag(np.insert((Eh[1:] - Eh[:-1]), 0, 0.0))
	x = Eh[:, None] / Ep[None, :]
	Enu = np.logspace(5, 12, 1000)
	with open(f'{reg}/pi.txt', 'w') as file:
		file.write(f'# `pi` Integrated Spectrum / 1/GeV - {datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")}\n')
		file.write(f'# Nucleus:\n')
		file.write(f'#     n = {n:.3} 1/cm**3\n')
		file.write(f'#     d = {d:.3} cm\n')
		file.write(f'# Energy / GeV ')
		file.write(f'# Spectrum / 1/GeV\n')
		Fph = meson_production(x, Ep, 'pi')
		Sh = Fph @ dp @ Sp
		cf = hadron_proton_cooling_factor(Eh, n, 'pi', d)
		Sh = np.nan_to_num(cf * Sh)
		Fhnu = meson_decay_neutrinos(Enu[:, None], Eh[None, :], 'pi')
		Snu = Fhnu @ dh @ Sh
		for row in zip(Enu, Snu):
			file.write(r'{0}   {1}'.format(*row))
			file.write(f'\n')
	with open(f'{reg}/K.txt', 'w') as file:
		file.write(f'# `K` Integrated Spectrum / 1/GeV - {datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")}\n')
		file.write(f'# Nucleus:\n')
		file.write(f'#     n = {n:.3} 1/cm**3\n')
		file.write(f'#     d = {d:.3} cm\n')
		file.write(f'# Energy / GeV ')
		file.write(f'# Spectrum / 1/GeV\n')
		Fph = meson_production(x, Ep, 'k')
		Sh = Fph @ dp @ Sp
		cf = hadron_proton_cooling_factor(Eh, n, 'k', d)
		Sh = np.nan_to_num(cf * Sh)
		Fhnu = meson_decay_neutrinos(Enu[:, None], Eh[None, :], 'k')
		Snu = Fhnu @ dh @ Sh
		for row in zip(Enu, Snu):
			file.write(r'{0}   {1}'.format(*row))
			file.write(f'\n')
	with open(f'{reg}/D0.txt', 'w') as file:
		file.write(f'# `D0` Integrated Spectrum / 1/GeV - {datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")}\n')
		file.write(f'# Nucleus:\n')
		file.write(f'#     n = {n:.3} 1/cm**3\n')
		file.write(f'#     d = {d:.3} cm\n')
		file.write(f'# Energy / GeV ')
		file.write(f'# Spectrum / 1/GeV\n')
		Fph = charmed_hadron_production(x, Ep, 'd0')
		Sh = Fph @ dp @ Sp
		cf = hadron_proton_cooling_factor(Eh, n, 'd0', d)
		Sh = np.nan_to_num(cf * Sh)
		Fhnu = charmed_hadron_decay_neutrinos(Enu[:, None], Eh[None, :], 'd0')
		Snu = Fhnu @ dh @ Sh
		for row in zip(Enu, Snu):
			file.write(r'{0}   {1}'.format(*row))
			file.write(f'\n')
	with open(f'{reg}/Dplus.txt', 'w') as file:
		file.write(f'# `D+` Integrated Spectrum / 1/GeV - {datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")}\n')
		file.write(f'# Nucleus:\n')
		file.write(f'#     n = {n:.3} 1/cm**3\n')
		file.write(f'#     d = {d:.3} cm\n')
		file.write(f'# Energy / GeV ')
		file.write(f'# Spectrum / 1/GeV\n')
		Fph = charmed_hadron_production(x, Ep, 'd+')
		Sh = Fph @ dp @ Sp
		cf = hadron_proton_cooling_factor(Eh, n, 'd+', d)
		Sh = np.nan_to_num(cf * Sh)
		Fhnu = charmed_hadron_decay_neutrinos(Enu[:, None], Eh[None, :], 'd+')
		Snu = Fhnu @ dh @ Sh
		for row in zip(Enu, Snu):
			file.write(r'{0}   {1}'.format(*row))
			file.write(f'\n')
	with open(f'{reg}/DplusS.txt', 'w') as file:
		file.write(f'# `D+s` Integrated Spectrum / 1/GeV - {datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")}\n')
		file.write(f'# Nucleus:\n')
		file.write(f'#     n = {n:.3} 1/cm**3\n')
		file.write(f'#     d = {d:.3} cm\n')
		file.write(f'# Energy / GeV ')
		file.write(f'# Spectrum / 1/GeV\n')
		Fph = charmed_hadron_production(x, Ep, 'd+s')
		Sh = Fph @ dp @ Sp
		cf = hadron_proton_cooling_factor(Eh, n, 'd+s', d)
		Sh = np.nan_to_num(cf * Sh)
		Fhnu = charmed_hadron_decay_neutrinos(Enu[:, None], Eh[None, :], 'd+s')
		Snu = Fhnu @ dh @ Sh
		for row in zip(Enu, Snu):
			file.write(r'{0}   {1}'.format(*row))
			file.write(f'\n')
	with open(f'{reg}/LAMplusC.txt', 'w') as file:
		file.write(f'# `Lam+c` Integrated Spectrum / 1/GeV - {datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")}\n')
		file.write(f'# Nucleus:\n')
		file.write(f'#     n = {n:.3} 1/cm**3\n')
		file.write(f'#     d = {d:.3} cm\n')
		file.write(f'# Energy / GeV ')
		file.write(f'# Spectrum / 1/GeV\n')
		Fph = charmed_hadron_production(x, Ep, 'lam+c')
		Sh = Fph @ dp @ Sp
		cf = hadron_proton_cooling_factor(Eh, n, 'lam+c', d)
		Sh = np.nan_to_num(cf * Sh)
		Fhnu = charmed_hadron_decay_neutrinos(Enu[:, None], Eh[None, :], 'lam+c')
		Snu = Fhnu @ dh @ Sh
		for row in zip(Enu, Snu):
			file.write(r'{0}   {1}'.format(*row))
			file.write(f'\n')

	print(f'\n# Default')
	print(f'# Nucleus:')
	print(f'#     n = {n:.3} 1/cm**3')
	print(f'#     d = {d:.3} cm\n')


nucleus_neutrino_spectrum('code/tabulate/nucleus/neutrinos')
