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
		cf = hadron_proton_cooling_factor(Eh, n, 'pi')
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
		cf = hadron_proton_cooling_factor(Eh, n, 'k')
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
		cf = hadron_proton_cooling_factor(Eh, n, 'd0')
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
		cf = hadron_proton_cooling_factor(Eh, n, 'd+')
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
		cf = hadron_proton_cooling_factor(Eh, n, 'd+s')
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
		cf = hadron_proton_cooling_factor(Eh, n, 'lam+c')
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







# Fph_pi = meson_production(x, Ep, 'pi')
# Fph_K = meson_production(x, Ep, 'k')
# Fph_D0 = charmed_hadron_production(x, Ep, 'd0')
# Fph_Dp = charmed_hadron_production(x, Ep, 'd+')
# Fph_DpS = charmed_hadron_production(x, Ep, 'd+s')
# Fph_LpC = charmed_hadron_production(x, Ep, 'lam+c')
# 
# dp = np.diag(np.insert((Ep[1:] - Ep[:-1]), 0, 0.0))
# 
# Sh_pi = Fph_pi @ dp @ Sp
# Sh_K = Fph_K @ dp @ Sp
# Sh_D0 = Fph_D0 @ dp @ Sp
# Sh_Dp = Fph_Dp @ dp @ Sp
# Sh_DpS = Fph_DpS @ dp @ Sp
# Sh_LpC = Fph_LpC @ dp @ Sp
# 
# cf_pi = hadron_proton_cooling_factor(Eh, n, 'pi')
# cf_K = hadron_proton_cooling_factor(Eh, n, 'k')
# cf_D0 = hadron_proton_cooling_factor(Eh, n, 'd0')
# cf_Dp = hadron_proton_cooling_factor(Eh, n, 'd+')
# cf_DpS = hadron_proton_cooling_factor(Eh, n, 'd+s')
# cf_LpC = hadron_proton_cooling_factor(Eh, n, 'lam+c')
# 
# Sh_pi = np.nan_to_num(cf_pi * Sh_pi)
# Sh_K = np.nan_to_num(cf_K * Sh_K)
# Sh_D0 = np.nan_to_num(cf_D0 * Sh_D0)
# Sh_Dp = np.nan_to_num(cf_Dp * Sh_Dp)
# Sh_DpS = np.nan_to_num(cf_DpS * Sh_DpS)
# Sh_LpC = np.nan_to_num(cf_LpC * Sh_LpC)
# 
# Sh_c = Sh_D0 + Sh_Dp + Sh_DpS + Sh_LpC
# 
# Enu = np.logspace(5, 12, 1000)
# 
# Fhnu_pi = meson_decay_neutrinos(Enu[:, None], Eh[None, :], 'pi')
# Fhnu_K = meson_decay_neutrinos(Enu[:, None], Eh[None, :], 'K')
# Fhnu_D0 = charmed_hadron_decay_neutrinos(Enu[:, None], Eh[None, :], 'd0')
# Fhnu_Dp = charmed_hadron_decay_neutrinos(Enu[:, None], Eh[None, :], 'd+')
# Fhnu_DpS = charmed_hadron_decay_neutrinos(Enu[:, None], Eh[None, :], 'd+s')
# Fhnu_LpC = charmed_hadron_decay_neutrinos(Enu[:, None], Eh[None, :], 'lam+c')
# 
# dh = np.diag(np.insert((Eh[1:] - Eh[:-1]), 0, 0.0))
# 
# Snu_pi = Fhnu_pi @ dh @ Sh_pi
# Snu_K = Fhnu_K @ dh @ Sh_K
# Snu_D0 = Fhnu_D0 @ dh @ Sh_D0
# Snu_Dp = Fhnu_Dp @ dh @ Sh_Dp
# Snu_DpS = Fhnu_DpS @ dh @ Sh_DpS
# Snu_LpC = Fhnu_LpC @ dh @ Sh_LpC
# 
# Snu_c = Snu_D0 + Snu_Dp + Snu_DpS + Snu_LpC
# 
# # plt.plot(Ep, Ep**2 * Sp, '-', label='p')
# 
# # plt.plot(Eh, Eh**2 * Sh_pi, '--', label='pi')
# # plt.plot(Eh, Eh**2 * Sh_K, '-.', label='K')
# # plt.plot(Eh, Eh**2 * Sh_c, ':', label='c')
# 
# plt.plot(Enu, Enu**2 * Snu_pi, label='pi')
# plt.plot(Enu, Enu**2 * Snu_K, label='K')
# plt.plot(Enu, Enu**2 * Snu_c, label='c')
# plt.plot(Enu, Enu**2 * Snu_D0, ':', label='D0')
# plt.plot(Enu, Enu**2 * Snu_Dp, ':', label='Dp')
# plt.plot(Enu, Enu**2 * Snu_DpS, ':', label='DpS')
# plt.plot(Enu, Enu**2 * Snu_LpC, ':', label='LpC')
# 

E, pi = np.genfromtxt('code/tabulate/nucleus/neutrinos/pi.txt', unpack=True)
E, K = np.genfromtxt('code/tabulate/nucleus/neutrinos/K.txt', unpack=True)
E, D0 = np.genfromtxt('code/tabulate/nucleus/neutrinos/D0.txt', unpack=True)
E, Dplus = np.genfromtxt('code/tabulate/nucleus/neutrinos/Dplus.txt', unpack=True)
E, DplusS = np.genfromtxt('code/tabulate/nucleus/neutrinos/DplusS.txt', unpack=True)
E, LAMplusC = np.genfromtxt('code/tabulate/nucleus/neutrinos/LAMplusC.txt', unpack=True)

c = D0 + Dplus + DplusS + LAMplusC

plt.fill_between(E, 3 * E**2 * c, E**2 * c / 3, color='none', facecolor='b', alpha=0.25)
plt.plot(E, E**2 * c, 'b')
plt.plot(E, E**2 * pi, 'r')
plt.plot(E, E**2 * K, 'k')

plt.xscale('log')
plt.yscale('log')

plt.xlim(1e6, 1e11)
plt.ylim(1e27, 1e34)

plt.legend()

plt.show()
