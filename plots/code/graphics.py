import numpy as np
import matplotlib.pyplot as plt

from code.functional import *



# fig, ax = plt.subplots(nrows=1, ncols=2)
# 
# fig.set_size_inches(7.0, 3.0)
# 
# s, p, pi, K = np.genfromtxt('code/tabulate/other/inelastic_scattering.txt', unpack=True)
# 
# ax[0].plot(s, p, 'b-', label=r'$p \kern-0.1pt p$')
# ax[0].plot(s, pi, 'r--', label=r'$\pi \kern-0.1pt p$')
# ax[0].plot(s, K, 'k:', label=r'$K \kern-0.2pt p$')
# 
# ax[0].set_xlabel(r'$s$ $\mathrel{/} \kern-0.1pt$ GeV$^2$')
# ax[0].set_ylabel(r'$\sigma_{h \kern-0.1pt p}$ $\mathrel{/}$ mb')
# 
# ax[0].set_xscale('log')
# 
# ax[0].set_xlim(1e2, 1e9)
# 
# ax[0].legend(handlelength=1.8)
# 
# x, y1, y2, y3 = np.genfromtxt('code/tabulate/other/sample_charm_hadron.txt', unpack=True)
# 
# ax[1].plot(x, x * y1, 'b-', label=r'$E_p = 10^{12}$ GeV')
# ax[1].plot(x, x * y2, 'r--', label=r'$E_p = 10^{10}$ GeV')
# ax[1].plot(x, x * y3, 'k:', label=r'$E_p = 10^{8}$ GeV')
# 
# ax[1].set_xlabel(r'$x_h \kern-0.8pt = \kern-0.5pt E_h \kern+0.5pt / E_p$')
# ax[1].set_ylabel(r'$x_h \kern+0.5pt d\sigma \kern-0.3pt / \kern-0.8pt dx_h \kern+0.2pt$ $\mathrel{/}$ mb')
# 
# ax[1].set_xscale('log')
# ax[1].set_yscale('log')
# 
# ax[1].set_xlim(1e-6, 1e0)
# ax[1].set_ylim(1e-4, 1e1)
# 
# ax[1].legend(handlelength=1.8)
# 
# plt.savefig('build/cross_sections.pdf')
# plt.savefig('build/cross_sections.png')
# plt.close()



t = np.genfromtxt('code/tabulate/magnetar/without/neutrinos/axes.txt', skip_footer=1)
E = np.genfromtxt('code/tabulate/magnetar/without/neutrinos/axes.txt', skip_header=15)

en = 1e9
i = (np.abs(E - en)).argmin()
en = E[i]

pi = np.genfromtxt('code/tabulate/magnetar/without/neutrinos/pi.txt')
K = np.genfromtxt('code/tabulate/magnetar/without/neutrinos/K.txt')
D0 = np.genfromtxt('code/tabulate/magnetar/without/neutrinos/D0.txt')
Dplus = np.genfromtxt('code/tabulate/magnetar/without/neutrinos/Dplus.txt')
DplusS = np.genfromtxt('code/tabulate/magnetar/without/neutrinos/DplusS.txt')
LAMplusC = np.genfromtxt('code/tabulate/magnetar/without/neutrinos/LAMplusC.txt')

pi = pi[i, :]
K = K[i, :]
D0 = D0[i, :]
Dplus = Dplus[i, :]
DplusS = DplusS[i, :]
LAMplusC = LAMplusC[i, :]

c = D0 + Dplus + DplusS + LAMplusC

N = c.max()

plt.plot(t, c / N, 'b-', label=r'Total Charm Decay', zorder=2)
plt.plot(t, D0 / N, 'b--', label=r'$D^0$ Decay', zorder=2)
plt.plot(t, Dplus / N, 'r-', label=r'$D^+$ Decay', zorder=1)
plt.plot(t, DplusS / N, 'r--', label=r'$D^+_s$ Decay', zorder=1)
plt.plot(t, LAMplusC / N, 'k-', label=r'$\Lambda^{\kern-0.5pt +}_{\kern+0.5pt c}$ Decay', zorder=0)

plt.xlabel(r'$t$ $\mathrel{/} \kern-0.8pt$ s')
plt.ylabel(r'$\dot{\phi}_\nu$ $\mathrel{/} \kern-0.7pt$ max$\bigl( \dot{\phi}^c_\nu \kern+0.2pt \bigr)$')

plt.xscale('log')
plt.yscale('log')

plt.xlim(2e2, 1e6)
plt.ylim(1e-4, 2e0)

plt.legend(handlelength=1.8, loc=2)

plt.savefig('build/magnetar_charm_decay_comparison_without.pdf')
plt.savefig('build/magnetar_charm_decay_comparison_without.png')
plt.close()

plt.fill_between(t, 3 * c / N, c / (3 * N), color='none', facecolor='b', alpha=0.25, zorder=0)

plt.plot(t, pi / N, 'r', label=r'Pion Decay', zorder=1)
plt.plot(t, K / N, 'k', label=r'Kaon Decay', zorder=1)
plt.plot(t, c / N, 'b', label=r'Charm Decay', zorder=0)

plt.xlabel(r'$t$ $\mathrel{/} \kern-0.8pt$ s')
plt.ylabel(r'$\dot{\phi}_\nu$ $\mathrel{/} \kern-0.7pt$ max$\bigl( \dot{\phi}^c_\nu \kern+0.2pt \bigr)$')

plt.xscale('log')
plt.yscale('log')

plt.xlim(4e1, 4e6)
plt.ylim(1e-6, 1e1)

plt.legend(handlelength=1.8, loc=2)

plt.savefig('build/magnetar_neutrino_spectrum_without.pdf')
plt.savefig('build/magnetar_neutrino_spectrum_without.png')
plt.close()



E, pi1, pi2, pi3 = np.genfromtxt('code/tabulate/magnetar/without/integrate/pi.txt', unpack=True)
E, K1, K2, K3 = np.genfromtxt('code/tabulate/magnetar/without/integrate/K.txt', unpack=True)
E, D01, D02, D03 = np.genfromtxt('code/tabulate/magnetar/without/integrate/D0.txt', unpack=True)
E, Dplus1, Dplus2, Dplus3 = np.genfromtxt('code/tabulate/magnetar/without/integrate/Dplus.txt', unpack=True)
E, DplusS1, DplusS2, DplusS3 = np.genfromtxt('code/tabulate/magnetar/without/integrate/DplusS.txt', unpack=True)
E, LAMplusC1, LAMplusC2, LAMplusC3 = np.genfromtxt('code/tabulate/magnetar/without/integrate/LAMplusC.txt', unpack=True)

c1 = D01 + Dplus1 + DplusS1 + LAMplusC1
c2 = D02 + Dplus2 + DplusS2 + LAMplusC2
c3 = D03 + Dplus3 + DplusS3 + LAMplusC3

N = (E**2 * c3).max()

plt.scatter([], [], color='r', label=r'Pion Decay')
plt.scatter([], [], color='k', label=r'Kaon Decay')
plt.scatter([], [], color='b', label=r'Charm Decay')

plt.fill_between(E, 3 * E**2 * c3 / N, E**2 * c3 / (3 * N), color='none', facecolor='b', alpha=0.25)

plt.plot(E, E**2 * c1 / N, 'b:')
plt.plot(E, E**2 * c2 / N, 'b--')
plt.plot(E, E**2 * c3 / N, 'b-')
plt.plot(E, E**2 * pi1 / N, 'r:')
plt.plot(E, E**2 * pi2 / N, 'r--')
plt.plot(E, E**2 * pi3 / N, 'r-')
plt.plot(E, E**2 * K1 / N, 'k:', label=r'$10^3 \kern+1.5pt$s $\kern+0.3pt -$ $10^4 \kern+1.5pt$s')
plt.plot(E, E**2 * K2 / N, 'k--', label=r'$10^4 \kern+1.5pt$s $\kern+0.3pt -$ $10^5 \kern+1.5pt$s')
plt.plot(E, E**2 * K3 / N, 'k-', label=r'$10^3 \kern+1.5pt$s $\kern+0.3pt -$ $10^7 \kern+1.5pt$s') 

plt.xlabel(r'$E_\nu$ $\mathrel{/} \kern-0.6pt$ GeV')
plt.ylabel(r'$E_\nu^2\phi_\nu$ $\mathrel{/} \kern-0.7pt$ max$\bigl( \kern-0.3pt E_\nu^2\phi^c_\nu \kern+0.2pt \bigr)$')

plt.xscale('log')
plt.yscale('log')

plt.xlim(3e5, 1e11)
plt.ylim(1.6e-4, 3e3)

plt.legend(handlelength=1.8, loc=1)

plt.savefig('build/magnetar_integrated_neutrino_spectrum_without.pdf')
plt.savefig('build/magnetar_integrated_neutrino_spectrum_without.png')
plt.close()



# t = np.genfromtxt('code/tabulate/magnetar/with/neutrinos/axes.txt', skip_footer=1)
# E = np.genfromtxt('code/tabulate/magnetar/with/neutrinos/axes.txt', skip_header=15)
# 
# en = 1e9
# i = (np.abs(E - en)).argmin()
# en = E[i]
# 
# pi = np.genfromtxt('code/tabulate/magnetar/with/neutrinos/pi.txt')
# K = np.genfromtxt('code/tabulate/magnetar/with/neutrinos/K.txt')
# D0 = np.genfromtxt('code/tabulate/magnetar/with/neutrinos/D0.txt')
# Dplus = np.genfromtxt('code/tabulate/magnetar/with/neutrinos/Dplus.txt')
# DplusS = np.genfromtxt('code/tabulate/magnetar/with/neutrinos/DplusS.txt')
# LAMplusC = np.genfromtxt('code/tabulate/magnetar/with/neutrinos/LAMplusC.txt')
# 
# pi = pi[i, :]
# K = K[i, :]
# D0 = D0[i, :]
# Dplus = Dplus[i, :]
# DplusS = DplusS[i, :]
# LAMplusC = LAMplusC[i, :]
# 
# c = D0 + Dplus + DplusS + LAMplusC
# 
# N = c.max()
# 
# plt.plot(t, c / N, 'b-', label=r'Total Charm Decay', zorder=2)
# plt.plot(t, D0 / N, 'b--', label=r'$D^0$ Decay', zorder=2)
# plt.plot(t, Dplus / N, 'r-', label=r'$D^+$ Decay', zorder=1)
# plt.plot(t, DplusS / N, 'r--', label=r'$D^+_s$ Decay', zorder=1)
# plt.plot(t, LAMplusC / N, 'k-', label=r'$\Lambda^{\kern-0.5pt +}_{\kern+0.5pt c}$ Decay', zorder=0)
# 
# plt.xlabel(r'$t$ $\mathrel{/} \kern-0.8pt$ s')
# plt.ylabel(r'$\dot{\phi}_\nu$ $\mathrel{/} \kern-0.7pt$ max$\bigl( \dot{\phi}^c_\nu \kern+0.2pt \bigr)$')
# 
# plt.xscale('log')
# plt.yscale('log')
# 
# plt.xlim(2e2, 1e6)
# plt.ylim(1e-4, 2e0)
# 
# plt.legend(handlelength=1.8, loc=3)
# 
# plt.savefig('build/magnetar_charm_decay_comparison_with.pdf')
# plt.savefig('build/magnetar_charm_decay_comparison_with.png')
# plt.close()
# 
# plt.fill_between(t, 3 * c / N, c / (3 * N), color='none', facecolor='b', alpha=0.25, zorder=0)
# 
# plt.plot(t, pi / N, 'r', label=r'Pion Decay', zorder=1)
# plt.plot(t, K / N, 'k', label=r'Kaon Decay', zorder=1)
# plt.plot(t, c / N, 'b', label=r'Charm Decay', zorder=0)
# 
# plt.xlabel(r'$t$ $\mathrel{/} \kern-0.8pt$ s')
# plt.ylabel(r'$\dot{\phi}_\nu$ $\mathrel{/} \kern-0.7pt$ max$\bigl( \dot{\phi}^c_\nu \kern+0.2pt \bigr)$')
# 
# plt.xscale('log')
# plt.yscale('log')
# 
# plt.xlim(4e1, 4e6)
# plt.ylim(1e-6, 1e1)
# 
# plt.legend(handlelength=1.8, loc=3)
# 
# plt.savefig('build/magnetar_neutrino_spectrum_with.pdf')
# plt.savefig('build/magnetar_neutrino_spectrum_with.png')
# plt.close()
# 
# 
# 
# E, pi1, pi2, pi3 = np.genfromtxt('code/tabulate/magnetar/with/integrate/pi.txt', unpack=True)
# E, K1, K2, K3 = np.genfromtxt('code/tabulate/magnetar/with/integrate/K.txt', unpack=True)
# E, D01, D02, D03 = np.genfromtxt('code/tabulate/magnetar/with/integrate/D0.txt', unpack=True)
# E, Dplus1, Dplus2, Dplus3 = np.genfromtxt('code/tabulate/magnetar/with/integrate/Dplus.txt', unpack=True)
# E, DplusS1, DplusS2, DplusS3 = np.genfromtxt('code/tabulate/magnetar/with/integrate/DplusS.txt', unpack=True)
# E, LAMplusC1, LAMplusC2, LAMplusC3 = np.genfromtxt('code/tabulate/magnetar/with/integrate/LAMplusC.txt', unpack=True)
# 
# c1 = D01 + Dplus1 + DplusS1 + LAMplusC1
# c2 = D02 + Dplus2 + DplusS2 + LAMplusC2
# c3 = D03 + Dplus3 + DplusS3 + LAMplusC3
# 
# N = (E**2 * c3).max()
# 
# plt.scatter([], [], color='r', label=r'Pion Decay')
# plt.scatter([], [], color='k', label=r'Kaon Decay')
# plt.scatter([], [], color='b', label=r'Charm Decay')
# 
# plt.fill_between(E, 3 * E**2 * c3 / N, E**2 * c3 / (3 * N), color='none', facecolor='b', alpha=0.25)
# 
# plt.plot(E, E**2 * c1 / N, 'b:')
# plt.plot(E, E**2 * c2 / N, 'b--')
# plt.plot(E, E**2 * c3 / N, 'b-')
# plt.plot(E, E**2 * pi1 / N, 'r:')
# plt.plot(E, E**2 * pi2 / N, 'r--')
# plt.plot(E, E**2 * pi3 / N, 'r-')
# plt.plot(E, E**2 * K1 / N, 'k:', label=r'$10^3 \kern+1.5pt$s $\kern+0.3pt -$ $10^4 \kern+1.5pt$s')
# plt.plot(E, E**2 * K2 / N, 'k--', label=r'$10^4 \kern+1.5pt$s $\kern+0.3pt -$ $10^5 \kern+1.5pt$s')
# plt.plot(E, E**2 * K3 / N, 'k-', label=r'$10^3 \kern+1.5pt$s $\kern+0.3pt -$ $10^7 \kern+1.5pt$s') 
# 
# plt.xlabel(r'$E_\nu$ $\mathrel{/} \kern-0.6pt$ GeV')
# plt.ylabel(r'$E_\nu^2\phi_\nu$ $\mathrel{/} \kern-0.7pt$ max$\bigl( \kern-0.3pt E_\nu^2\phi^c_\nu \kern+0.2pt \bigr)$')
# 
# plt.xscale('log')
# plt.yscale('log')
# 
# plt.xlim(3e5, 1e11)
# plt.ylim(4.3e-4, 6e1)
# 
# plt.legend(handlelength=1.8, loc=1)
# 
# plt.savefig('build/magnetar_integrated_neutrino_spectrum_with.pdf')
# plt.savefig('build/magnetar_integrated_neutrino_spectrum_with.png')
# plt.close()
# 
# 
# 
# E, pi = np.genfromtxt('code/tabulate/nucleus/neutrinos/pi.txt', unpack=True)
# E, K = np.genfromtxt('code/tabulate/nucleus/neutrinos/K.txt', unpack=True)
# E, D0 = np.genfromtxt('code/tabulate/nucleus/neutrinos/D0.txt', unpack=True)
# E, Dplus = np.genfromtxt('code/tabulate/nucleus/neutrinos/Dplus.txt', unpack=True)
# E, DplusS = np.genfromtxt('code/tabulate/nucleus/neutrinos/DplusS.txt', unpack=True)
# E, LAMplusC = np.genfromtxt('code/tabulate/nucleus/neutrinos/LAMplusC.txt', unpack=True)
# 
# c = D0 + Dplus + DplusS + LAMplusC
# 
# N = (E**2 * c).max()
# 
# plt.plot(E, E**2 * c / N, 'b-', label=r'Total Charm Decay', zorder=2)
# plt.plot(E, E**2 * D0 / N, 'b--', label=r'$D^0$ Decay', zorder=2)
# plt.plot(E, E**2 * Dplus / N, 'r-', label=r'$D^+$ Decay', zorder=1)
# plt.plot(E, E**2 * DplusS / N, 'r--', label=r'$D^+_s$ Decay', zorder=1)
# plt.plot(E, E**2 * LAMplusC / N, 'k-', label=r'$\Lambda^{\kern-0.5pt +}_{\kern+0.5pt c}$ Decay', zorder=0)
# 
# plt.xlabel(r'$E_\nu$ $\mathrel{/} \kern-0.6pt$ GeV')
# plt.ylabel(r'$E_\nu^2\phi_\nu$ $\mathrel{/} \kern-0.7pt$ max$\bigl( \kern-0.3pt E_\nu^2\phi^c_\nu \kern+0.2pt \bigr)$')
# 
# plt.xscale('log')
# plt.yscale('log')
# 
# plt.xlim(3e5, 1e11)
# plt.ylim(5e-4, 2e0)
# 
# plt.legend(handlelength=1.8, loc=3)
# 
# plt.savefig('build/nucleus_charm_decay_comparison.pdf')
# plt.savefig('build/nucleus_charm_decay_comparison.png')
# plt.close()
# 
# plt.fill_between(E, 3 * E**2 * c / N, E**2 * c / (3 * N), color='none', facecolor='b', alpha=0.25)
# 
# plt.plot(E, E**2 * pi / N, 'r', label=r'Pion Decay', zorder=1)
# plt.plot(E, E**2 * K / N, 'k', label=r'Kaon Decay', zorder=1)
# plt.plot(E, E**2 * c / N, 'b', label=r'Charm Decay', zorder=0) 
# 
# plt.xlabel(r'$E_\nu$ $\mathrel{/} \kern-0.6pt$ GeV')
# plt.ylabel(r'$E_\nu^2\phi_\nu$ $\mathrel{/} \kern-0.7pt$ max$\bigl( \kern-0.3pt E_\nu^2\phi^c_\nu \kern+0.2pt \bigr)$')
# 
# plt.xscale('log')
# plt.yscale('log')
# 
# plt.xlim(3e5, 1e11)
# plt.ylim(5e-4, 2e3)
# 
# plt.legend(handlelength=1.8, loc=1)
# 
# plt.savefig('build/nucleus_neutrino_spectrum.pdf')
# plt.savefig('build/nucleus_neutrino_spectrum.png')
# plt.close()
