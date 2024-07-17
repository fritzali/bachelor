import numpy as np
import matplotlib.pyplot as plt

from functional import *



x = np.logspace(-7, 0, 1000)
y12 = charmed_hadron_differential_production(x, 1e12, 'd0')
y10 = charmed_hadron_differential_production(x, 1e10, 'd0')
y8 = charmed_hadron_differential_production(x, 1e8, 'd0')

plt.plot(x, x * y12, 'b', label=r'$E_p = 10^{12}$ GeV')
plt.plot(x, x * y10, 'r', label=r'$E_p = 10^{10}$ GeV')
plt.plot(x, x * y8, 'k', label=r'$E_p = 10^{8}$ GeV')

plt.xlabel(r'$x_h \kern-0.8pt = \kern-0.5pt E_h \kern+0.5pt / E_p$')
plt.ylabel(r'$x_h \kern+0.5pt d\sigma \kern-0.3pt / \kern-0.8pt dx_h \kern+0.2pt$ $\mathrel{/}$ mb')

plt.xscale('log')
plt.yscale('log')

plt.xlim(1e-6, 1e0)
plt.ylim(1e-4, 1e1)

plt.legend()

plt.savefig('build/charm_hadron_cross_section.pdf')
plt.savefig('build/charm_hadron_cross_section.png')
plt.close()



s = np.logspace(1, 10, 1000)
yp = inelastic_hadron_proton_scattering(s, 'p')
ypi = inelastic_hadron_proton_scattering(s, 'pi')
yk = inelastic_hadron_proton_scattering(s, 'k')

plt.plot(s, yp, 'b', label=r'$p \kern-0.1pt p$')
plt.plot(s, ypi, 'r', label=r'$\pi \kern-0.1pt p$')
plt.plot(s, yk, 'k', label=r'$K \kern-0.2pt p$')

plt.xlabel(r'$s$ $\mathrel{/} \kern-0.1pt$ GeV$^2$')
plt.ylabel(r'$\sigma_{h \kern-0.1pt p}$ $\mathrel{/}$ mb')

plt.xscale('log')

plt.xlim(1e2, 1e9)

plt.legend()

plt.savefig('build/hadron_proton_scattering.pdf')
plt.savefig('build/hadron_proton_scattering.png')
plt.close()



t = np.genfromtxt('code/tabulate/magnetar/neutrinos/axes.txt', skip_footer=1)
E = np.genfromtxt('code/tabulate/magnetar/neutrinos/axes.txt', skip_header=15)

en = 1e9
i = (np.abs(E - en)).argmin()
en = E[i]

pi = np.genfromtxt('code/tabulate/magnetar/neutrinos/pi.txt')
K = np.genfromtxt('code/tabulate/magnetar/neutrinos/K.txt')
D0 = np.genfromtxt('code/tabulate/magnetar/neutrinos/D0.txt')
Dplus = np.genfromtxt('code/tabulate/magnetar/neutrinos/Dplus.txt')
DplusS = np.genfromtxt('code/tabulate/magnetar/neutrinos/DplusS.txt')
LAMplusC = np.genfromtxt('code/tabulate/magnetar/neutrinos/LAMplusC.txt')

pi = pi[i, :]
K = K[i, :]
D0 = D0[i, :]
Dplus = Dplus[i, :]
DplusS = DplusS[i, :]
LAMplusC = LAMplusC[i, :]

c = D0 + Dplus + DplusS + LAMplusC

N = c.max()

plt.plot(t, c / N, 'b-', label=r'Total Charm Decay', zorder=5)
plt.plot(t, D0 / N, 'b--', label=r'$D^0$ Decay', zorder=4)
plt.plot(t, Dplus / N, 'r-', label=r'$D^+$ Decay', zorder=3)
plt.plot(t, DplusS / N, 'r--', label=r'$D^+_s$ Decay', zorder=2)
plt.plot(t, LAMplusC / N, 'k-', label=r'$\Lambda^{\kern-0.5pt +}_{\kern+0.5pt c}$ Decay', zorder=1)

plt.xlabel(r'$t$ $\mathrel{/} \kern-0.8pt$ s')
plt.ylabel(r'$\dot{\phi}_\nu$ $\mathrel{/} \kern-0.7pt$ max$\bigl( \dot{\phi}^c_\nu \kern+0.2pt \bigr)$')

plt.xscale('log')
plt.yscale('log')

plt.xlim(2e2, 1e6)
plt.ylim(1e-4, 2e0)

plt.legend(handlelength=1.8, loc=2)

plt.savefig('build/charm_decay_comparison.pdf')
plt.savefig('build/charm_decay_comparison.png')
plt.close()

plt.fill_between(t, 3 * c / N, c / (3 * N), color='b', alpha=0.25, zorder=0)

plt.plot(t, pi / N, 'r', label=r'Pion Decay', zorder=3)
plt.plot(t, K / N, 'k', label=r'Kaon Decay', zorder=2)
plt.plot(t, c / N, 'b', label=r'Charm Decay', zorder=1)

plt.xlabel(r'$t$ $\mathrel{/} \kern-0.8pt$ s')
plt.ylabel(r'$\dot{\phi}_\nu$ $\mathrel{/} \kern-0.7pt$ max$\bigl( \dot{\phi}^c_\nu \kern+0.2pt \bigr)$')

plt.xscale('log')
plt.yscale('log')

plt.xlim(4e1, 4e6)
plt.ylim(1e-6, 1e1)

plt.legend(handlelength=1.8, loc=2)

plt.savefig('build/neutrino_spectrum.pdf')
plt.savefig('build/neutrino_spectrum.png')
plt.close()



E, pi1, pi2, pi3 = np.genfromtxt('code/tabulate/magnetar/integrate/pi.txt', unpack=True)
E, K1, K2, K3 = np.genfromtxt('code/tabulate/magnetar/integrate/K.txt', unpack=True)
E, D01, D02, D03 = np.genfromtxt('code/tabulate/magnetar/integrate/D0.txt', unpack=True)
E, Dplus1, Dplus2, Dplus3 = np.genfromtxt('code/tabulate/magnetar/integrate/Dplus.txt', unpack=True)
E, DplusS1, DplusS2, DplusS3 = np.genfromtxt('code/tabulate/magnetar/integrate/DplusS.txt', unpack=True)
E, LAMplusC1, LAMplusC2, LAMplusC3 = np.genfromtxt('code/tabulate/magnetar/integrate/LAMplusC.txt', unpack=True)

c1 = D01 + Dplus1 + DplusS1 + LAMplusC1
c2 = D02 + Dplus2 + DplusS2 + LAMplusC2
c3 = D03 + Dplus3 + DplusS3 + LAMplusC3

N = (E**2 * c3).max()

# plt.plot(E, E**2 * pi1 / N, 'r:')
# plt.plot(E, E**2 * pi2 / N, 'r--')
plt.plot(E, E**2 * pi3 / N, 'r-', label=r'Pion Decay')
# plt.plot(E, E**2 * K1 / N, 'k:')
# plt.plot(E, E**2 * K2 / N, 'k--')
plt.plot(E, E**2 * K3 / N, 'k-', label=r'Kaon Decay')
# plt.plot(E, E**2 * c1 / N, 'b:')
# plt.plot(E, E**2 * c2 / N, 'b--')
plt.plot(E, E**2 * c3 / N, 'b-', label=r'Charm Decay')

plt.xlabel(r'$t$ $\mathrel{/} \kern-0.8pt$ s')
plt.ylabel(r'$E_\nu^2\phi_\nu$ $\mathrel{/} \kern-0.7pt$ max$\bigl( \kern-0.3pt E_\nu^2\phi^c_\nu \kern+0.2pt \bigr)$')

plt.xscale('log')
plt.yscale('log')

plt.xlim(3e5, 1e11)
plt.ylim(1e-4, 3e3)

plt.legend(handlelength=1.8, loc=1)

plt.savefig('build/integrated_neutrino_spectrum.pdf')
plt.savefig('build/integrated_neutrino_spectrum.png')
plt.close()
