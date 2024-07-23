import numpy as np
import matplotlib.pyplot as plt

from code.functional import *


plt.figure(figsize=(5.4, 3.1))

h, y = np.genfromtxt('code/tabulate/other/events.txt', unpack=True, dtype=None)

i = [3, 5, 0, 2, 8, 7, 1, 4, 10, 9, 6, 11, 14, 15, 16, 20, 21, 17, 18, 19, 13, 12, 22, 24, 25, 23]

plt.bar(range(len(y)), y[i], tick_label=h[i], color='b', linewidth=0, width=0.6, alpha=0.6, label=r'\textsc{sibyll} 2.3c')

plt.ylabel(r'yield $\mathrel{/} \kern-0.25pt$ event')

plt.yscale('log')

plt.gca().tick_params(axis='x', which='minor', bottom=False)

plt.xlim(-0.7, 25.7)

plt.legend(loc=1)

plt.savefig('build/event_generator.pdf')
plt.savefig('build/event_generator.png')
plt.close()



s, p, pi, K = np.genfromtxt('code/tabulate/other/sample_inelastic_scattering.txt', unpack=True)

plt.figure(figsize=(5.0, 3.2))

plt.plot(s, p, 'b', label=r'$p \kern-0.1pt p$')
plt.plot(s, pi, 'r', label=r'$\pi \kern-0.1pt p$')
plt.plot(s, K, 'k', label=r'$K \kern-0.2pt p$')

plt.xlabel(r'$s$ $\mathrel{/} \kern-0.1pt$ GeV$^2$')
plt.ylabel(r'$\sigma_{h \kern-0.1pt p}$ $\mathrel{/}$ mb')

plt.xscale('log')

plt.xlim(1e2, 1e9)

plt.legend(handlelength=1.5, loc=2)

plt.savefig('build/hadron_proton_scattering.pdf')
plt.savefig('build/hadron_proton_scattering.png')
plt.close()



x, y1, y2, y3 = np.genfromtxt('code/tabulate/other/sample_charm_hadron.txt', unpack=True)

plt.figure(figsize=(5.0, 3.2))

plt.plot(x, x * y1, 'b', label=r'$E_p = 10^{12}$ GeV')
plt.plot(x, x * y2, 'r', label=r'$E_p = 10^{10}$ GeV')
plt.plot(x, x * y3, 'k', label=r'$E_p = 10^{8}$ GeV')

plt.xlabel(r'$x_h \kern-0.8pt = \kern-0.5pt E_h \kern+0.5pt / E_p$')
plt.ylabel(r'$x_h \kern+0.5pt d\sigma \kern-0.3pt / \kern-0.8pt dx_h \kern+0.2pt$ $\mathrel{/}$ mb')

plt.xscale('log')
plt.yscale('log')

plt.xlim(1e-6, 1e0)
plt.ylim(1e-4, 1e1)

plt.legend(handlelength=1.5, loc=1)

plt.savefig('build/charm_hadron_cross_section.pdf')
plt.savefig('build/charm_hadron_cross_section.png')
plt.close()



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

plt.legend(handlelength=1.5, loc=2)

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

plt.legend(handlelength=1.5, loc=2)

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

plt.legend(handlelength=1.5, loc=1)

plt.savefig('build/magnetar_integrated_neutrino_spectrum_without.pdf')
plt.savefig('build/magnetar_integrated_neutrino_spectrum_without.png')
plt.close()



t = np.genfromtxt('code/tabulate/magnetar/with/neutrinos/axes.txt', skip_footer=1)
E = np.genfromtxt('code/tabulate/magnetar/with/neutrinos/axes.txt', skip_header=15)

en = 1e9
i = (np.abs(E - en)).argmin()
en = E[i]

pi = np.genfromtxt('code/tabulate/magnetar/with/neutrinos/pi.txt')
K = np.genfromtxt('code/tabulate/magnetar/with/neutrinos/K.txt')
D0 = np.genfromtxt('code/tabulate/magnetar/with/neutrinos/D0.txt')
Dplus = np.genfromtxt('code/tabulate/magnetar/with/neutrinos/Dplus.txt')
DplusS = np.genfromtxt('code/tabulate/magnetar/with/neutrinos/DplusS.txt')
LAMplusC = np.genfromtxt('code/tabulate/magnetar/with/neutrinos/LAMplusC.txt')

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

plt.legend(handlelength=1.5, loc=3)

plt.savefig('build/magnetar_charm_decay_comparison_with.pdf')
plt.savefig('build/magnetar_charm_decay_comparison_with.png')
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

plt.legend(handlelength=1.5, loc=3)

plt.savefig('build/magnetar_neutrino_spectrum_with.pdf')
plt.savefig('build/magnetar_neutrino_spectrum_with.png')
plt.close()



E, pi1, pi2, pi3 = np.genfromtxt('code/tabulate/magnetar/with/integrate/pi.txt', unpack=True)
E, K1, K2, K3 = np.genfromtxt('code/tabulate/magnetar/with/integrate/K.txt', unpack=True)
E, D01, D02, D03 = np.genfromtxt('code/tabulate/magnetar/with/integrate/D0.txt', unpack=True)
E, Dplus1, Dplus2, Dplus3 = np.genfromtxt('code/tabulate/magnetar/with/integrate/Dplus.txt', unpack=True)
E, DplusS1, DplusS2, DplusS3 = np.genfromtxt('code/tabulate/magnetar/with/integrate/DplusS.txt', unpack=True)
E, LAMplusC1, LAMplusC2, LAMplusC3 = np.genfromtxt('code/tabulate/magnetar/with/integrate/LAMplusC.txt', unpack=True)

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
plt.ylim(4.3e-4, 6e1)

plt.legend(handlelength=1.5, loc=1)

plt.savefig('build/magnetar_integrated_neutrino_spectrum_with.pdf')
plt.savefig('build/magnetar_integrated_neutrino_spectrum_with.png')
plt.close()



E, pi = np.genfromtxt('code/tabulate/nucleus/neutrinos/pi.txt', unpack=True)
E, K = np.genfromtxt('code/tabulate/nucleus/neutrinos/K.txt', unpack=True)
E, D0 = np.genfromtxt('code/tabulate/nucleus/neutrinos/D0.txt', unpack=True)
E, Dplus = np.genfromtxt('code/tabulate/nucleus/neutrinos/Dplus.txt', unpack=True)
E, DplusS = np.genfromtxt('code/tabulate/nucleus/neutrinos/DplusS.txt', unpack=True)
E, LAMplusC = np.genfromtxt('code/tabulate/nucleus/neutrinos/LAMplusC.txt', unpack=True)

c = D0 + Dplus + DplusS + LAMplusC

N = (E**2 * c).max()

plt.plot(E, E**2 * c / N, 'b-', label=r'Total Charm Decay', zorder=2)
plt.plot(E, E**2 * D0 / N, 'b--', label=r'$D^0$ Decay', zorder=2)
plt.plot(E, E**2 * Dplus / N, 'r-', label=r'$D^+$ Decay', zorder=1)
plt.plot(E, E**2 * DplusS / N, 'r--', label=r'$D^+_s$ Decay', zorder=1)
plt.plot(E, E**2 * LAMplusC / N, 'k-', label=r'$\Lambda^{\kern-0.5pt +}_{\kern+0.5pt c}$ Decay', zorder=0)

plt.xlabel(r'$E_\nu$ $\mathrel{/} \kern-0.6pt$ GeV')
plt.ylabel(r'$E_\nu^2\phi_\nu$ $\mathrel{/} \kern-0.7pt$ max$\bigl( \kern-0.3pt E_\nu^2\phi^c_\nu \kern+0.2pt \bigr)$')

plt.xscale('log')
plt.yscale('log')

plt.xlim(3e5, 1e11)
plt.ylim(5e-4, 2e0)

plt.legend(handlelength=1.5, loc=3)

plt.savefig('build/nucleus_charm_decay_comparison.pdf')
plt.savefig('build/nucleus_charm_decay_comparison.png')
plt.close()

plt.fill_between(E, 3 * E**2 * c / N, E**2 * c / (3 * N), color='none', facecolor='b', alpha=0.25)

plt.plot(E, E**2 * pi / N, 'r', label=r'Pion Decay', zorder=1)
plt.plot(E, E**2 * K / N, 'k', label=r'Kaon Decay', zorder=1)
plt.plot(E, E**2 * c / N, 'b', label=r'Charm Decay', zorder=0) 

plt.xlabel(r'$E_\nu$ $\mathrel{/} \kern-0.6pt$ GeV')
plt.ylabel(r'$E_\nu^2\phi_\nu$ $\mathrel{/} \kern-0.7pt$ max$\bigl( \kern-0.3pt E_\nu^2\phi^c_\nu \kern+0.2pt \bigr)$')

plt.xscale('log')
plt.yscale('log')

plt.xlim(3e5, 1e11)
plt.ylim(5e-4, 2e3)

plt.legend(handlelength=1.5, loc=1)

plt.savefig('build/nucleus_neutrino_spectrum.pdf')
plt.savefig('build/nucleus_neutrino_spectrum.png')
plt.close()
