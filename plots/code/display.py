import numpy as np
import matplotlib.pyplot as plt

from code.functional import *


RUB = '#17365C'
TUDO = '#83B818'
ACC = '#D98207'



k_x, k_sig = np.genfromtxt('code/tabulate/other/present_k_differential.txt')
pi_x, pi_sig = np.genfromtxt('code/tabulate/other/present_pi_differential.txt')

E = 10154863318.74878
pp_sig_inel = inelastic_hadron_proton_scattering(E, 'p')
x = np.logspace(np.log10(min(k_x)), np.log10(max(pi_x)), 1000)

k_cr = E * pp_sig_inel * meson_production(x, E, 'k')
pi_cr = E * pp_sig_inel * meson_production(x, E, 'pi')

plt.plot(k_x, k_x * k_sig, c=RUB, label=r'Schroller')
plt.plot(x, x * k_cr, c=TUDO, label=r'Kelner{\kern+0.5pt}*')

plt.xscale('log')
plt.yscale('log')

plt.xlim(1e-11, 1e0)
plt.ylim(1e-1, plt.ylim()[1])

plt.xticks(np.logspace(-10, 0, 6))

plt.xlabel(r'$x_E = \kern+0.1pt E_K \kern+0.4pt / E_p$')
plt.ylabel(r'$x_E \kern+0.4pt d \kern-0.3pt \sigma / d \kern-0.3pt x_E \kern+0.1pt$ $\mathrel{/}$ $\symup{mb}$')

plt.title(r'$pp \rightarrow K^+ X$ $( E_p \kern-0.3pt = \kern-0.2pt 10^{10}$ $\symup{GeV} \kern+0.7pt )$', pad=7)

plt.legend(loc=2)

plt.savefig('build/present_kaon_weighted.pdf')
plt.savefig('build/present_kaon_weighted.png')
plt.close()

plt.plot(pi_x, pi_x * pi_sig, c=RUB, label=r'Schroller')
plt.plot(x, x * pi_cr, c=TUDO, label=r'Kelner{\kern+0.5pt}*')

plt.xscale('log')
plt.yscale('log')

plt.xlim(1e-11, 1e0)
plt.ylim(1e0, plt.ylim()[1])

plt.xticks(np.logspace(-10, 0, 6))

plt.xlabel(r'$x_E = \kern+0.1pt E_\pi \kern+0.4pt / E_p$')
plt.ylabel(r'$x_E \kern+0.4pt d \kern-0.3pt \sigma / d \kern-0.3pt x_E \kern+0.1pt$ $\mathrel{/}$ $\symup{mb}$')

plt.title(r'$pp \rightarrow \pi^+ X$ $( E_p \kern-0.3pt = \kern-0.2pt 10^{10}$ $\symup{GeV} \kern+0.7pt )$', pad=7)

plt.legend(loc=2)

plt.savefig('build/present_pion_weighted.pdf')
plt.savefig('build/present_pion_weighted.png')
plt.close()



plt.figure(figsize=(5.5, 3.5))

h, y = np.genfromtxt('code/tabulate/other/present_events.txt', unpack=True, dtype=None)

i = [3, 5, 0, 2, 8, 7, 1, 4, 10, 9, 6, 11, 14, 15, 16, 20, 21, 17, 18, 19, 13, 12, 22, 24, 25, 23]

plt.bar(range(len(y)), y[i], tick_label=h[i], color=TUDO, linewidth=0, width=0.75, label=r'\textsc{sibyll} 2.3c')

plt.ylabel(r'yield / \kern-0.15pt event')

plt.yscale('log')

plt.gca().tick_params(axis='x', which='minor', bottom=False)

plt.xlim(-0.625, 25.625)

plt.legend(loc=1)

plt.savefig('build/present_event_generator.pdf')
plt.savefig('build/present_event_generator.png')
plt.close()



plt.figure(figsize=(5.5, 3.5))

s, p, pi, K = np.genfromtxt('code/tabulate/other/sample_inelastic_scattering.txt', unpack=True)

plt.plot(s, p, label=r'$p p$', c=TUDO)
plt.plot(s, pi, label=r'$\pi^+ \kern-0.5pt p$', c=ACC)
plt.plot(s, K, label=r'$K^+ \kern-0.5pt p$', c='k')

plt.xlabel(r'$s$ $\mathrel{/} \kern-0.1pt$ $\symup{GeV}^2$')
plt.ylabel(r'$\sigma$ $\mathrel{/}$ $\symup{mb}$')

plt.xscale('log')

plt.xlim(1e2, 1e9)

plt.legend(loc=2)

plt.savefig('build/present_hadron_scattering.pdf')
plt.savefig('build/present_hadron_scattering.png')
plt.close()



plt.figure(figsize=(5.7, 3.8))

x, y1, y2, y3 = np.genfromtxt('code/tabulate/other/sample_charm_hadron.txt', unpack=True)

E = 1e10
pp_sig_inel = inelastic_hadron_proton_scattering(E, 'p')
x = np.logspace(-6.1, -0.001, 1000)

pi = E * pp_sig_inel * meson_production(x, E, 'pi')
K = E * pp_sig_inel * meson_production(x, E, 'k')
D = charmed_hadron_differential_production(x, E, 'd0')

plt.plot(x, x * pi, label=r'$h = \pi^+$', c=ACC)
plt.plot(x, x * K, label=r'$h = K^+$', c='k')
plt.plot(x, x * D, label=r'$h = D^0$', c=TUDO)

plt.xlabel(r'$x_E = \kern+0.1pt E_h \kern+0.4pt / E_p$')
plt.ylabel(r'$x_E \kern+0.4pt d \kern-0.3pt \sigma / d \kern-0.3pt x_E \kern+0.1pt$ $\mathrel{/}$ $\symup{mb}$')

plt.title(r'$pp \rightarrow hX$ $( E_p \kern-0.3pt = \kern-0.2pt 10^{10}$ $\symup{GeV} \kern+0.7pt )$', pad=7)

plt.xscale('log')
plt.yscale('log')

plt.xlim(1e-6, 1e0)
plt.ylim(1e-4, 1e3)

plt.legend(loc=3)

plt.savefig('build/present_hadron_weighted.pdf')
plt.savefig('build/present_hadron_weighted.png')
plt.close()



plt.figure(figsize=(5.5, 3.5))

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

plt.plot(t, c / N, '-', label=r'Total Charm Decay', zorder=2, c=TUDO)
plt.plot(t, D0 / N, '--', label=r'$D^0$ Decay', zorder=2, c=TUDO)
plt.plot(t, Dplus / N, '-', label=r'$D^+$ Decay', zorder=1, c=ACC)
plt.plot(t, DplusS / N, '--', label=r'$D^+_s$ Decay', zorder=1, c=ACC)
plt.plot(t, LAMplusC / N, 'k-', label=r'$\Lambda^{\kern-0.5pt +}_{\kern+0.5pt c}$ Decay', zorder=0)

plt.xlabel(r'$t$ $\mathrel{/}$ $\symup{s}$')
plt.ylabel(r'$\dot{\phi}_{\kern-0.3pt \nu}$ $\mathrel{/} \kern-0.7pt$ $\symup{max} \kern+0.5pt \bigl( \dot{\phi}^c_\nu \kern+0.2pt \bigr)$')

plt.xscale('log')
plt.yscale('log')

plt.xlim(2e2, 1e6)
plt.ylim(1e-4, 2e0)

plt.legend(loc=2)

plt.savefig('build/present_magnetar_charm_without.pdf')
plt.savefig('build/present_magnetar_charm_without.png')
plt.close()

plt.fill_between(t, 3 * c / N, c / (3 * N), color='none', facecolor=TUDO, alpha=0.3, zorder=0)

plt.plot(t, pi / N, label=r'Pion Decay', zorder=1, c=ACC)
plt.plot(t, K / N, label=r'Kaon Decay', zorder=1, c='k')
plt.plot(t, c / N, label=r'Charm Decay', zorder=0, c=TUDO)

plt.xlabel(r'$t$ $\mathrel{/}$ $\symup{s}$')
plt.ylabel(r'$\dot{\phi}_{\kern-0.3pt \nu}$ $\mathrel{/} \kern-0.7pt$ $\symup{max} \kern+0.5pt \bigl( \dot{\phi}^c_\nu \kern+0.2pt \bigr)$')

plt.xscale('log')
plt.yscale('log')

plt.xlim(4e1, 4e6)
plt.ylim(1e-6, 1e1)

plt.legend(loc=2)

plt.savefig('build/present_magnetar_flux_without.pdf')
plt.savefig('build/present_magnetar_flux_without.png')
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

plt.scatter([], [], label=r'Pion Decay', c=ACC)
plt.scatter([], [], label=r'Kaon Decay', c='k')
plt.scatter([], [], label=r'Charm Decay', c=TUDO)

plt.fill_between(E, 3 * E**2 * c3 / N, E**2 * c3 / (3 * N), color='none', facecolor=TUDO, alpha=0.3)

plt.plot(E, E**2 * c1 / N, ':', c=TUDO)
plt.plot(E, E**2 * c2 / N, '--', c=TUDO)
plt.plot(E, E**2 * c3 / N, '-', c=TUDO)
plt.plot(E, E**2 * pi1 / N, ':', c=ACC)
plt.plot(E, E**2 * pi2 / N, '--', c=ACC)
plt.plot(E, E**2 * pi3 / N, '-', c=ACC)
plt.plot(E, E**2 * K1 / N, 'k:', label=r'$10^3 \kern+1.5pt \symup{s}$ $\kern+0.8pt -$ $10^4 \kern+1.5pt \symup{s}$')
plt.plot(E, E**2 * K2 / N, 'k--', label=r'$10^4 \kern+1.5pt \symup{s}$ $\kern+0.8pt -$ $10^5 \kern+1.5pt \symup{s}$')
plt.plot(E, E**2 * K3 / N, 'k-', label=r'$10^3 \kern+1.5pt \symup{s}$ $\kern+0.8pt -$ $10^7 \kern+1.5pt \symup{s}$') 

plt.xlabel(r'$E_\nu$ $\mathrel{/}$ $\symup{GeV}$')
plt.ylabel(r'$E_\nu^2\phi_\nu$ $\mathrel{/}$ $\symup{max} \kern+0.5pt \bigl( \kern-0.3pt E_\nu^2\phi^c_\nu \kern+0.2pt \bigr)$')

plt.xscale('log')
plt.yscale('log')

plt.xlim(3e5, 1e11)
plt.ylim(1.6e-4, 3e3)

plt.legend(loc=1)

plt.savefig('build/present_magnetar_fluence_without.pdf')
plt.savefig('build/present_magnetar_fluence_without.png')
plt.close()



plt.figure(figsize=(5.5, 3.5))

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

plt.plot(t, c / N, '-', label=r'Total Charm Decay', zorder=2, c=TUDO)
plt.plot(t, D0 / N, '--', label=r'$D^0$ Decay', zorder=2, c=TUDO)
plt.plot(t, Dplus / N, '-', label=r'$D^+$ Decay', zorder=1, c=ACC)
plt.plot(t, DplusS / N, '--', label=r'$D^+_s$ Decay', zorder=1, c=ACC)
plt.plot(t, LAMplusC / N, 'k-', label=r'$\Lambda^{\kern-0.5pt +}_{\kern+0.5pt c}$ Decay', zorder=0)

plt.xlabel(r'$t$ $\mathrel{/}$ $\symup{s}$')
plt.ylabel(r'$\dot{\phi}_{\kern-0.3pt \nu}$ $\mathrel{/} \kern-0.7pt$ $\symup{max} \kern+0.5pt \bigl( \dot{\phi}^c_\nu \kern+0.2pt \bigr)$')

plt.xscale('log')
plt.yscale('log')

plt.xlim(2e2, 1e6)
plt.ylim(1e-4, 2e0)

plt.legend(loc=3)

plt.savefig('build/present_magnetar_charm_with.pdf')
plt.savefig('build/present_magnetar_charm_with.png')
plt.close()

plt.fill_between(t, 3 * c / N, c / (3 * N), color='none', facecolor=TUDO, alpha=0.3, zorder=0)

plt.plot(t, pi / N, label=r'Pion Decay', zorder=1, c=ACC)
plt.plot(t, K / N, label=r'Kaon Decay', zorder=1, c='k')
plt.plot(t, c / N, label=r'Charm Decay', zorder=0, c=TUDO)

plt.xlabel(r'$t$ $\mathrel{/}$ $\symup{s}$')
plt.ylabel(r'$\dot{\phi}_{\kern-0.3pt \nu}$ $\mathrel{/} \kern-0.7pt$ $\symup{max} \kern+0.5pt \bigl( \dot{\phi}^c_\nu \kern+0.2pt \bigr)$')

plt.xscale('log')
plt.yscale('log')

plt.xlim(4e1, 4e6)
plt.ylim(1e-6, 1e1)

plt.legend(loc=3)

plt.savefig('build/present_magnetar_flux_with.pdf')
plt.savefig('build/present_magnetar_flux_with.png')
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

plt.scatter([], [], label=r'Pion Decay', c=ACC)
plt.scatter([], [], label=r'Kaon Decay', c='k')
plt.scatter([], [], label=r'Charm Decay', c=TUDO)

plt.fill_between(E, 3 * E**2 * c3 / N, E**2 * c3 / (3 * N), color='none', facecolor=TUDO, alpha=0.3)

plt.plot(E, E**2 * c1 / N, ':', c=TUDO)
plt.plot(E, E**2 * c2 / N, '--', c=TUDO)
plt.plot(E, E**2 * c3 / N, '-', c=TUDO)
plt.plot(E, E**2 * pi1 / N, ':', c=ACC)
plt.plot(E, E**2 * pi2 / N, '--', c=ACC)
plt.plot(E, E**2 * pi3 / N, '-', c=ACC)
plt.plot(E, E**2 * K1 / N, 'k:', label=r'$10^3 \kern+1.5pt \symup{s}$ $\kern+0.8pt -$ $10^4 \kern+1.5pt \symup{s}$')
plt.plot(E, E**2 * K2 / N, 'k--', label=r'$10^4 \kern+1.5pt \symup{s}$ $\kern+0.8pt -$ $10^5 \kern+1.5pt \symup{s}$')
plt.plot(E, E**2 * K3 / N, 'k-', label=r'$10^3 \kern+1.5pt \symup{s}$ $\kern+0.8pt -$ $10^7 \kern+1.5pt \symup{s}$') 

plt.xlabel(r'$E_\nu$ $\mathrel{/}$ $\symup{GeV}$')
plt.ylabel(r'$E_\nu^2\phi_\nu$ $\mathrel{/}$ $\symup{max} \kern+0.5pt \bigl( \kern-0.3pt E_\nu^2\phi^c_\nu \kern+0.2pt \bigr)$')

plt.xscale('log')
plt.yscale('log')

plt.xlim(3e5, 1e11)
plt.ylim(4.3e-4, 6e1)

plt.legend(loc=1)

plt.savefig('build/present_magnetar_fluence_with.pdf')
plt.savefig('build/present_magnetar_fluence_with.png')
plt.close()



plt.figure(figsize=(5.5, 3.5))

E, pi = np.genfromtxt('code/tabulate/nucleus/neutrinos/pi.txt', unpack=True)
E, K = np.genfromtxt('code/tabulate/nucleus/neutrinos/K.txt', unpack=True)
E, D0 = np.genfromtxt('code/tabulate/nucleus/neutrinos/D0.txt', unpack=True)
E, Dplus = np.genfromtxt('code/tabulate/nucleus/neutrinos/Dplus.txt', unpack=True)
E, DplusS = np.genfromtxt('code/tabulate/nucleus/neutrinos/DplusS.txt', unpack=True)
E, LAMplusC = np.genfromtxt('code/tabulate/nucleus/neutrinos/LAMplusC.txt', unpack=True)

c = D0 + Dplus + DplusS + LAMplusC

N = (E**2 * c).max()

plt.plot(E, E**2 * c / N, '-', label=r'Total Charm Decay', zorder=2, c=TUDO)
plt.plot(E, E**2 * D0 / N, '--', label=r'$D^0$ Decay', zorder=2, c=TUDO)
plt.plot(E, E**2 * Dplus / N, '-', label=r'$D^+$ Decay', zorder=1, c=ACC)
plt.plot(E, E**2 * DplusS / N, '--', label=r'$D^+_s$ Decay', zorder=1, c=ACC)
plt.plot(E, E**2 * LAMplusC / N, 'k-', label=r'$\Lambda^{\kern-0.5pt +}_{\kern+0.5pt c}$ Decay', zorder=0)

plt.xlabel(r'$E_\nu$ $\mathrel{/}$ $\symup{GeV}$')
plt.ylabel(r'$E_\nu^2\phi_\nu$ $\mathrel{/}$ $\symup{max} \kern+0.5pt \bigl( \kern-0.3pt E_\nu^2\phi^c_\nu \kern+0.2pt \bigr)$')

plt.xscale('log')
plt.yscale('log')

plt.xlim(3e5, 1e11)
plt.ylim(5e-4, 2e0)

plt.legend(loc=3)

plt.savefig('build/present_nucleus_charm.pdf')
plt.savefig('build/present_nucleus_charm.png')
plt.close()

plt.figure(figsize=(5.5, 3.5))

plt.fill_between(E, 3 * E**2 * c / N, E**2 * c / (3 * N), color='none', facecolor=TUDO, alpha=0.3)

plt.plot(E, E**2 * pi / N, label=r'Pion Decay', zorder=1, c=ACC)
plt.plot(E, E**2 * K / N, label=r'Kaon Decay', zorder=1, c='k')
plt.plot(E, E**2 * c / N, label=r'Charm Decay', zorder=0, c=TUDO)

plt.xlabel(r'$E_\nu$ $\mathrel{/}$ $\symup{GeV}$')
plt.ylabel(r'$E_\nu^2\phi_\nu$ $\mathrel{/}$ $\symup{max} \kern+0.5pt \bigl( \kern-0.3pt E_\nu^2\phi^c_\nu \kern+0.2pt \bigr)$')

plt.xscale('log')
plt.yscale('log')

plt.xlim(3e5, 1e11)
plt.ylim(5e-4, 2e3)

plt.legend(loc=1)

plt.savefig('build/present_nucleus_fluence.pdf')
plt.savefig('build/present_nucleus_fluence.png')
plt.close()
