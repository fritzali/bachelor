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

plt.xlabel(r'$x_E = \kern-0.5pt E_h \kern+0.5pt / E_p$')
plt.ylabel(r'$x_E \kern+1.5pt d\sigma \kern-0.3pt / \kern-0.8pt dx_E \kern+1.2pt$ $\mathrel{/}$ mb')

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
