import matplotlib.pyplot as plt
import numpy as np

from code.functional import *

n = 1e15
d = 1e15

Ep = np.logspace(5, 13, 1001)
Sp = 1 / Ep**2

od = proton_proton_optical_depth(Ep, n, d)

Sp = od * Sp

Eh = np.logspace(5, 13, 1000)

x = Eh[:, None] / Ep[None, :]

Fph = meson_production(x, Ep, 'pi')

dp = np.diag(np.insert((Ep[1:] - Ep[:-1]), 0, 0.0))

Sh = Fph @ dp @ Sp

cf = hadron_proton_cooling_factor(Eh, n, 'pi')

Sh = cf * Sh

Enu = np.logspace(5, 12, 999)

Fhnu = meson_decay_neutrinos(Enu[:, None], Eh[None, :], 'pi')

print(Fhnu)

dh = np.diag(np.insert((Eh[1:] - Eh[:-1]), 0, 0.0))

Snu = Fhnu @ dh @ Sh

plt.plot(Ep, Sp / Sp[100])
plt.plot(Eh, Sh / Sh[100])
plt.plot(Enu, Snu / Snu[100])

plt.xscale('log')
plt.yscale('log')

plt.show()
