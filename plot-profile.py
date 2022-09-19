#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

plt.style.use('bmh')

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

# densidade eletronica de perfil exponencial
Ne_expprf = lambda r: 6.5956e+04 * np.exp(-10.54 * r)

# densidade eletronica interpolada do Bahcall
r, Ne_interp = np.loadtxt('elecdens.txt', unpack=True)
eV_cm = 50677.30716548338   # eV . cm
R_sun = 6.957e+10   # centimetros
G_FxN_A = 7.0240967108658126    # Fermi x Avogadro em eV^{-2}
Ne_interp *= np.sqrt(2) * G_FxN_A * R_sun * (1/eV_cm)**2

plt.plot(r, Ne_expprf(r), label=r'exponencial')
plt.plot(r, Ne_interp, label=r'interpolado')
plt.xlabel(r'$r = R / R_{\odot}$', fontsize=20)
plt.ylabel(r'$\sqrt{2} G_F N_e$', fontsize=20)
plt.legend(fontsize=16)
plt.title(r'Comparação entre os perfis de densidade eletrônica')
plt.savefig("./perfis-elecdens.png", dpi=300, format='png', bbox_inches="tight")
#plt.show()


#print("integral =", np.trapz(nu, r))
# O calculo acima resulta em 0.9999674861119427
