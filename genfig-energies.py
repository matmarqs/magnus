#!/usr/bin/env python3

import numpy as np

#magnus = 2
magnus = 2
profile = '-interp'
#profile = ''
#tipo = 'log'
tipo = ''


from matplotlib import pyplot as plt
from matplotlib import rc
plt.style.use('bmh')
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

#fl = "{:15.5e}"
titulo = f"Probabilidade de sobrevivência de neutrinos solares (M{magnus}{profile})"
#energy = 6.44   # in MeV

energy, prob = np.loadtxt(f"output/energiesXsurv_prob-{tipo}Magnus{magnus}{profile}.txt", unpack=True)

plt.plot(energy, prob, label=r'$P_{ee}$')
#plt.plot(ti, electron  , label=r'$P(\nu_e \to \nu_e)$')
#plt.plot(ti, mu        , label=r'$P(\nu_e \to \nu_\mu)$')
#plt.plot(ti, tau       , label=r'$P(\nu_e \to \nu_\tau)$')
#plt.plot(ti, 1-electron, label=r'$P(\nu_e \to \nu_\mu)$')
plt.xlabel(r'$E_\nu$ em MeV', fontsize=20)
if tipo == 'log':
    plt.xscale('log')
plt.ylabel(r'probabilidade', fontsize=20)
#plt.ylim(0, 1.05)
plt.legend(fontsize=14)
plt.title(titulo)
plt.savefig(f"fig/energies-surv-{tipo}Magnus{magnus}{profile}.png", dpi=300, format='png', bbox_inches="tight")
