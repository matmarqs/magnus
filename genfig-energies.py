#!/usr/bin/env python3

import numpy as np

from matplotlib import pyplot as plt
from matplotlib import rc
plt.style.use('bmh')
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

fl = "{:15.5e}"
titulo = r'Probabilidade de sobrevivência de neutrinos solares'
#energy = 6.44   # in MeV

energy, prob = np.loadtxt('output/energiesXsurv_prob.txt', unpack=True)

energy = np.log10(energy)

plt.plot(energy, prob, label=r'$P_{ee}$')
#plt.plot(ti, electron  , label=r'$P(\nu_e \to \nu_e)$')
#plt.plot(ti, mu        , label=r'$P(\nu_e \to \nu_\mu)$')
#plt.plot(ti, tau       , label=r'$P(\nu_e \to \nu_\tau)$')
#plt.plot(ti, 1-electron, label=r'$P(\nu_e \to \nu_\mu)$')
plt.xlabel(r'log$(E_\nu)$ em MeV', fontsize=20)
plt.ylabel(r'probabilidade', fontsize=20)
plt.ylim(0, 1.05)
plt.legend(fontsize=14)
plt.title(titulo)
plt.savefig('fig/energies-surv.png', dpi=300, format='png', bbox_inches="tight")
