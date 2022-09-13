#!/usr/bin/env python3

import numpy as np

from matplotlib import pyplot as plt
from matplotlib import rc
plt.style.use('bmh')
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

fl = "{:15.5e}"
titulo = r'Probabilidade de transição $\nu_e \to \nu_\alpha$ no interior do Sol para 3 gerações'
energy = 6.44   # in MeV

i, ti, electron = np.loadtxt('output-surv.txt', unpack=True)

plt.plot(ti, electron  , label=r'$P(\nu_e \to \nu_e)$')
#plt.plot(ti, mu        , label=r'$P(\nu_e \to \nu_\mu)$')
#plt.plot(ti, tau       , label=r'$P(\nu_e \to \nu_\tau)$')
#plt.plot(ti, 1-electron, label=r'$P(\nu_e \to \nu_\mu)$')
plt.xlabel(r'$r = R/R_\odot$', fontsize=20)
plt.ylabel(r'Probabilidade', fontsize=20)
plt.ylim(0, 1.05)
plt.legend(fontsize=14)
plt.title(titulo)
plt.savefig('fig/3nu-energy_%.2f.png' % (energy), dpi=300, format='png', bbox_inches="tight")
