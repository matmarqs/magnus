#!/usr/bin/env python3

import sys
import numpy as np

my_dpi = 96
import matplotlib.pyplot as plt
from matplotlib import rc
plt.style.use('bmh')
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

#from labellines import labelLines

energy, pp, b, n, o, f, be, pep, hep = np.loadtxt(sys.argv[1], unpack=True)

#fl = "{:15.5e}"
titulo = f"Média da probabilidade de sobrevivência para diferentes reações"
#energy = 6.44   # in MeV

#plt.plot(energy, prob, label=r'$P_{ee}$')
#plt.plot(ti, electron  , label=r'$P(\nu_e \to \nu_e)$')
#plt.plot(ti, mu        , label=r'$P(\nu_e \to \nu_\mu)$')
#plt.plot(ti, tau       , label=r'$P(\nu_e \to \nu_\tau)$')
#plt.plot(ti, 1-electron, label=r'$P(\nu_e \to \nu_\mu)$')

plt.plot(energy, pp, label=r'pp')
plt.plot(energy, b, label=r'$^8$B')
plt.plot(energy, n, label=r'$^{13}$N')
plt.plot(energy, o, label=r'$^{15}$O')
plt.plot(energy, f, label=r'$^{17}$F')
plt.plot(energy, be, label=r'$^7$Be')
plt.plot(energy, pep, label=r'pep')
plt.plot(energy, hep, label=r'hep')
#xvals = [0.058, 0.032, 0.029, 0.03, 0.031, 0.082, 0.058, 0.058]
#labelLines(plt.gca().get_lines(), align=False, xvals=xvals, fontsize=13)
plt.xscale('log')
plt.xlabel(r'$E$ (MeV)', fontsize=20)
plt.ylabel(r'$\langle P(\nu_e \to \nu_e) \rangle$', fontsize=20)
#plt.ylim(0, 1.05)
plt.legend(fontsize=14)
plt.title(titulo)
plt.savefig(f"surv_prob.png", dpi=300, format='png', bbox_inches="tight")
#plt.show()
