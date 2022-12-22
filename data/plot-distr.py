#!/usr/bin/env python3

#import sys
import numpy as np

from scipy import integrate

import matplotlib.pyplot as plt
from matplotlib import rc
plt.style.use('bmh')
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
from labellines import labelLines

r, pp, b, n, o, f, be, pep, hep = np.loadtxt('fractions.txt', unpack=True)

pp = pp  / integrate.simpson(pp,  r)
b = b  / integrate.simpson(b,  r)
n = n / integrate.simpson(n, r)
o = o / integrate.simpson(o, r)
f = f / integrate.simpson(f, r)
be = be / integrate.simpson(be, r)
pep = pep / integrate.simpson(pep, r)
hep = hep / integrate.simpson(hep, r)

plt.plot(r, pp, label=r'pp')
plt.plot(r, b, label=r'$^8$B')
plt.plot(r, n, label=r'$^{13}$N')
plt.plot(r, o, label=r'$^{15}$O')
plt.plot(r, f, label=r'$^{17}$F')
plt.plot(r, be, label=r'$^7$Be')
plt.plot(r, pep, label=r'pep')
plt.plot(r, hep, label=r'hep')
xvals = [0.058, 0.032, 0.029, 0.03, 0.031, 0.082, 0.058, 0.058]
labelLines(plt.gca().get_lines(), align=False, xvals=xvals, fontsize=13)

plt.xlabel(r'$r = R / R_{\odot}$', fontsize=20)
plt.ylabel(r'densidade de probabilidade', fontsize=16)
plt.legend(fontsize=16)
plt.title(r'Densidade de produção radial para neutrinos de diferentes reações')
plt.savefig("./raio.png", dpi=300, format='png', bbox_inches="tight")
#plt.show()


#print("integrate.simpson =", np.trapz(nu, r))
# O calculo acima resulta em 0.9999674861119427
