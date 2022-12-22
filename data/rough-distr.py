#!/usr/bin/env python3

import sys
import numpy as np

r, p = np.loadtxt(sys.argv[1], unpack=True)

n_r = len(r)
num_zones   = 45
num_zones_2 = int(num_zones/2)
z = 0
while (num_zones * z < n_r):
    print(" {:.5f}   {:.7e}".format(r[num_zones_2*(2*z+1)], sum(p[num_zones*z:num_zones*(z+1)])))
    z += 1

#print(sum(p))
#for i in range(num):
#    print(" {:.7f}   {:.10e}".format(R[i], 10.0 ** lnNe[i]))
