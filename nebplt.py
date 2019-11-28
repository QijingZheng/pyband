#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from subprocess import call

plt.style.use('ggplot')

############################################################
# Generating and Loading data
# call vtst scripts to generate the MEP.
call('nebbarrier.pl')
call('nebspline.pl')

# scatters of the barrier
mep_d = np.loadtxt('neb.dat')
# splines of the barrier
mep_s = np.loadtxt('spline.dat')
############################################################
# Plotting

nrows = 1; ncols = 1
fig = plt.figure(
    figsize = (4.8, 2.4)
)
axes = [
    plt.subplot(nrows, ncols, ii+1)
    for ii in range(ncols * nrows)
]

axes[0].plot(mep_s[:,1], mep_s[:,2], ls='-', lw=0.5, color='r', alpha=0.8)
axes[0].plot(mep_d[:,1], mep_d[:,2], ls='none',
             marker='*', ms=6, mew=0, mfc='b')

axes[0].set_xlabel('Reaction Coordinate', labelpad=5)
axes[0].set_ylabel('Energy [eV]', labelpad=5)

plt.tight_layout()
# plt.show()
plt.savefig('mep.png', dpi=300)

try:
    from subprocess import call
    call('feh -xdF mep.png'.split())
except:
    pass

