#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import numpy as np
from math import pi, sqrt
import argparse

############################################################
# The following data are extracted from phonopy.units

PlanckConstant = 4.13566733e-15         # [eV s]
Hbar = PlanckConstant/(2*pi)            # [eV s]
Avogadro = 6.02214179e23
SpeedOfLight = 299792458                # [m/s]
EV = 1.60217733e-19                     # [J]
kB = 1.3806504e-23 / EV                 # [eV/K] 8.6173383e-05

THzToEv = PlanckConstant * 1e12         # [eV]
THzToCm = 1.0e12 / (SpeedOfLight * 100) # [cm^-1] 33.356410
CmToEv = THzToEv / THzToCm              # [eV] 1.2398419e-4
EvTokJmol = EV / 1000 * Avogadro        # [kJ/mol] 96.4853910
############################################################

# The conversion matrix, the row and colomn corresponds to
# ['eV', 'kJ/mol', 'cm-1', 'nm', 'THz', 'fs', 'K']

conv_func = [
    lambda x: x,                    # eV to eV
    lambda x: x * EvTokJmol,        # eV to kJ/mol
    lambda x: x / CmToEv,           # eV to cm-1
    lambda x: 1E7 / (x / CmToEv),   # eV to nm
    lambda x: x / THzToEv,          # eV to THz
    lambda x: 1000 / (x / THzToEv), # eV to fs
    lambda x: x / kB,               # eV to Kelvin
]

conv_func_inv = [
    lambda x: x,                  # eV to eV
    lambda x: x / EvTokJmol,      # kJ/mol to eV
    lambda x: x * CmToEv,         # cm-1 to eV
    lambda x: 1E7 * CmToEv / x,   # nm to eV
    lambda x: x * THzToEv,        # THz to eV
    lambda x: 1000 * THzToEv / x, # eV to cm-1
    lambda x: x * kB,             # Kelvin to eV
]

if __name__ == "__main__":
    available_units = ['eV', 'kJ/mol', 'cm-1', 'nm', 'THz', 'fs', 'K']

    arg = argparse.ArgumentParser(add_help=True)

    arg.add_argument('-u', dest='unit', action='store', type=str,
                     default='eV',
                     choices=available_units,
                     help='The unit of the input numbers.')
    arg.add_argument('values', metavar='V', type=str,
                     nargs='+',
                     help='The input values')

    p = arg.parse_args()

    print(''.join(["{:^12s}".format(u) for u in available_units]))
    print('=' * 12 * len(available_units))
    for val in p.values:
        unit_idx = [u.lower() for u in available_units].index(p.unit.lower())
        # first convert to eV
        fval = conv_func_inv[unit_idx](float(val))
        print(''.join([
            "{:^12.6G}".format(C(fval)) for C in conv_func])
        )
