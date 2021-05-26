#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import ase
import sys
import numpy as np
from ase.io import read, write
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary, ZeroRotation


def init_vel(cml):
    '''
    '''
    arg = parse_cml_args(cml)

    # read in the initial structure
    init_pos = read(arg.poscar)

    # set the momenta corresponding to T
    MaxwellBoltzmannDistribution(
        init_pos, temperature_K=arg.temperature, force_temp=True
    )
    # set the center-of-mass to 0
    Stationary(init_pos)
    # set the total angular momentum to 0
    ZeroRotation(init_pos)

    ################################################################################
    # No need to do this, use "force_temp" parameter in MaxwellBoltzmannDistribution
    ################################################################################

    # scale the temperature to T
    # vel = init_pos.get_velocities()
    # Tn = init_pos.get_temperature()
    # vel *= np.sqrt(arg.temperature / Tn)
    # init_pos.set_velocities(vel)

    # write the structure
    write(arg.out, init_pos, vasp5=True, direct=True)

    # units in VASP and ASE are different
    vel = init_pos.get_velocities() * ase.units.fs
    # np.savetxt('init_vel.dat', vel, fmt='%20.16f')
    # append the velocities to the POSCAR
    with open(arg.out, 'a+') as pos:
        pos.write('\n')
        pos.write(
            '\n'.join([
                ''.join(["%20.16f" % x for x in row])
                for row in vel
            ])
        )


def parse_cml_args(cml):
    '''
    CML parser.
    '''
    arg = argparse.ArgumentParser(add_help=True)

    arg.add_argument('-i', dest='poscar', action='store', type=str,
                     default='POSCAR',
                     help='The slab POSCAR.')
    arg.add_argument('-o', dest='out', action='store', type=str,
                     default='out.vasp',
                     help='Default output filename.')
    arg.add_argument('-t', '--temperature', dest='temperature',
                     action='store', type=float,
                     default=300,
                     help='The temperature.')

    return arg.parse_args(cml)


if __name__ == '__main__':
    init_vel(sys.argv[1:])
