#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import sys, argparse
from ase.io import read, write
from ase.collections import g2
from ase.build import molecule

def del_mol(cml):
    arg = parse_cml_args(cml)

    slab = read(arg.poscar)
    remaing_atoms = [ii for ii in range(slab.get_number_of_atoms())
                     if (ii+1) not in arg.atoms]
    new = slab[remaing_atoms]

    # add vacuum 
    if arg.vacuum:
        if arg.xvacuum == 'a':
            axis = (0, 1, 2)
        else:
            axis = 'xyz'.index(arg.xvacuum)
        new.center(vacuum=arg.vacuum / 2., axis=axis)
    else:
        # By default, center the slab in the z direction and keep the original vacuum
        # length.
        new.center(axis=2)

    write(arg.out, new, vasp5=True, direct=True,
          label=open(arg.poscar).readline().strip())

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
    arg.add_argument('-a', '--atoms', dest='atoms', action='store', type=int,
                     default=[], nargs='*',
                     help='The indices of the atoms that you want to delete.')
    arg.add_argument('-v', '--vacuum', dest='vacuum', action='store', type=float,
                     default=None, 
                     help='Set new vacuum length.')
    arg.add_argument('--xvacuum', dest='xvacuum', action='store', type=str,
                     default='z', choices=['x', 'y', 'z', 'a'],
                     help='Vacuum direction.')
    arg.add_argument('--no-sort-pos', dest='sort_pos', action='store_false',
                     help='Sort the coordinates.')

    return arg.parse_args(cml)

if __name__ == '__main__':
    del_mol(sys.argv[1:])
