#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import sys, argparse
from ase.io import read, write

def parse_cml_args(cml):
    '''
    CML parser.
    '''
    arg = argparse.ArgumentParser(add_help=True)

    arg.add_argument('-i', dest='poscar', action='store', type=str,
                     default='POSCAR',
                     help='The POSCAR based on which to make supercell.')
    arg.add_argument('-o', dest='out', action='store', type=str,
                     default=None,
                     help='Default output filename.')
    arg.add_argument('-s', '--size', dest='size', action='store', type=int,
                     default=[1, 1, 1], nargs=3,
                     help='The supercell size.')
    arg.add_argument('-n', dest='new_sym_order', action='store',
                     type=str, default=None, nargs='*',
                     help='New order of the chemical symbols.')
    arg.add_argument('--no-sort-pos', dest='sort_pos', action='store_false',
                     help='Sort the coordinates.')
    arg.add_argument('-v', '--vacuum', dest='vacuum', action='store', type=float,
                     default=None, 
                     help='Set new vacuum length.')
    arg.add_argument('--ivacuum', dest='ivacuum', action='store', type=str,
                     default='z', choices=['x', 'y', 'z'],
                     help='Vacuum direction.')

    return arg.parse_args(cml)

def mk_supercell(cml):
    arg = parse_cml_args(cml)

    pc = read(arg.poscar)
    sc = pc * arg.size
    org_chem_symbols = np.array(sc.get_chemical_symbols())

    # add vacuum 
    if arg.vacuum:
        sc.center(vacuum=arg.vacuum / 2., axis='xyz'.index(arg.ivacuum))

    org_atom_index = np.arange(len(sc), dtype=int)

    # Sort first by z-coordinates then by y and x.
    if arg.sort_pos:
        rpos = np.round(sc.positions, 4)
        new_atom_index = np.lexsort(
                (rpos[:, 0], rpos[:, 1], rpos[:, 2])
                )
        sc = sc[new_atom_index]

    # New order of chemical symbols
    if arg.new_sym_order:
        assert set(arg.new_sym_order) == set(org_chem_symbols)
        chem_sym_order = arg.new_sym_order
    else:
        # Just stick to the original order.
        chem_sym_order = []
        for ss in org_chem_symbols:
            if not ss in chem_sym_order:
                chem_sym_order.append(ss)

    # Re-arrange the atoms according to the new order of chemical symbols
    new_atom_index = [ii for ss in chem_sym_order
                         for ii in org_atom_index[sc.symbols == ss]]
    sc = sc[new_atom_index]

    if arg.out:
        write(arg.out, sc, vasp5=True, direct=True,
              label=open(arg.poscar).readline().strip())
    else:
        write('out_' + 'x'.join(["%d" % x for x in arg.size]) + '.vasp',
              sc, vasp5=True, direct=True,
              label=open(arg.poscar).readline().strip())

if __name__ == '__main__':
    mk_supercell(sys.argv[1:])
