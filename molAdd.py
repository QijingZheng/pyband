#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import sys, argparse
from ase.io import read, write
from ase.collections import g2
from ase.build import molecule

def add_mol(cml):
    arg = parse_cml_args(cml)

    slab = read(arg.slab)

    assert arg.molecule is not None, "Please specify the name of the adsorbates."
    mol = molecule(arg.molecule, pbc=slab.pbc, cell=slab.cell)
    mol.center()

    print (' '.join(mol.get_chemical_symbols()))

    if arg.rot:
        for za in arg.rot:
            if za[0].lower() in 'xyz':
                axis = za[0].lower()
                angle = za[1:]
            else:
                axis = 'z'
                angle = za

            try:
                angle = float(angle)
            except:
                print("Please enter valid rotation parameter, e.g. z90!")
                raise
            mol.rotate(angle, axis)

    if arg.rotx:
        mol.rotate(arg.rotx, 'x')
    if arg.roty:
        mol.rotate(arg.roty, 'y')
    if arg.rotz:
        mol.rotate(arg.rotz, 'z')

    mol.positions += slab.positions[arg.atom_index] - mol.positions[arg.mol_index]
    mol.positions += [arg.offset[0], arg.offset[1], arg.height]
    new = slab + mol

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

    org_atom_index = np.arange(len(new), dtype=int)

    # Sort first by z-coordinates then by y and x.
    if arg.sort_pos:
        rpos = np.round(new.positions, 4)
        new_atom_index = np.lexsort(
                (rpos[:, 0], rpos[:, 1], rpos[:, 2])
                )
        new = new[new_atom_index]

    # Just stick to the original order.
    chem_sym_order = []
    for ss in slab.get_chemical_symbols() + mol.get_chemical_symbols():
        if not ss in chem_sym_order:
            chem_sym_order.append(ss)

    # Re-arrange the atoms according to the new order of chemical symbols
    new_atom_index = [ii for ss in chem_sym_order
                         for ii in org_atom_index[new.symbols == ss]]
    new = new[new_atom_index]

    write(arg.out, new, vasp5=True, direct=True,
          label=open(arg.slab).readline().strip())

def parse_cml_args(cml):
    '''
    CML parser.
    '''
    arg = argparse.ArgumentParser(add_help=True)

    arg.add_argument('-i', dest='slab', action='store', type=str,
                     default='POSCAR',
                     help='The slab onto which the adsorbate should be added.')
    arg.add_argument('-o', dest='out', action='store', type=str,
                     default='out.vasp',
                     help='Default output filename.')
    arg.add_argument('-m', '--molecule', dest='molecule', action='store', type=str,
                     choices=g2.names, default=None,
                     help='The chemical formula of the adsorbed molecule, e.g.  H2O.')
    arg.add_argument('--height', dest='height', action='store',
                     type=float, default=1.0, 
                     help='Height of the the adsorbate above the surface.')
    arg.add_argument('--mol_index', dest='mol_index', action='store',
                     type=int, default=0, 
                     help='Index of the atom in the molecule to be positioned above the adsorbed site.')
    arg.add_argument('-a', '--atom_index', dest='atom_index', action='store',
                     type=int, default=0, 
                     help='Index of the atom onto which the adsorbate should be added.')
    arg.add_argument('--offset', dest='offset', action='store',
                     type=float, default=(0.0, 0.0), 
                     help='Offset the adsorbate.')
    arg.add_argument('--rot', dest='rot', action='store',
                     type=str, default=[], nargs='+',
                     help='Rotate around a specified axis by an angle, e.g.  --rot z90 --- rotate 90 degrees around z-axis.')
    arg.add_argument('--rotx', dest='rotx', action='store',
                     type=float, default=None, 
                     help='Rotation around x-axis.')
    arg.add_argument('--roty', dest='roty', action='store',
                     type=float, default=None, 
                     help='Rotation around y-axis.')
    arg.add_argument('--rotz', dest='rotz', action='store',
                     type=float, default=None, 
                     help='Rotation around z-axis.')
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
    add_mol(sys.argv[1:])
