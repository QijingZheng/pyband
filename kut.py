#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, argparse
import numpy as np
from ase.atoms import Atoms
from ase.build import surface
from ase.io import read, write
from ase.constraints import FixAtoms

def parse_cml_args(cml):
    '''
    CML parser.
    '''
    arg = argparse.ArgumentParser(add_help=True)

    arg.add_argument('-i', dest='poscar', action='store', type=str,
                     default='POSCAR',
                     help='The primitive cell base on which to cut the slab.')
    arg.add_argument('-o', dest='out', action='store', type=str,
                     default=None,
                     help='Default output filename.')
    arg.add_argument('--hkl', dest='hkl', action='store', type=int,
                     default=[1, 1, 1], nargs=3,
                     help='Surface normal in Miller indices (h,k,l).')
    arg.add_argument('-n', dest='nlayer', action='store', type=int,
                     default=3,
                     help='Number of layers of the slab.')
    arg.add_argument('-l', dest='layer_thickness', action='store', type=float,
                     default=0.1,
                     help='Default thickness of each atomic layer.')
    arg.add_argument('-s', dest='new_sym_order', action='store',
                     type=str, default=None, nargs='*',
                     help='New order of the chemical symbols.')
    arg.add_argument('-d', dest='delete_atomic_layer', action='store',
                     type=int, default=None, nargs='*',
                     help='Delete the unwanted atomic layer.')
    arg.add_argument('-f', dest='fix_atomic_layers', action='store',
                     type=int, default=None, nargs='*',
                     help='Fix the atoms in the specified atomic layers.')
    arg.add_argument('-v', '--vacuum', dest='vacuum', action='store', type=float,
                     default=15., 
                     help='Set new vacuum length.')
    arg.add_argument('--ivacuum', dest='ivacuum', action='store', type=str,
                     default='z', choices=['x', 'y', 'z'],
                     help='Vacuum direction.')

    return arg.parse_args(cml)

def find_natoms_layers(pos_z, layer_thickness):
    '''
    '''
    natoms = pos_z.size
    natoms_per_layers = []
    indices_for_layers = []
    remaining_atoms = np.ones(natoms, dtype=bool)

    while np.sum(natoms_per_layers) != natoms:
        layer_center = pos_z[remaining_atoms].min()
        layer_indices = np.arange(natoms, dtype=int)[
                (layer_center - layer_thickness/2. < pos_z)
                &
                (layer_center + layer_thickness/2. > pos_z)
            ]

        remaining_atoms[layer_indices] = False
        indices_for_layers.append(layer_indices)
        natoms_per_layers.append(len(layer_indices))
        
        if np.sum(remaining_atoms) == 0:
            break

    n_atomic_layers = len(natoms_per_layers)

    return n_atomic_layers, natoms_per_layers, indices_for_layers


def cut_slab(cml):
    '''
    '''
    args = parse_cml_args(cml)

    primitive_cell = read(args.poscar)
    slab = surface(primitive_cell, args.hkl, layers=args.nlayer)

    pos_z = slab.positions.copy()[:, 'xyz'.index(args.ivacuum)]
    n_atomic_layers, natoms_per_layers, indices_for_layers = \
    find_natoms_layers(pos_z, args.layer_thickness)
    print("{} atomic layers found!".format(n_atomic_layers))

    if args.delete_atomic_layer:
        kept_atoms = [
                ii for jj in range(n_atomic_layers)
                   for ii in indices_for_layers[jj]
                   if jj not in args.delete_atomic_layer
                ]
        slab = slab[kept_atoms]

        pos_z = slab.positions.copy()[:, 'xyz'.index(args.ivacuum)]
        n_atomic_layers, natoms_per_layers, indices_for_layers = \
        find_natoms_layers(pos_z, args.layer_thickness)

    if args.fix_atomic_layers:
        C = FixAtoms(indices=[ii for jj in range(n_atomic_layers)
                                 for ii in indices_for_layers[jj] 
                                 if jj in args.fix_atomic_layers])
        slab.set_constraint(C)


    slab.center(args.vacuum / 2, axis='xyz'.index(args.ivacuum))
    org_chem_symbols = np.array(slab.get_chemical_symbols())
    org_atom_index = np.arange(len(slab), dtype=int)

    # Sort first by z-coordinates then by y and x.
    rpos = np.round(slab.positions, 4)
    new_atom_index = np.lexsort(
            (rpos[:, 0], rpos[:, 1], rpos[:, 2])
            )
    slab = slab[new_atom_index]

    # New order of chemical symbols
    if args.new_sym_order:
        assert set(args.new_sym_order) == set(org_chem_symbols)
        chem_sym_order = args.new_sym_order
    else:
        # Just stick to the original order.
        chem_sym_order = []
        for ss in primitive_cell.get_chemical_symbols():
            if not ss in chem_sym_order:
                chem_sym_order.append(ss)

    # Re-arrange the atoms according to the new order of chemical symbols
    new_atom_index = [ii for ss in chem_sym_order
                         for ii in org_atom_index[slab.symbols == ss]]
    slab = slab[new_atom_index]

    if args.out is None:
        outF = 'out_{}{}{}_{}.vasp'.format(
                args.hkl[0], args.hkl[1], args.hkl[2], args.nlayer
                )
        write(outF, slab, vasp5=True)
    else:
        outF = args.out
        write(outF, slab)

if __name__ == '__main__':
    cut_slab(sys.argv[1:])
