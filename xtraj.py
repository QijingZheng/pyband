#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, argparse
import numpy as np
from ase.io import read, write

def ase2axsf(trajs, ofile='traj.axsf'):
    '''
    Save the trajectory to the Xcrysden variable-cell-axsf format.
    '''
    nsteps = len(trajs)

    with open(ofile, 'w+') as out:
        out.write("ANIMSTEPS {}\nCRYSTAL\n".format(nsteps))

        for ii in range(nsteps):
            snap = trajs[ii]
            natoms = snap.get_number_of_atoms()

            if ii == 0:
                out.write('PRIMVEC\n')
                out.write(
                    '\n'.join([''.join(["%20.16f" % xx for xx in row])
                                                       for row in snap.cell])
        )
            out.write("\nPRIMCOORD {}\n".format(ii+1))
            out.write("{} 1\n".format(natoms))
            np.savetxt(out, np.c_[snap.get_atomic_numbers(), snap.positions],
                       fmt='%5d %22.16f %22.16f %22.16f'
                    )

def xdatcar2traj(cml):
    arg = parse_cml_args(cml)

    if arg.snaps:
        trajs = [read(f) for f in arg.snaps]
    else:
        trajs = read(arg.inputFile, index=':')

    if arg.outFmt == 'xyz':
        write('{}.xyz'.format(arg.outPrefix), trajs)
    elif arg.outFmt == 'pdb':
        write('{}.pdb'.format(arg.outPrefix), trajs)
    else:
        ase2axsf(trajs, ofile='{}.axsf'.format(arg.outPrefix))

def parse_cml_args(cml):
    '''
    CML parser.
    '''
    arg = argparse.ArgumentParser(add_help=True)

    arg.add_argument('-i', dest='inputFile', action='store', type=str,
                     default='XDATCAR',
                     help='The input file, either OUTCAR or XDATCAR.')
    arg.add_argument('-f', dest='outFmt', action='store', type=str,
                     default='axsf', choices=['xyz', 'pdb', 'axsf'],
                     help='The format of the trajectory file.')
    arg.add_argument('-o', dest='outPrefix', action='store', type=str,
                     default='traj',
                     help='The prefix of the output file.')
    arg.add_argument('-l', dest='snaps', action='store', type=str,
                     default=[], nargs='+',
                     help='List of structure files.')

    return arg.parse_args(cml)

if __name__ == '__main__':
    xdatcar2traj(sys.argv[1:])
