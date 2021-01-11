#!/usr/bin/env python3
'''
A script to plot work function using LOCPOT.
Requirement: python3, numpy, matplotlib, ase
Author: @Ionizing
Date: 22:46, Jan 11th, 2021.
'''

from argparse import ArgumentParser
import numpy as np
import matplotlib.pyplot as plt
from ase.calculators.vasp.vasp import VaspChargeDensity


def locpot_mean(fname="LOCPOT", axis='z', savefile='locpot.dat'):
    '''
    Reads the LOCPOT file and calculate the average potential along `axis`.
     @in: See function argument.
    @out:
          - xvals: grid data along selected axis;
          - mean: averaged potential corresponding to `xvals`.
    '''
    locd = VaspChargeDensity(fname)
    cell = locd.atoms[0].cell
    latlens = np.linalg.norm(cell, axis=1)
    vol = np.linalg.det(cell)

    iaxis = ['x', 'y', 'z'].index(axis.lower())
    axes = [0, 1, 2]
    axes.remove(iaxis)
    axes = tuple(axes)

    locpot = locd.chg[0]
    # must multiply with cell volume, similar to CHGCAR
    mean = np.mean(locpot, axes) * vol

    xvals = np.linspace(0, latlens[iaxis], locpot.shape[iaxis])

    # save to 'locpot.dat'
    np.savetxt(savefile, np.c_[xvals, mean],
               fmt='%13.5f', header='Distance(A) Potential(eV)')
    return (xvals, mean)


def parse_cml_arguments():
    parser = ArgumentParser(
        description='A tool to plot work function according to LOCPOT', add_help=True)
    parser.add_argument('-a', '--axis', type=str, action='store',
                        help='Which axis to be calculated: x, y or z. Default by z', default='z', choices=['x', 'y', 'z'])
    parser.add_argument('input', nargs='?', type=str,
                        help='The input file name, default by LOCPOT', default='LOCPOT')
    parser.add_argument('-w', '--write', type=str, action='store',
                        help='Save raw work function data to file, default by locpot.dat', default='locpot.dat')
    parser.add_argument('-o', '--output', type=str, action='store',
                        help='Output image file name, default by Workfunction.png', default='Workfunction.png')
    parser.add_argument('--dpi', type=int, action='store',
                        help='DPI of output image, default by 400', default=400)
    parser.add_argument('--title', type=str, action='store',
                        help='Title in output image. If none, no title is added, default is None', default=None)
    return parser.parse_args()
    pass


if '__main__' == __name__:
    args = parse_cml_arguments()
    x, y = locpot_mean(args.input, args.axis, args.write)

    plt.plot(x, y, color='k')
    plt.xlabel('Distance(A)')
    plt.ylabel('Potential(eV)')
    plt.grid(color='gray', ls='-.')
    plt.xlim(0, np.max(x))
    plt.ylim(np.min(y)-0.5, np.max(y)+0.5)

    if args.title:
        plt.title(args.title)
    plt.savefig(args.output, dpi=args.dpi)
