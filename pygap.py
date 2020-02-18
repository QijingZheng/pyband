#!/usr/bin/env python

import os, sys
import numpy as np
from optparse import OptionParser

############################################################
def get_bandinfo_from_outcar(inf='OUTCAR'):
    '''
    extract band energies from OUTCAR.
    '''
    outcar = [line for line in open(inf) if line.strip()]

    # nkpts = nband = ispin = Lvkpts = ibasis = Efermi = LineEfermi = None
    for ii, line in enumerate(outcar):
        if 'NKPTS =' in line:
            nkpts = int(line.split()[3])
            nband = int(line.split()[-1])

        if 'ISPIN  =' in line:
            ispin = int(line.split()[2])

        if "k-points in reciprocal lattice and weights" in line:
            Lvkpts = ii + 1

        if 'reciprocal lattice vectors' in line:
            ibasis = ii + 1

        if 'E-fermi' in line:
            Efermi = float(line.split()[2])
            LineEfermi = ii + 1
            # break

    # assert nkpts and nband and ispin and Lvkpts and ibasis and Efermi and LineEfermi, \
    #         "File seems not to be an OUTCAR format..."

    # k-points vectors and weights
    tmp = np.array([line.split() for line in outcar[Lvkpts:Lvkpts+nkpts]],
                   dtype=float)
    vkpts = tmp[:,:3]

    # for ispin = 2, there are two extra lines "spin component..."
    N = (nband + 2) * nkpts * ispin + (ispin - 1) * 2
    bands = []
    # vkpts = []
    for line in outcar[LineEfermi:LineEfermi + N]:
        if 'spin component' in line or 'band No.' in line:
            continue
        if 'k-point' in line:
            # vkpts += [line.split()[3:]]
            continue
        bands.append(float(line.split()[1]))

    bands = np.array(bands, dtype=float).reshape((ispin, nkpts, nband))
    vkpts = np.array(vkpts)

    return Efermi, bands, vkpts

def find_band_info(inf='OUTCAR', ratio=0.2, zero=None):
    '''
    Find the band information, e.g. VBM and CBM indexes etc.
    '''

    efermi, bands, vkpts = get_bandinfo_from_outcar(inf)
    if zero is not None:
        efermi = zero
    nspin, nkpts, nbands = bands.shape

    band_index = np.arange(nbands, dtype=int)

    band_energy_max = np.max(bands, axis=1)
    band_energy_min = np.min(bands, axis=1)

    fermi_cross_band = (band_energy_min < efermi) & (efermi < band_energy_max)

    band_info = []
    sys_info = {"NKPTS":nkpts, "NBANDS":nbands, "NSPIN":nspin, "Efermi":efermi}

    ivbm = icbm = -1
    for ii in range(nspin):
        # Fermi level does NOT cross any band, definitely a semiconductor.
        if not np.any(fermi_cross_band[ii]):
            bmax = band_energy_max[ii]
            bmin = band_energy_min[ii]

            s1 = (bmax[:-1] < efermi) & (efermi < bmax[1:])
            s2 = (bmin[:-1] < efermi) & (efermi < bmin[1:])
            
            ivbm_1 = list(s1).index(True)
            ivbm_2 = list(s2).index(True)
            
            assert ivbm_1 == ivbm_2
            ivbm = ivbm_1
            icbm = ivbm_1 + 1

        # Fermi level cross a few bands, maybe dopedsemiconductor or metal
        else:
            # find out the bands that cross the Fermi level
            xband_index = band_index[fermi_cross_band[ii]]

            e_xband_max = bands[ii,:,xband_index].max()
            e_xband_min = bands[ii,:,xband_index].min()
            e_xband_rng = e_xband_max - e_xband_min

            # the relative postion of Fermi level in the bands
            fermi_pos   = (efermi - e_xband_min) / e_xband_rng

            # Fermi level is at the band edges
            if (fermi_pos < ratio):
                # Fermi level near CBM
                icbm = xband_index.min()
                ivbm = icbm - 1
            elif (fermi_pos > (1-ratio)):
                # Fermi level near VBM
                ivbm = xband_index.max()
                icbm = ivbm + 1
            else:
                # metal
                pass

        if (icbm >= 0) and (ivbm >= 0):
            evbm = bands[ii,:,ivbm].max()
            vkpt_index = np.argsort(bands[ii,:,ivbm])[-1]
            kvbm = vkpts[vkpt_index]

            ecbm = bands[ii,:,icbm].min()
            ckpt_index = np.argsort(bands[ii,:,icbm])[0]
            kcbm = vkpts[ckpt_index]

            if vkpt_index == ckpt_index:
                which_gap = 'Direct_Gap'
            else:
                which_gap = 'inDirect_Gap'
        else:
            evbm = 0; ecbm = 0
            kcbm = kvbm = [0, 0, 0]
            which_gap = "Metal"

        # print ivbm, icbm, evbm, ecbm, kvbm, kcbm
        band_info.append(
                dict((("IVBM", ivbm + 1), ("ICBM", icbm + 1),
                      ("EVBM", evbm), ("ECBM", ecbm),
                      ("VBM_KPT_IND", vkpt_index), ("CBM_KPT_IND", ckpt_index),
                      ("KVBM", kvbm), ("KCBM", kcbm),
                      ("GAP",  ecbm - evbm),
                      ('NOTE', which_gap)
                      )))
        # print band_info[ii]

    format_band_info(sys_info, band_info)

############################################################
def format_band_info(sys_info, band_info):
    '''
    Output the band information.
    '''
    nspin = len(band_info)
    label = ['IND', 'ENG', 'KPT']
    lines  = ''

    lines += '{:6s} = {:4d}; '.format("NSPIN", sys_info['NSPIN'])
    lines += '{:6s} = {:4d};\n'.format("NKPTS", sys_info['NKPTS'])
    lines += '{:6s} = {:4d}; '.format("NBANDS",sys_info['NBANDS'])
    lines += '{:6s} = {:8.4f}\n'.format("Efermi", sys_info['Efermi'])

    if nspin == 2:
        lines += "-" * 54 + '\n'
        lines += '{:^10s}{:^22s}{:^22s}'.format('', 'SPIN_UP', 'SPIN_DN') + '\n'
        lines += "-" * 54 + '\n'
    else:
        lines += "-" * 32 + '\n'
    
    for band_label in ['CBM', 'VBM']:
        for ii, prefix in enumerate(['i', 'e', 'k']):
            if ii == 1:
                lines += '{:^5s}'.format(band_label)
            else:
                lines += '{:^5s}'.format('')

            lines += '{:^5s}'.format(label[ii])
            for ispin in range(nspin):
                info = band_info[ispin]

                k = (prefix + band_label).upper()
                if k.startswith('K'):
                    tmp = '{:6.4f} {:6.4f} {:6.4f}'.format(info[k][0], info[k][1], info[k][2])
                    tmp = '{:^22s}'.format(tmp)
                if k.startswith('I'):
                    tmp = '{:^22d}'.format(info[k])
                if k.startswith('E'):
                    tmp = '{:^22.5f}'.format(info[k])
                lines += tmp
            lines += '\n'

        if nspin == 2:
            lines += ' ' * 10 + '_' * 44 + '\n'
        else:
            lines += ' ' * 10 + '_' * 22 + '\n'

    lines +=  "{:^10s}".format('GAP')
    for ispin in range(nspin):
        lines += "{:^22.5f}".format(band_info[ispin]["GAP"])
    lines += '\n'

    lines +=  "{:^10s}".format('')
    for ispin in range(nspin):
        lines += "{:^22s}".format(band_info[ispin]["NOTE"])
    lines += '\n'

    if nspin == 2:
        if ('_Gap' in band_info[0]["NOTE"]) and ('_Gap' in band_info[1]["NOTE"]):
            vbm_erg_spin = np.array([xx["EVBM"] for xx in band_info])
            cbm_erg_spin = np.array([xx["ECBM"] for xx in band_info])
            vbm_kpt_spin = np.array([xx["VBM_KPT_IND"] for xx in band_info])
            cbm_kpt_spin = np.array([xx["CBM_KPT_IND"] for xx in band_info])

            vsort = np.argsort(vbm_erg_spin)
            csort = np.argsort(cbm_erg_spin)

            vbm_spin_erg_max = vbm_erg_spin[vsort[-1]]
            cbm_spin_erg_min = cbm_erg_spin[csort[0]]
            vbm_spin_kpt_ind = vbm_kpt_spin[vsort[-1]]
            cbm_spin_kpt_ind = cbm_kpt_spin[csort[0]]

            total_gap = cbm_spin_erg_min - vbm_spin_erg_max
            if vbm_spin_kpt_ind == cbm_spin_kpt_ind:
                note = 'Direct_Gap'
            else:
                note = 'inDirect_Gap'

            # lines += ' ' * 10 + '_' * 44 + '\n'
            lines += ' ' * 10 + "{:^44.5f}".format(total_gap) + '\n'
            lines += ' ' * 10 + "{:^44s}".format(note) + '\n'

    if nspin == 2:
        lines += "-" * 54
    else:
        lines += "-" * 32

    print lines

def command_line_arg():
    usage = "usage: %prog [options] OUTCAR1 OUTCAR2..."  
    par = OptionParser(usage=usage)

    par.add_option('-r', '--ratio',
            action='store', type="float",
            dest='ratio', default=0.2,
            help='')
    par.add_option('-z', '--zero',
            action='store', type="float",
            dest='zero', default=None,
            help='')

    return  par.parse_args( )

############################################################
if __name__ == '__main__':
    opts, args = command_line_arg()

    if (len(args) == 0):
        if os.path.isfile('OUTCAR'):
            args.append('OUTCAR')

    for inf in args:
        if os.path.isfile(inf):
            print inf, "->"
            find_band_info(inf, opts.ratio, opts.zero)
