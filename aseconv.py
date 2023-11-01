#!/usr/bin/env python

from os.path import splitext
from ase.io import read, write
from optparse import OptionParser

############################################################
def command_line_arg():
    usage = "USAGE: %prog -i input_format -o output_format InputFiles"  
    par = OptionParser(usage=usage)

    par.add_option('-i', '--inp',
            action='store', type="string",
            dest='in_format', default=None,
            help='Format of input file.')

    par.add_option('-o', '--out',
            action='store', type="string",
            dest='out_format', default='vasp',
            help='Format of output file.')

    par.add_option('-p', '--prefix',
            action='store', type="string",
            dest='prefix', default=None,
            help='Prefix for output files.')

    return  par.parse_args( )


############################################################
if __name__ == '__main__':
    opts, args = command_line_arg()

    if opts.out_format == 'vasp':
        keyArgs = {'vasp5': True, 'direct': True}
    else:
        keyArgs = {}

    for ii, inF in enumerate(args):
        geo = read(inF, format=opts.in_format)
        if opts.prefix:
            outF = "{:s}_{:02d}.{:s}".format(opts.prefix, ii, opts.out_format)
        else:
            outF = "{:s}.{:s}".format(splitext(inF)[0], opts.out_format)

        write(outF, geo, **keyArgs)
