""" Plot a K2 TPF and overlay a sky survey image."""

from __future__ import division, print_function, absolute_import

import argparse
import matplotlib.pyplot as plt

from tpf import TargetPixelFile

def k2pix():
    """ Script to plot a K2 TPF and overlay a sky survey image."""
    parser = argparse.ArgumentParser(
        description="Plots a co-added Target Pixel File (TPF) from NASA's "
                    "Kepler/K2/TESS with a survey image overlaid.")
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='')
    parser.add_argument('--output', metavar='FILENAME',
                        type=str, default=None,
                        help='.gif or .mp4 output filename (default: gif with the same name as the input file)')
    parser.add_argument('--stretch', type=str, default='log',
                        help='type of contrast stretching: "linear", "sqrt", '
                             '"power", "log", or "asinh" (default is "log")')
    parser.add_argument('--min_cut', type=float, default=None,
                        help='minimum cut level (default: use min_percent)')
    parser.add_argument('--max_cut', type=float, default=None,
                        help='maximum cut level (default: use max_percent)')
    parser.add_argument('--min_percent', metavar='%', type=float, default=1.,
                        help='minimum cut percentile (default: 1.0)')
    parser.add_argument('--max_percent', metavar='%', type=float, default=95.,
                        help='maximum cut percentile (default: 95)')
    parser.add_argument('--cmap', type=str,
                        default='gray', help='matplotlib color map name '
                                             '(default: gray)')
    parser.add_argument('tpf_filename', nargs='+',
                        help='path to one or more Target Pixel Files (TPF)')

    datagroup = parser.add_mutually_exclusive_group()
    datagroup.add_argument('--raw', action='store_true',
                           help="show the uncalibrated pixel counts ('RAW_CNTS')")
    datagroup.add_argument('--background', action='store_true',
                           help="show the background flux ('FLUX_BKG')")
    datagroup.add_argument('--cosmic', action='store_true',
                           help="show the cosmic rays ('COSMIC_RAYS')")

    args = parser.parse_args(args)

    # # What is the data column to show?
    # if args.raw:
    #     data_col = 'RAW_CNTS'
    # elif args.background:
    #     data_col = 'FLUX_BKG'
    # elif args.cosmic:
    #     data_col = 'COSMIC_RAYS'
    # else:
    #     data_col = 'FLUX'

    for fn in args.tpf_filename:
        try:
            tpf = TargetPixelFile(fn, verbose=args.verbose)
            # and here goes the call to actually make the figure
        except Exception as e:
            if args.verbose:
                raise e
            else:
                print('ERROR: {}'.format(e))

# Example use
if __name__ == '__main__':
    # fn = ('http://archive.stsci.edu/missions/kepler/target_pixel_files/'
    #       '0007/000757076/kplr000757076-2010174085026_lpd-targ.fits.gz')
    fn = ("https://archive.stsci.edu/missions/k2/target_pixel_files/c5/"
          "211900000/72000/ktwo211972478-c05_lpd-targ.fits.gz")
    tpf = TargetPixelFile(fn)
    # Then run with all default arguments to generate image
    tpf.create_figure()
    plt.show()
