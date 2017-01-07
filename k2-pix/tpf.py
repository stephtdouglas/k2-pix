#!/usr/bin/env python
# -*- coding: utf-8 -*-
""""Create movies or animated gifs from Kepler Target Pixel Files.

Author: Geert Barentsen, Stephanie Douglas, Erin Maier
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

__all__ = ["TargetPixelFile"]

import argparse
import imageio
import numpy as np
from tqdm import tqdm
import warnings

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pl
from matplotlib.image import imsave
import matplotlib.patheffects as path_effects
from matplotlib.colors import NoNorm

from astropy.io import fits
from astropy.time import Time
from astropy import log
from astropy import visualization
from astropy.wcs import WCS

###
# Kepler constants
###

# Quality flags, from the Kepler Archive Manual
KEPLER_QUALITY_FLAGS = {
    "1": "Attitude tweak",
    "2": "Safe mode",
    "4": "Coarse point",
    "8": "Earth point",
    "16": "Zero crossing",
    "32": "Desaturation event",
    "64": "Argabrightening",
    "128": "Cosmic ray",
    "256": "Manual exclude",
    "1024": "Sudden sensitivity dropout",
    "2048": "Impulsive outlier",
    "4096": "Argabrightening",
    "8192": "Cosmic ray",
    "16384": "Detector anomaly",
    "32768": "No fine point",
    "65536": "No data",
    "131072": "Rolling band",
    "262144": "Rolling band",
    "524288": "Possible thruster firing",
    "1048576": "Thruster firing"
}

###
# Fix TPF header WCS keywords
###

def k2_ConvertHeaderWCS(tpf_header):
    """
    Reads the WCS for the raw pixel counts into astropy.wcs
    and returns the WCS.
    """
    funny_keywords = {'1CTYP4': 'CTYPE1',
                      '2CTYP4': 'CTYPE2',
                      '1CRPX4': 'CRPIX1',
                      '2CRPX4': 'CRPIX2',
                      '1CRVL4': 'CRVAL1',
                      '2CRVL4': 'CRVAL2',
                      '1CUNI4': 'CUNIT1',
                      '2CUNI4': 'CUNIT2',
                      '1CDLT4': 'CDELT1',
                      '2CDLT4': 'CDELT2',
                      '11PC4': 'PC1_1',
                      '12PC4': 'PC1_2',
                      '21PC4': 'PC2_1',
                      '22PC4': 'PC2_2'}
    mywcs = {}
    for oldkey, newkey in funny_keywords.iteritems():
        mywcs[newkey] = tpf_header[oldkey]

    return WCS(mywcs).wcs_pix2world(1, 1, 1)

###
# Helper classes
###

class BadKeplerFrame(Exception):
    """Raised if a frame is empty."""
    pass


class BadCadenceRange(Exception):
    """Raised if a the range of frames is invalid."""
    pass


class TargetPixelFile(object):
    """Represent a Target Pixel File (TPC) from the Kepler spacecraft.

    Parameters
    ----------
    filename : str
        Path of the pixel file.

    cache : boolean
        If the file name is a URL, this specifies whether or not to save the
        file locally in Astropy’s download cache. (default: True)

    verbose : boolean
        If True, print debugging output.
    """
    def __init__(self, filename, cache=True, verbose=False):
        self.filename = filename
        self.hdulist = fits.open(filename, cache=cache)
        self.header = self.hdulist[1].header
        self.no_frames = len(self.hdulist[1].data['FLUX'])
        self.verbose = verbose

    @property
    def objectname(self):
        try:
            return self.hdulist[0].header['OBJECT']
        except KeyError:
            return ''

    @property
    def ra(self):
        return self.hdulist[0].header['RA_OBJ']

    @property
    def dec(self):
        return self.hdulist[0].header['DEC_OBJ']

    @property
    def wcs(self):
        return k2_ConvertHeaderWCS(self.header)

    @property
    def flux(self):
        flx = np.array(self.hdulist[1].data["FLUX"])

        bad_frames = np.zeros(len(flx),bool)

        # delete bad frames
        for frameno in np.arange(self.no_frames):
            flags = self.quality_flags(frameno)
            if len(flags)>0:
                bad_frames[frameno] = True

            infinites = np.isfinite(flx[frameno])
            fin_med = np.nanmedian(flx[frameno][infinites==False])
            flx[frameno][infinites==False] = fin_med

        good_flux = np.delete(flx, np.where(bad_frames)[0], axis=0)
        return good_flux

    def flux_binned(self):
        flux_series = self.flux
        return np.sum(flux_series,axis=0)

    def quality_flags(self, frameno):
        """Returns a list of strings describing the quality flags raised."""
        quality = self.hdulist[1].data['QUALITY'][frameno]
        flags = []
        for flag in KEPLER_QUALITY_FLAGS.keys():
            if quality & int(flag) > 0:
                flags.append(KEPLER_QUALITY_FLAGS[flag])
        return flags


    def cut_levels(self, min_percent=1., max_percent=95., data_col='FLUX'):
        """Determine the cut levels for contrast stretching.

        Returns
        -------
        vmin, vmax : float, float
            Min and max cut levels.
        """

        # Get co-added flux
        sample = self.flux_binned()

        # Scale image
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', message="(.*)invalid value(.*)")
            vmin, vmax = np.percentile(sample[sample > 0],
                                       [min_percent, max_percent])
        return vmin, vmax

    def create_figure(self, stretch='log', vmin=1, vmax=None,
                      cmap='gray', data_col='FLUX'):
        """Returns a matplotlib Figure object that visualizes a frame.

        Parameters
        ----------
        frameno : int
            Image number in the target pixel file.

        vmin : float, optional
            Minimum cut level (default: 0).

        vmax : float, optional
            Maximum cut level (default: 5000).

        cmap : str, optional
            The matplotlib color map name.  The default is 'gray',
            can also be e.g. 'gist_heat'.

        raw : boolean, optional
            If `True`, show the raw pixel counts rather than
            the calibrated flux. Default: `False`.

        Returns
        -------
        image : array
            An array of unisgned integers of shape (x, y, 3),
            representing an RBG colour image x px wide and y px high.
        """
        # Get the flux data to visualize
        flx = self.flux_binned()
        # print(np.shape(flx))

        # calculate cut_levels
        if vmax is None:
            vmin, vmax = self.cut_levels()

        # Determine the figsize
        shape = list(flx.shape)
        # print(shape)
        # Create the figureand display the flux image using matshow
        fig = pl.figure(figsize=shape)
        # Display the image using matshow
        ax = fig.add_subplot(1, 1, 1)
        if self.verbose:
            print('{} vmin/vmax = {}/{} (median={})'.format(data_col, vmin, vmax, np.nanmedian(flx)))

        if stretch == 'linear':
            stretch_fn = visualization.LinearStretch()
        elif stretch == 'sqrt':
            stretch_fn = visualization.SqrtStretch()
        elif stretch == 'power':
            stretch_fn = visualization.PowerStretch(1.0)
        elif stretch == 'log':
            stretch_fn = visualization.LogStretch()
        elif stretch == 'asinh':
            stretch_fn = visualization.AsinhStretch(0.1)
        else:
            raise ValueError('Unknown stretch: {0}.'.format(stretch))

        transform = (stretch_fn +
                     visualization.ManualInterval(vmin=vmin, vmax=vmax))
        ax.imshow((255*transform(flx)).astype(int), aspect='auto',
                   origin='lower', interpolation='nearest',
                   cmap=cmap, norm=NoNorm())
        ax.set_xticks([])
        ax.set_yticks([])
        ax.axis('off')
        fig.subplots_adjust(left=0.0, right=1.0, top=1.0, bottom=0.0)
        fig.canvas.draw()
        return fig
