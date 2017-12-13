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

from astropy.io import fits
from astropy import log
from astropy import visualization
from astropy.wcs import WCS
import astropy.units as u
from astropy import coordinates as coords

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
    for oldkey, newkey in funny_keywords.items():
        mywcs[newkey] = tpf_header[oldkey]

    return WCS(mywcs)

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
        if not self.verbose:
            import warnings
            warnings.filterwarnings("ignore")

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
    def position(self):
        return coords.SkyCoord(self.ra,self.dec,unit=u.deg,frame="icrs")

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
