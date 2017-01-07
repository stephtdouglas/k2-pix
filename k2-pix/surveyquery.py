import numpy as np
import scipy as sp
from astropy import coordinates as coords
from k2phot import tpf_io
from astropy.io import fits
from astropy.wcs import WCS

def k2_ConvertHeaderWCS(tpf_file):
	"""
	Reads the WCS for the raw pixel counts into astropy.wcs
	and returns the WCS.
	Thanks Geert Barentsen, AAS229 Hack Together Day 
	"""
	f = fits.open(tpf_file)
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
		mywcs[newkey] = f[1].header[oldkey] 
	return WCS(mywcs).wcs_pix2world(1, 1, 1)

# Function to Get RA & Dec for the Object
def k2_getRADec(tpf_file):
	table, times, pixels, maskmap, maskheader, kpmag = tpf_io.get_data(tpf_file)
	SC_pos = coords.SkyCoord(maskheader["RA_OBJ"], maskheader["DEC_OBJ"], unit=u.deg)
	return SC_pos

# Get 2Mass K-Band Image
    twomass_images, pix_2mass, hdr_2mass = None, None, None
    try:
        twomass_images = SkyView.get_images(position=pos, survey=['DSS2 Red'],pixels=500)
        # print(twomass_images)
        if len(twomass_images)>0:
            pix_2mass = twomass_images[0][0].data
            hdr_2mass = twomass_images[0][0].header
    except astroquery.exceptions.TimeoutError,urllib.HTTPError:
        pix_2mass, hdr_2mass = None, None

# Get SDSS R-Band Image


# Function to Cut Survey Image


