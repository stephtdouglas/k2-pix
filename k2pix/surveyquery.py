import numpy as np
import astroquery
from astropy.io import fits
from astropy.wcs import WCS
from astroquery.skyview import SkyView
import matplotlib.pyplot as plt
import urllib

# Function to Get a Survey Image Provided the Object Position & Survey Name
# May need a downgrade of HTML5Lib: pip install --upgrade html5lib==1.0b8
def getSVImg(SC_ObjPos, SkyViewSurvey):
    """
    Input: a SkyCoord Position Vector and a Survey Type that is Compatable with NASA SkyView
    Output: A Numpy Pixel Array and Header Array
    """
    img_survey, pix_survey, hdr_survey = None, None, None
    try:
        img_survey = SkyView.get_images(position=SC_ObjPos, survey=[SkyViewSurvey],
                                        coordinates="J2000", pixels=30)
        if len(img_survey) > 0:
            pix_survey = img_survey[0][0].data
            hdr_survey = img_survey[0][0].header
    except (astroquery.exceptions.TimeoutError, urllib.error.HTTPError):
        pix_survey, hdr_survey = None, None
    return pix_survey, hdr_survey

if __name__ == "__main__":
    # Demo to Showcase that the SkyView image is working.
    Survey = 'DSS'
    ObjPos = 'Dumbbell Nebula'
    pixels, header = getSVImg(ObjPos, Survey)
    plt.imshow(pixels)
    plt.show()
