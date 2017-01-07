from __future__ import division, print_function, absolute_import

import imageio
import numpy as np
from tqdm import tqdm
import warnings

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.image import imsave
import matplotlib.patheffects as path_effects
from matplotlib.colors import NoNorm

from astropy import log
from astropy import visualization
from astropy.wcs import WCS

import surveyquery

# Figure class
class K2Fig(object):
    """Figure showing K2 target pixel stamp and sky survey image."""

    def __init__(self,TPF):
        self.TPF = TPF
        self.verbose = self.TPF.verbose

    def cut_levels(self, min_percent=1., max_percent=95., data_col='FLUX'):
            """Determine the cut levels for contrast stretching.

            Returns
            -------
            vmin, vmax : float, float
                Min and max cut levels.
            """

            # Get co-added flux
            # update to use TPF
            sample = self.TPF.flux_binned()

            # Scale image
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore', message="(.*)invalid value(.*)")
                vmin, vmax = np.percentile(sample[sample > 0],
                                           [min_percent, max_percent])
            return vmin, vmax



    # Set up the figure and axes using astropy WCS
    def create_figure(self, output_filename, survey, stretch='log', vmin=1, vmax=None, min_percent=1, max_percent=95,
                      cmap='gray', contour_color='red', data_col='FLUX'):
        """Returns a matplotlib Figure object that visualizes a frame.

        Parameters
        ----------

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
        # Update to use TPF
        flx = self.TPF.flux_binned()
        # print(np.shape(flx))

        # calculate cut_levels
        if vmax is None:
            vmin, vmax = self.cut_levels(min_percent,max_percent,data_col)

        # Determine the figsize
        shape = list(flx.shape)
        # print(shape)
        # Create the figure and display the flux image using matshow
        fig = plt.figure(figsize=shape)
        # Display the image using matshow

        # Update to generate axes using WCS axes instead of plain axes
        ax = plt.subplot(projection=self.TPF.wcs)
        ax.set_xlabel('RA')
        ax.set_ylabel('Dec')

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

        current_ylims = ax.get_ylim()
        current_xlims = ax.get_xlim()

        pixels, header = surveyquery.getSVImg(self.TPF.position, survey)
        levels = np.linspace(np.min(pixels),np.percentile(pixels,95),10)
        ax.contour(pixels,transform=ax.get_transform(WCS(header)),
                    levels=levels,color=contour_color)

        ax.set_xlim(current_xlims)
        ax.set_ylim(current_ylims)

        fig.canvas.draw()
        plt.savefig(output_filename, bbox_inches='tight', dpi=300)
        return fig
