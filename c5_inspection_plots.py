
# coding: utf-8

# In[16]:

from __future__ import division, print_function, absolute_import
import sys
import itertools
import time

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.ticker
import matplotlib.cm as cm
from matplotlib.font_manager import FontProperties
import matplotlib.gridspec as gridspec
from matplotlib import colors
import urllib2 as urllib

import pywcsgrid2
import astroquery
from astroquery.sdss import SDSS
from astroquery.skyview import SkyView
from astropy.wcs import WCS
import astropy.io.ascii as at
import astropy.io.fits as fits
import astropy.table as table
from astropy import coordinates as coords
import astropy.units as u

import K2fov.projection as proj
import K2fov.fov as fov
from K2fov.K2onSilicon import angSepVincenty,getRaDecRollFromFieldnum

import k2phot
from k2phot import centroid
from k2phot import tpf_io
from k2spin import prot
from hypra.utils import cat_match, cat_io, k2utils
from hypra.plot import color_mag
from hypra.plot.plot_k2 import plot_chips, setup_k2_axes
import convertmass

period_file_base = "../k2_tables/c5_k2sc_output_2016-05-17_"
nfiles = 14
# period_file_base = "../k2_tables/c5_k2sc_output_extra_2016-08-01_"
# nfiles = 5

cmap = plt.cm.viridis_r

for i in range(nfiles):
    #print(i+1)
    period_file = "{0}{1}.csv".format(period_file_base,i+1)
    new_res = at.read(period_file)

    peak_file = "{0}{1}_allpeaks.csv".format(period_file_base,i+1)
    new_peaks = at.read(peak_file)
    if i==0:
        res = new_res
        peaks = new_peaks
    else:
        res = table.vstack([res,new_res])
        peaks = table.vstack([peaks,new_peaks])

# print(res.dtype)
c5_epic = res["EPIC"]

def setup_figure():

    fig = plt.figure(figsize=(14,16))

    gs_base = gridspec.GridSpec(1,2,width_ratios=[1,3])

    gs_sky = gridspec.GridSpecFromSubplotSpec(3,1,subplot_spec=gs_base[0])

    sky_axes = [plt.subplot(gs_sky[i]) for i in range(3)]

    lc_n = 7

    gs_lc = gridspec.GridSpecFromSubplotSpec(lc_n,1,subplot_spec=gs_base[1])

    lc_axes = [plt.subplot(gs_lc[i]) for i in range(lc_n)]

    return fig, sky_axes, lc_axes


# In[17]:

def k2sc_io(filename):
    """ Read in a K2SC light curve file, and return the time and flux.

    Flux is calculated as the completely detrended flux plus the
    time-dependent trend (i.e., only the position-dependent trend is removed),
    and then normalized.

    Inputs
    ------
    filename: string
        a valid K2SC file, including full or relative path

    Returns:
    --------
    time, flux, white: arrays

    """
    with fits.open(filename) as hdu:
        #print(hdu.info())

        tabl = hdu[1].data
        good = np.isfinite(tabl["flux"]) & (np.isfinite(tabl["trend_t"]))
        med_trend = np.median(tabl["trend_t"][good])
        time = tabl["time"][good]
        flux = tabl["flux"][good] + tabl["trend_t"][good] - med_trend

        white = tabl["flux"][good]

    return time,flux,white


# In[35]:

def plot_phased(ax, t, f, period, power,color):

    phased_t = t % period / period

    ax.plot(phased_t,f,'.',color=color,label="Prot={0:.2f}, Power={1:.2f}".format(float(period),float(power)))
    ax.set_xlim(0,1)
    ax.set_yticklabels([])
    ax.legend(loc="best",numpoints=1,borderaxespad=0)


# In[73]:

def plot_lcs(axes, EPIC):

    lc_file = "/home/stephanie/data/c5_k2sc/hlsp_k2sc_k2_llc_{0}-c05_kepler_v1_lc.fits".format(EPIC)
    t, f, w = k2sc_io(lc_file)

    i = np.where(res["EPIC"]==EPIC)[0]
    if len(i>0):
        i = i[0]

    ylims = np.percentile(f[np.isfinite(f)],[0.05,99.05])

    color1 = cmap(0.35)#plt.cm.inferno(0.5)
    color2 = cmap(0.85)#plt.cm.inferno(0.75)
    color3 = cmap(0.6)
    label_fontsize = "small"

    # Periodogram
    ls_out = prot.run_ls(t,f,np.ones_like(f),0.1,prot_lims=[0.1,70],run_bootstrap=False)
    periods, pgram = ls_out[2], ls_out[3]

    if res[i]["num_sig"]>=2:
        peak_locs = np.where(peaks["EPIC"]==EPIC)[0]
        if len(peak_locs)>2:
            sort_locs = peak_locs[np.argsort(peaks["power"][peak_locs])]
            # print(peaks["period"][sort_locs], peaks["power"][sort_locs])
            third_period = peaks["period"][sort_locs[-3]]
            third_power = peaks["power"][sort_locs[-3]]
    else:
        third_period, third_power = None, None

    axes[0].plot(periods,pgram,'k-')
    axes[0].set_xlim(0.1,70)
    axes[0].set_xscale("log")
    axes[0].plot(res[i]["sig_period"],res[i]["sig_power"]*1.1,'v',
                mfc=color1,mec="none",ms=11)
    axes[0].plot(res[i]["sec_period"],res[i]["sec_power"]*1.1,'v',
                mfc=color2,mec="none",ms=11)
    if third_period is not None:
        axes[0].plot(third_period, third_power*1.1,'v',
                    mfc=color3,mec="none",ms=11)
    if res[i]["sig_power"]>0:
        ymax = res[i]["sig_power"]*1.15
    else:
        ymax = max(pgram)
    axes[0].set_ylim((0,ymax))
    axes[0].axhline(res[i]["threshold"],color='grey',ls="-.",
                    label="threshold {0:2f}".format(float(res[i]["threshold"])))
    #axes[0].set_xlabel("Period (d)",y=0.95)
    #axes[0].tick_params(labelbottom=False,labeltop=True)
    axes[0].set_xticklabels(["","0.1","1","10"])
    axes[0].set_ylabel("Power",fontsize=label_fontsize)
    axes[0].set_xlabel("Period (d)",fontsize=label_fontsize)

    # Full light curve
    axes[1].plot(t,f,'k.')
    axes[1].set_ylim(ylims)
    axes[1].set_xlim(t[0],t[-1])
    axes[1].set_ylabel("Light curve",fontsize=label_fontsize)
    axes[1].set_yticklabels([])
    # axes[1].tick_params(labelbottom=False)
    axes[1].set_xlabel("Time (d)",fontsize=label_fontsize)

    # White noise
    axes[2].plot(t,w,'k.')
    axes[2].set_ylim(np.percentile(w,[0.5,99.5]))
    axes[2].set_xlim(t[0],t[-1])
    axes[2].set_ylabel("White noise",fontsize=label_fontsize)
    axes[2].set_yticklabels([])
    # axes[2].tick_params(labelbottom=False)
    axes[2].set_xlabel("Time (d)",fontsize=label_fontsize)

    # Time-dependent trend
    axes[3].plot(t,f - w,'k.')
    axes[3].set_ylabel("Time-dep",fontsize=label_fontsize)
    axes[3].set_xlim(t[0],t[-1])
    axes[3].set_yticklabels([])
    # axes[3].tick_params(labelbottom=False)
    axes[3].set_xlabel("Time (d)",fontsize=label_fontsize)


    # Phase folded 1
    if res[i]["sig_period"]>0:
        plot_phased(axes[4],t,f,res[i]["sig_period"],res[i]["sig_power"],color=color1)
        axes[4].set_ylim(ylims)
        axes[4].set_ylabel("P1",fontsize=label_fontsize)

        if res[i]["sig_period"]>=2:
            repeats = np.arange(t[0],t[-1],res[i]["sig_period"])
            for r in repeats:
                axes[1].axvline(r,ls="--",color=color1,lw=2)
        # axes[4].tick_params(labelbottom=False)
        axes[4].set_xlabel("Phase",fontsize=label_fontsize)
    else:
        axes[4].set_axis_off()


    # Phase folded 2
    if res[i]["sec_period"]>0:
        plot_phased(axes[5],t,f,res[i]["sec_period"],res[i]["sec_power"],color=color2)
        axes[5].set_ylim(ylims)
        axes[5].set_ylabel("P2",fontsize=label_fontsize)

        # print(res[i]["sec_period"])
        if res[i]["sec_period"]>=2:
            repeats = np.arange(t[0],t[-1],res[i]["sec_period"])
            # print(repeats)
            for r in repeats:
                # print(r)
                axes[1].axvline(r,ls=":",color=color2,lw=2)
        axes[5].set_xlabel("Phase",fontsize=label_fontsize)
    else:
        # axes[5].set_yticks([])
        # axes[5].set_yticklabels([])
        axes[5].set_axis_off()

    if third_period is not None:
        plot_phased(axes[6],t,f,third_period,third_power,color=color3)
        axes[6].set_ylim(ylims)
        axes[6].set_ylabel("P3",fontsize=label_fontsize)
        axes[6].set_xlim(0,1)
        axes[6].set_xlabel("Phase",fontsize=label_fontsize)
        #
        # if third_period]>=2:
        #     repeats = np.arange(t[0],t[-1],third_period)
        #     # print(repeats)
        #     for r in repeats:
        #         # print(r)
        #         axes[1].axvline(r,ls="-.",color=color3,lw=2)
    else:
        # axes[6].set_yticks([])
        # axes[6].set_yticklabels([])
        axes[6].set_axis_off()


    plt.subplots_adjust(hspace=0.6)


# In[56]:

def stamp(img, maskmap, ax=None, cmap="cubehelix"):
    """Plot a single pixel stamp."""

    if ax is None:
        fig = plt.figure(figsize=(8,8))
        ax = plt.subplot(111)

    ax.matshow(img, cmap=cmap, origin="lower", norm=colors.LogNorm())
    ax.set_xlim(-0.5,maskmap.shape[1]-0.5)
    ax.set_ylim(-0.5,maskmap.shape[0]-0.5)
    # I remain unconvinced that these are labelled right...
    # but the coordinates are plotting right, and this matches the DSS images
    # (except East and West are flipped)!
    ax.set_ylabel("Y")
    ax.set_xlabel("X")

    return ax


# In[57]:

def extract_wcs(dataheader):
    """
    Genderate WCS for a K2 TPF, since astropy has decided not to cooperate.

    dataheader is the header of extension #1 (zero-indexing!) of a K2 TPF
    """
    w4 = WCS(naxis=2)

    w4.wcs.crpix = [dataheader["1CRPX4"],dataheader["2CRPX4"]]
    w4.wcs.cdelt = np.array([dataheader["1CDLT4"],dataheader["2CDLT4"]])
    w4.wcs.crval = [dataheader["1CRVL4"],dataheader["2CRVL4"]]
    w4.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w4.wcs.pc = [[dataheader["11PC4"],dataheader["12PC4"]],
                 [dataheader["21PC4"],dataheader["22PC4"]]]
    return w4


# In[69]:

def plot_sky(axes, EPIC):


    # Get K2 image
    tpf_file = "/home/stephanie/data/c5_tpf/ktwo{0}-c05_lpd-targ.fits.gz".format(EPIC)
    tabl, times, pixels, maskmap, maskheader, kpmag = tpf_io.get_data(tpf_file)
    coadd = np.sum(pixels,axis=0)
    pos = coords.SkyCoord(maskheader["RA_OBJ"], maskheader["DEC_OBJ"],
                          unit=u.deg)

    with fits.open(tpf_file) as hdu2:
        dataheader = hdu2[1].header

#     keysel = np.empty(13, "S6")
#     keysel[:] = "binary"
    w2 = extract_wcs(dataheader)

    # Get 2MASS K image
    twomass_images, pix_2mass, hdr_2mass = None, None, None
    try:
        twomass_images = SkyView.get_images(position=pos, survey=['DSS2 Red'],pixels=500)
        # print(twomass_images)
        if len(twomass_images)>0:
            pix_2mass = twomass_images[0][0].data
            hdr_2mass = twomass_images[0][0].header
    except astroquery.exceptions.TimeoutError,urllib.HTTPError:
        pix_2mass, hdr_2mass = None, None

    # Get SDSS r image
    pix,hdr,im = None,None,None
    xid = SDSS.query_region(pos)
    if xid is not None:
        one_match = table.Table(xid[0])
        im = SDSS.get_images(matches=one_match, band='r')
        pix = im[0][0].data
        hdr = im[0][0].header

    # Plot K2 image
    # Set up the GridHelper to merge the axes
    grid_helper = pywcsgrid2.GridHelper(wcs=w2)

    # Plot the pixel stamp as usual, except with the WCS
    ax1 = pywcsgrid2.subplot(441, grid_helper=grid_helper,
                             aspect=1, adjustable="box-forced")
    ax1.matshow(coadd, origin="lower", cmap=cmap, norm=colors.LogNorm())
    if pix is not None:
        median = np.median(pix)
        stdev = np.std(pix)
        levels = np.linspace(median + stdev, np.max(pix), 20)
        #grey = plt.cm.Greys(0.3)
        ax1[hdr].contour(pix,cmap=plt.cm.Greys, levels=levels)

        # Plot the SDSS image rotated into the same frame as the pixel stamp
        ax2 = pywcsgrid2.subplot(445, grid_helper=grid_helper,
                                 aspect=1, adjustable="box-forced",
                                 sharex=ax1, sharey=ax1)
        vmax = np.percentile(pix,99.9)
        ax2[hdr].imshow_affine(pix, origin="lower", cmap=cmap,#'Greys',
                               norm=colors.LogNorm(vmax=vmax,clip=True))
        median2 = np.median(coadd)
        stdev2 = np.std(coadd)
        levels2 = np.linspace(median, np.max(coadd), 5)

    if pix_2mass is not None:

        # Matplotlib was choking on negative values when it saved the figure
        pos_min = np.min(pix_2mass[pix_2mass>0])
        vmin = np.max([pos_min,np.percentile(pix_2mass,0.01)])
        vmax = np.percentile(pix_2mass,99.9)
        pix_2mass.flags.writeable=True
        pix_2mass[pix_2mass<=0] = pos_min

        # Plot the 2MASS image rotated into the same frame as the pixel stamp
        ax4 = pywcsgrid2.subplot(449, grid_helper=grid_helper,
                                 aspect=1, adjustable="box-forced",
                                 sharex=ax1, sharey=ax1)
        ax4[hdr_2mass].imshow_affine(pix_2mass, origin="lower",
                                    cmap=cmap#'Greys'
                                    ,norm=colors.LogNorm(vmax=vmax,
                                    vmin=vmin,clip=True))


    # Plot chip image - not working?
    ax3 = plt.subplot(4,4,13)
    ax3.plot(pos.ra.value,pos.dec.value,'*',color=cmap(0.9),ms=24)
    setup_k2_axes(ax3, extents=[125.5,135,15.5,24])
    plot_chips(ax3, 5)
    # ax3.set_yticks(np.arange(8,25,4))
    ax3.set_yticks(np.arange(16,25,2))
    # ax3.set_xticks(np.arange(122,139,4))
    # ax3.set_xticks(np.arange(122,139,2),minor=True)
    ax3.set_xlabel(r"$\alpha_{2000}$",fontsize="large")
    ax3.set_ylabel(r"$\delta_{2000}$",fontsize="large")

    plt.subplots_adjust(hspace=0.15)

    del(im,pix,hdr)
    del(twomass_images, pix_2mass, hdr_2mass)

if __name__=="__main__":

    # fig, sky_axes, lc_axes = setup_figure()
    # plot_sky(sky_axes,211897926)
    # plot_lcs(lc_axes,211897926)
    # plt.tight_layout()
    #
    # plt.savefig("../k2_plots/test_inspection_211897926.png")

    print(len(c5_epic),len(np.unique(c5_epic)))


    for i,ep in enumerate(c5_epic):

#         if i<=419:
#             continue

        ep = 211972478
        i = np.where(c5_epic==ep)[0]

        print(i,ep)
        fig, sky_axes, lc_axes = setup_figure()
        plot_sky(sky_axes,ep)
        plot_lcs(lc_axes,ep)
        plt.suptitle("EPIC {0}".format(ep))
        # plt.tight_layout()

        plt.savefig("../k2_plots/k2sc_inspect/inspection_{0}.png".format(ep))
        # plt.savefig("../k2_plots/test_inspection_{0}.png".format(ep))
        plt.close("all")

        # Pause 1 second in an attempt to reduce server time-out
        time.sleep(1)

        if i>=1:
            break
