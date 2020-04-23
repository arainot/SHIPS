############################
# Date: 06/04/2020
# Title: PSF-fitting script for SHIPS
# Author: J. Bodensteiner (2019). Edits: A. Rainot (2020)
# Description: Use this script to extract a spectrum of sources in IRDIS images with the SHIPS pipeline.
# VIP version: 0.9.11 (Rainot edit.)
# Python version: 3 ONLY
############################

# Load libraries
import lmfit
import pandas as pd
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.stats import sigma_clipped_stats
from astropy.visualization import simple_norm
from astropy.modeling.fitting import LevMarLSQFitter
from photutils.psf import photometry, DAOGroup, EPSFBuilder
from photutils import CircularAperture, CircularAnnulus, aperture_photometry
from astropy.nddata import NDData
from photutils.psf import extract_stars
import matplotlib.pyplot as plt
from matplotlib.patches import Circle


##############################################################################
# 0. class and function definition
##############################################################################
class Star:
    def __init__(self, xcoord, ycoord, star_id, mag='None'):
        self.xcoord = xcoord  # pixel x coordinate
        self.ycoord = ycoord  # pixel y coordinate
        self.star_id = star_id  # star id
        self.magnitude = mag  # magnitude from input list


############################################################################
# write a fits file containing a 2D image
def write_2Dimage(header, image, outfilename):
    hdul_new = fits.HDUList()
    hdul_new.append(fits.PrimaryHDU(data=image, header=header))
    hdul_new.writeto(outfilename)

# round odd PSF size to the lower odd number
def round_down_to_odd(f):
    return np.ceil(f) // 2 * 2 - 1

# Gaussian function
def gaussian(x, height, center, std):
    return (height * np.exp(-1 * (x - center)**2 / (2*std**2)))


# fit a Gaussian
def f_single(params, wavelength, fluxes):
    h = params['h']
    std = params['std']
    cen = params['cen']
    g = gaussian(wavelength, h, cen, std)
    error = (fluxes - g)**2
    return error


# get PSF from the image of the central star
def prep_psf(psf_data, x_star, y_star, n_resample=6, plot=False):
    # create bkg subtracted cutouts around stars used for creation of ePSF
    stars_tbl = Table()
    stars_tbl['x'] = [x_star]  # coordinates of central star
    stars_tbl['y'] = [y_star]

    # subtract a median background from the PSF image
    mean_val, median_val, std_val = sigma_clipped_stats(psf_data, sigma=2.)
    psf_data -= median_val

    # extract 61 x 61 px cutouts around the stars => whole image
    stars = extract_stars(NDData(data=psf_data), stars_tbl, size=61)

    # initialize EPSF Builder object with required oversampling factor
    n_resample = 6
    epsf_builder = EPSFBuilder(oversampling=n_resample, maxiters=10,
                               progress_bar=False)
    # build the actual ePSF from the cutouts
    epsf, fitted_stars = epsf_builder(stars)

    # plot cuts through the ePSF at different x positions
    len_x = len(epsf.data[0, :])
    x = np.linspace(0, 100, len_x)

    params = lmfit.Parameters()
    params.add('h', 0.03, min=0.025, max=0.045, vary=True)
    params.add('std', 5, min=1, max=20, vary=True)
    params.add('cen', 50, min=47, max=53, vary=True)

    xvals = [int(len(x)/4), int(len(x)/2.1), int(len(x)/2)]
    cutthrough = epsf.data[int(len(x)/2), :]

    minimizer = lmfit.Minimizer(f_single, params, fcn_args=(x, cutthrough))
    result = minimizer.minimize()

    h, gauss_std = result.params['h'], result.params['std']
    cen = result.params['cen']

    if plot is True:
        # plot the star (if its only one)
        fig_stars, ax_stars = plt.subplots()
        norm = simple_norm(stars[0], 'log', percent=99.)
        colbar_in = ax_stars.imshow(stars[0], origin='lower', cmap='inferno',
                                    norm=norm)
        fig_stars.colorbar(colbar_in)
        fig_stars.suptitle('Star used to fit the PSF')
        # plt.show()

        fig_psfima, ax_psf_ima = plt.subplots()
        colbar = ax_psf_ima.imshow(epsf.data, origin='lower', cmap='inferno')
        fig_psfima.suptitle("Fitted ePSF")
        fig_psfima.colorbar(colbar)

        fig_psf, ax_psf = plt.subplots()

        for xval in xvals:
            cutthrough = epsf.data[xval, :]
            lab = 'x = ' + str(xval)
            ax_psf.plot(x, cutthrough, label=lab)

        g = gaussian(x, h, cen, gauss_std)
        ax_psf.plot(x, g, label='Gaussian fit')

        ax_psf.set_xlabel('y coordinate')
        ax_psf.set_ylabel('ePSF')
        ax_psf.legend()

        plt.show()

    return epsf, gauss_std

# Define the PSF-fitting function
def psf_fitting(array, input_psf, fwhm, comp_pos, psf_pos, full_output=False):

    """ Applies PSF fitting to given sources in order to retrieve positions & flux parameters with errors.

    Parameters
    ----------
    array : array_like
        Input 2D frame.
    input_psf : array_like
        2d array psf.
    fwhm : float
        FWHM of the PSF.
    comp_pos : float or array_like
        Positions (X,Y) of companions.
    psf_pos : float or array_like
        Positions (X,Y) of the center of the PSF.
    full_output : True or False
        Plots the fit on positions and the residual images.
    Returns
    -------
    x_comp : array_like
        Best-fit X positions for every input companion.
    xerr_comp : array_like
        Best-fit X positional errors for every input companion.
    y_comp : array_like
        Best-fit Y positions for every input companion.
    yerr_comp : array_like
        Best-fit Y positional errors for every input companion.
    f_comp : array_like
        Best-fit flux photometry for every input companion.
    ferr_comp : array_like
        Best-fit flux photometry errors for every input companion.
    if full_output = True
        contrast_mag : array_like
            Contrast magnitude of sources
        contrast_mag : array_like
            Contrast magnitude errors of sources
    """

    pos_comp = np.array(comp_pos) # Create a usuable array for the companion positions
    x_centralstar, y_centralstar = psf_pos[0], psf_pos[1]  # pixel
    psf_data = input_psf
    epsf, gauss_std = prep_psf(psf_data, x_centralstar, y_centralstar, plot=False)

    ############################################################################
    # definitions for the photometry
    ############################################################################
    aper_rad = 1.*fwhm
    fshape = int(round_down_to_odd(len(input_psf[0])))
    phot = photometry.BasicPSFPhotometry(group_maker=DAOGroup(15.),
                                         psf_model=epsf,
                                         bkg_estimator=None,
                                         fitter=LevMarLSQFitter(),
                                         fitshape=(fshape),
                                         aperture_radius=aper_rad)


    ############################################################################
    # Measure the flux of the central star
    ############################################################################
    pos = Table(names=['x_0', 'y_0'], data=[[x_centralstar], [y_centralstar]])
    result_tab = phot.do_photometry(image=psf_data, init_guesses=pos)
    flux_star, fluxerr_star = result_tab['flux_fit'][0], result_tab['flux_unc'][0]

    print("######################################################################")
    print("The guess flux is : %2.4f" % result_tab['flux_0'])
    print("The measured flux of the central star is: % 2.4f + /- % 2.4f"
          % (flux_star, fluxerr_star))
    print("######################################################################")
    residual_image = phot.get_residual_image()

    # plot the input image and the residual image for the central star
    fig_res = plt.figure(figsize=(8, 8))

    ax1 = fig_res.add_subplot(1, 2, 1)
    norm = simple_norm(psf_data, 'sqrt', percent=95.)
    ax1.imshow(psf_data, cmap='viridis', aspect=1, origin='lower', norm=norm)
    ax1.set_title('Input Image')
    # mark input position of central star
    e = Circle(xy=(x_centralstar, y_centralstar), radius=aper_rad)
    e.set_facecolor('none')
    e.set_edgecolor(color='k')
    ax1.add_artist(e)

    ax2 = fig_res.add_subplot(1, 2, 2)
    ax2.imshow(residual_image, cmap='viridis', aspect=1, origin='lower')
    ax2.set_title('Residual Image')
    # mark fitted position of central star
    e = Circle(xy=(result_tab['x_fit'], result_tab['y_fit']), radius=aper_rad)
    e.set_facecolor('none')
    e.set_edgecolor(color='crimson')
    ax2.add_artist(e)
    plt.show()

    ############################################################################
    ############################################################################
    ############################################################################
    # Load in the actual data
    ############################################################################
    image = array

    ##############################################################################
    # Coordinates of companions
    # xpos, ypos, ids = table['x'], table['y'], table['star_id']
    xpos = pos_comp[:,0]
    ypos = pos_comp[:,1]
    ids = np.linspace(1,len(pos_comp),len(pos_comp))
    # make a list of star objects of class star from input file
    starlist = list()
    for i in range(1, len(xpos)):
        s = Star(xpos[i], ypos[i], ids[i])
        starlist.append(s)

    # plot them all on the image
    cols = ['crimson', 'cyan', 'C1', 'deepskyblue', 'limegreen', 'magenta',
            'limegreen', 'red', 'cyan', 'C1', 'deepskyblue',
            'crimson', 'cyan', 'C1', 'deepskyblue', 'limegreen', 'magenta',
            'limegreen', 'red', 'cyan', 'C1', 'deepskyblue']

    fig, ax = plt.subplots(figsize=(8, 8))
    norm = simple_norm(image, 'sqrt', percent=95.)
    im1 = ax.imshow(image, aspect=1, origin='lower', cmap='Greys', norm=norm)

    for s in range(len(starlist)):
        star = starlist[s]
        e = Circle(xy=(star.xcoord, star.ycoord), radius=2.5)
        e.set_facecolor('none')
        e.set_edgecolor(color=cols[s],)
        ax.add_artist(e)
        ax.annotate(str(star.star_id), (star.xcoord, star.ycoord), (3, 10),
                    color=cols[s], textcoords='offset points', fontsize=12)
    plt.xlabel("x [px]", fontsize=13)
    plt.ylabel("y [px]", fontsize=13)
    plt.show()

    # define plotting window and
    wind = 25
    rad = aper_rad

    # define arrays to store values
    x_comp = np.zeros_like(xpos)
    xerr_comp = np.zeros_like(xpos)
    y_comp = np.zeros_like(xpos)
    yerr_comp = np.zeros_like(xpos)
    f_comp = np.zeros_like(xpos)
    ferr_comp = np.zeros_like(xpos)
    c_comp = np.zeros_like(xpos) # Contrast magnitudes
    cerr_comp = np.zeros_like(xpos) # Contrast magnitude errors

    # LOOP OVER STARS
    for i in range(len(xpos)-1):
        star = starlist[i]
        x_pos, y_pos = [], []
        x_pos.append(star.xcoord)
        y_pos.append(star.ycoord)

        # center of the frames
        cy, cx = round(image.shape[0]/2), round(image.shape[1]/2)
        # star-center distances
        xy = (star.xcoord-cx, star.ycoord-cy)
        # radial distance from center
        r_pos = np.sqrt(xy[0]**2 + xy[1]**2)

        annulus = CircularAnnulus((star.xcoord, star.ycoord), 15, 25)
        phot_annulus = aperture_photometry(image, annulus, method='subpixel',
                                           subpixels=5)
        mean_bkg = phot_annulus['aperture_sum'] / annulus.area

        # aperture = CircularAperture((star.xcoord, star.ycoord), fwhm/2.)
        # annulus = CircularAnnulus((cx, cy), r_pos-fwhm, r_pos+fwhm)
        # phot_aperture = aperture_photometry(image, aperture, method='subpixel',
        #                                     subpixels=5)
        # phot_annulus = aperture_photometry(image, annulus, method='subpixel',
        #                                    subpixels=5)
        # mean_background = (phot_annulus['aperture_sum'] -
        #                    phot_aperture['aperture_sum'])/(annulus.area -
        #                                                    aperture.area)
        # print(mean_background)

        # prepare the image => correct local background around the source
        for h in range(image.shape[0]):
            for j in range(image.shape[1]):
                image[h][j] -= mean_bkg[0]

        pos = Table(names=['x_0', 'y_0'], data=[x_pos, y_pos])
        result_tab = phot.do_photometry(image=image, init_guesses=pos)
        residual_image = phot.get_residual_image()

        x_comp[i] = result_tab['x_fit'][0]
        xerr_comp[i] = result_tab['x_0_unc'][0]
        y_comp[i] = result_tab['y_fit'][0]
        yerr_comp[i] = result_tab['y_0_unc'][0]
        f_comp[i] = result_tab['flux_fit'][0]
        ferr_comp[i] = result_tab['flux_unc'][0]

        flux_comp = result_tab['flux_fit'][0]
        fluxerr_comp = result_tab['flux_unc'][0]

        flux_ratio = flux_comp / flux_star
        flux_ratio_err = flux_ratio * ((fluxerr_comp / flux_comp)**2 +
                                       (fluxerr_star / flux_star)**2)**0.5

        contrast_mag = -2.5 * np.log10(flux_ratio)
        contrast_mag_err = 2.5 / np.log(10) * ((fluxerr_comp / flux_comp)**2 +
                                               (fluxerr_star / flux_star)**2)**0.5
        print("-----------------------------------------------------------------")
        print("Star " + str(star.star_id))
        print("The measured xpos is: %2.4f +/- %2.4f" % (x_comp[i], xerr_comp[i]))
        print("The measured ypos is: %2.4f +/- %2.4f" % (y_comp[i], yerr_comp[i]))
        print("The measured flux is: %2.4f +/- %2.4f" % (flux_comp, fluxerr_comp))

        if full_output == True:

            print("The contrast magnitude is: %2.4f +/- %2.4f" % (contrast_mag, contrast_mag_err))
            c_comp[i] = contrast_mag
            cerr_comp[i] = contrast_mag_err
            # plot the result
            fig = plt.figure(figsize=(8, 8))

            ax = fig.add_subplot(1, 2, 1)
            norm = simple_norm(image, 'sqrt', percent=95.)
            ax.imshow(image, cmap='viridis', aspect=1, origin='lower', norm=norm)
            ax.set_title('Star ' + str(star.star_id))
            e = Circle(xy=(star.xcoord, star.ycoord), radius=aper_rad)
            e.set_facecolor('none')
            e.set_edgecolor(color='k')
            ax.add_artist(e)

            ax.set_xlim([star.xcoord+wind, star.xcoord-wind])
            ax.set_ylim([star.ycoord+wind, star.ycoord-wind])
            ax.set_xlabel('x coordinate [px]')
            ax.set_ylabel('y coordinate [px]')

            ax2 = fig.add_subplot(1, 2, 2)
            norm = simple_norm(image, 'sqrt', percent=95.)

            ax2.imshow(residual_image, cmap='viridis', aspect=1,
                       origin='lower', norm=norm)
            e = Circle(xy=(result_tab['x_fit'], result_tab['y_fit']), radius=aper_rad)
            e.set_facecolor('none')
            e.set_edgecolor(color='crimson')
            ax2.add_artist(e)

            e = Circle(xy=(star.xcoord, star.ycoord), radius=aper_rad)
            e.set_facecolor('none')
            e.set_edgecolor(color='k')
            ax2.add_artist(e)

            ax2.set_xlim([star.xcoord+wind, star.xcoord-wind])
            ax2.set_ylim([star.ycoord+wind, star.ycoord-wind])

            ax2.set_title('Residual Image')
            ax2.set_xlabel('x coordinate [px]')
            ax2.set_ylabel('y coordinate [px]')

            plt.show()
        print("-----------------------------------------------------------------")
    if full_output == True:
        return x_comp,xerr_comp,y_comp,yerr_comp,f_comp,ferr_comp,flux_star,fluxerr_star,c_comp,cerr_comp
    else:
        return x_comp,xerr_comp,y_comp,yerr_comp,f_comp,ferr_comp,flux_star,fluxerr_star
