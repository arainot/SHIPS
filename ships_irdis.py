############################
# Date: 08/08/2019
# Title: Running script for SHIPS for IRDIS data
# Description: Use this script to run SHIPS for IRDIS data. In this script you'll find all the necessary parameters to run SHIPS. ONLY SPHERE-DC DATA FOR NOW. VIP is used.
# VIP version: 0.9.11 (Rainot edit.)
# Python version: 3 ONLY
############################

# Set up your parameters

## Define images to analyse
cube_filepath = '/Users/alan/Documents/PhD/Data/SPHERE/IRDIS/QZCar/ird_convert_dc-IRD_SCIENCE_REDUCED_MASTER_CUBE-center_im.fits'
#cube_filepath = '/Users/alan/Desktop/cube_free_E.fits'
wavelength_filepath = '/Users/alan/Documents/PhD/Data/SPHERE/IRDIS/QZCar/ird_convert_dc-IRD_SCIENCE_LAMBDA_INFO-lam.fits'
angles_filepath = '/Users/alan/Documents/PhD/Data/SPHERE/IRDIS/QZCar/ird_convert_dc-IRD_SCIENCE_PARA_ROTATION_CUBE-rotnth.fits'
psf_filepath = '/Users/alan/Documents/PhD/Data/SPHERE/IRDIS/QZCar/ird_convert_recenter_dc5-IRD_SCIENCE_PSF_MASTER_CUBE-median_unsat.fits'
# wavelength_filepath = '/Users/alan/Documents/PhD/Data/SPHERE/Hugues_data/IRDIS/CEN3/ird_convert_recenter_dc-IRD_SCIENCE_LAMBDA_INFO-lam.fits'
# cube_filepath = '/Users/alan/Documents/PhD/Data/SPHERE/Hugues_data/IRDIS/CEN3/ird_convert_recenter_dc-IRD_SCIENCE_REDUCED_MASTER_CUBE-center_im.fits'
# angles_filepath = '/Users/alan/Documents/PhD/Data/SPHERE/Hugues_data/IRDIS/CEN3/ird_convert_recenter_dc-IRD_SCIENCE_PARA_ROTATION_CUBE-rotnth.fits'
# psf_filepath = '/Users/alan/Documents/PhD/Data/SPHERE/Hugues_data/IRDIS/CEN3/ird_convert_recenter_dc-IRD_SCIENCE_PSF_MASTER_CUBE-median_unsat.fits'

## Photometry
comp_pos = ([490,455],[559,579],[687,629],[569,701],[310,485],[695,418],[295,471],[653,706],[517,242],[226,568],[345,776],[846,490],[897,506],[70,532],[818,179],[451,59],[48,397],[730,44],[647,15]) # Companion position in pixels (X,Y)
#comp_pos = ([492,456],[559,579],[688,630],[570,702],[311,486],[696,419],[296,472],[654,707],[517,241],[227,569],[344,775],[847,491],[899,507],[72,533],[819,180],[451,60],[49,396],[732,44],[648,16]) # Companion position in pixels (X,Y)
psf_pos = (32, 32) # PSF position in pixels (X,Y)
radial_dist = [60.16643583,81.60882305,210.78899402,197.121377136,201.68291946,211,219.,240.63665556,269.09106265,290.11,313.16,334.845228417,384.760748992,442.00180288,451.91383567,456.697517309,477.277760379,515.37656136072,515.7130985344468] # Radial distance of companion in pixels
position_angle = [295.9,311.8] # Position angle of companion in degrees
noise_aperture_pos_comp = (512,512) # Position in pixels of the circular annulus aperture for noise measurement in the case of the companion
noise_aperture_pos_psf = (12,22) # Position in pixels of the circular annulus aperture for noise measurement in the case of the PSF

## Computing power
ncores = 4 # Number of cores you are willing to share for the computation

## Do you want to see the image?
see_cube = False # Original cube
see_collapsed_cube = False # Collapsed cube
see_psf_norm = False # Normalised PSF
see_cube_centre = False # Check if the image is centered correctly
size_psf = 31
## PCA
ncomp_pca = 1 # Number of principal components for PCA
opti_pca = False # Optimise the number of PCA components?
source = (501,525) # Source where to optimise the PCA

## SNR maps
snr_maps = False # Would you like to make and save an SNR map to disk?
snr_map_file = '/home/alan/data/Backup_macbook/SPHERE/IRDIS/QZCar/SNRmap_VIP.fits' # Finish the file with .fits

## Detection
adi_frame = False # Would you like to apply ADI on the frame?
adi_plot = False # Would you like to see the resulting plot?
adi_min_scale = -1 # Minimum colour scale for the ADI plot
adi_max_scale = 3 # Maximum colour scale for the ADI plot
detection = False # Would you like the algorithm to detect sources for you? !! WARNING: this is a simple detection !!
detect_sigma = 5 # What sigma limit would you like for the detection?

## Contrast curves
contrast_curves = False # True or False !! computationally intensive !!
n_branches = 1 # Number of branches for contrast curves

## Photometric errors of PSF
psf_errors = True # Compute the photometric errors of the central star's PSF
psf_errors_save = True # Save the errors to a file?
psf_errors_file = '/Users/alan/Documents/PhD/Data/SPHERE/IRDIS/QZCar/PSF_errors.txt' # Filepath to save the PSF errors

## Spectrum extraction with Simplex Nelder-Mead optimisation
extract_spec = False # Will start the simplex Nelder-Mead optimisation for spectrum extraction
ann_width = 3 # Annulus width of Simplex
aper_radius = 3 # Aperture Radius of PCA
save_spec = False # Save the spectrum to ascii file
sspec_file_K1 = '/Users/alan/Documents/PhD/Data/SPHERE/IRDIS/QZCar/VIP_simplex_K1.txt' # Filepath to save the Simplex spectrum for the K1 band
sspec_file_K2 = '/Users/alan/Documents/PhD/Data/SPHERE/IRDIS/QZCar/VIP_simplex_K2.txt' # Filepath to save the Simplex spectrum for the K2 band

## Spectrum extraction with MCMC
extract_mcmc = False # Will compute the MCMC for all sources !! This takes ~22h per source and is very computer intensive !!
source = 'QZCar' # Give name for your primary star
mcmc_path = '/Users/alan/Documents/PhD/Data/SPHERE/IRDIS/QZCar/spectra/' # Directory where MCMC results will be stored

## Reading MCMC results
read_mcmc = False # Do you wish to read the MCMC results?
source = 'QZCar' # Give name for your source
mcmc_path = '/Users/alan/Documents/PhD/Data/SPHERE/IRDIS/QZCar/MCMC/' # Directory where MCMC results are stored

## Load calibrated FASTWIND models of the central star
fastwind = False # Use FASTWIND model spectra for the star
fastwind_path = '/Users/alan/Nextcloud/PhD/Thesis/SPHERE/spectra/fastwind/qzcarAa1/' # Directory where the FASTWIND flux are
rad_fast = 22.1 # Radius of model star
dist_fast = 100. # Distance to consider for the flux of the calibrated spectrum in Ro

## Compute calibrated spectrum of companion
calib_spec = False # Do you wish to calibrate the spectrum of the companions?
save_calib_spec = False # Would you like to save the calibrated spectrum & associated error?
calib_star_spec_path = '/Users/alan/Nextcloud/PhD/Thesis/SPHERE/spectra/fastwind/qzcar_fastwind_spec.txt' # Path to calibrated spectrum of central star
sspec_file = '/Users/alan/Documents/PhD/Data/SPHERE/IRDIS/QZCar/VIP_simplex.txt' # Path to spectrum file
cspec_file = '/Users/alan/Documents/PhD/Data/SPHERE/IRDIS/QZCar/VIP_calib_spectra.txt' # Path to calibrated spectrum

## Magnitude contrasts
mag_contr = False # Do you to calculate the magnitude contrasts for your sources?
print_mag_contr = False # Do you wish to print the magnitude contrasts to the screen?

## Absolute magnitude !! Work in Progress !!
abs_mag = False # Would you like to calculate the absolute magnitudes of your sources?
print_abs_mag = False # Do you wish to print the absolute magnitudes to the screen?
star_mag_Y = 5.75 # Magnitude of central star Y band
star_mag_J = 5.551 # Magnitude of central star J band
star_mag_H = 5.393 # Magnitude of central star H band
star_mag_V = 6.24 # Magnitude of central star V band
star_dist = 2300. # Distance to central star in parsec

# ---------------------------------------------------------------------------

# Running script (DO NOT MODIFY)

# Some definitions

## Load libraries
import __init__
import sys
import matplotlib
import vip_hci
from hciplot import plot_frames, plot_cubes
from vip_hci.metrics.fakecomp import cube_inject_companions, frame_inject_companion, normalize_psf
import numpy as np
import scipy
import astropy.io.fits as fits
from astropy import wcs
import photutils
from photutils import CircularAperture
from photutils import CircularAnnulus
import matplotlib.pyplot as plt
import glob
#import pyklip.instruments.SPHERE as SPHERE
import math as mh
#import pyklip.parallelized as parallelized
from scipy.integrate import quad, dblquad

## Define constants
c = 299792458. # Speed of light
Ro = 6.957e8 # Solar Radius
sr2pc = 44334448.0068964 # Convert steraradians to parsec
pxscale = 0.01225 # IRDIS pixel scale in arcsec/pixel
PA = np.array(position_angle) + 90 # Correct for VIP unconventional rotation axis

## Open image files
cube = vip_hci.fits.open_fits(cube_filepath)
wl = vip_hci.fits.open_fits(wavelength_filepath)
angs = vip_hci.fits.open_fits(angles_filepath)
psf = vip_hci.fits.open_fits(psf_filepath)

## Define some Parameters
psf = np.median(psf, axis=1) # Take the median value of all psf observations
psf_scaled = np.zeros_like(psf) # The psf will need to be scaled
flevel = np.zeros_like(cube[:,0,0,0]) # Flux level for the companion
flevel = np.array(flevel) # Redefinition - why?
centy, centx = vip_hci.var.frame_center(cube[0,0]) # Find the coordinate centers of the frame

## Check the RAW data cubes
if see_cube == True:
    ds9 = vip_hci.Ds9Window()
    ds9.display(cube[0,0])

## Get FWHM of images & normalised PSF
psf_med = vip_hci.preproc.cosmetics.cube_crop_frames(psf, size_psf, xy=(32, 32), verbose=True, force=True) # Resize the PSF
psf_norm, maxflux, fwhm = vip_hci.metrics.normalize_psf(psf_med, fwhm='fit',size=None, threshold=None,mask_core=None, model='gauss',imlib='opencv',interpolation='lanczos4',force_odd=True,full_output=True,verbose=False) # maxflux is a dummy variable

### Plot it
if see_psf_norm == True:
    plot_frames(psf_norm[1], grid=True, size_factor=4)

## Check if the cube is centred correctly by plotting
if see_cube_centre == True:
    plot_frames(vip_hci.preproc.frame_crop(cube[0,0], 50), grid=True, size_factor=4)

## Optimise the number of PCA components
if opti_pca == True:
    vip_hci.pca.pca(cube[0], angs, fwhm=fwhm[0], source_xy=(501,525),mask_center_px=None, ncomp=(1, 41, 2))
    sys.exit("PCA optimised. To continue, please input the PCA value in the script and skip this process.")

## Detection with VIP, for now only with the first wavelength
if adi_frame == True:
    ### Compute the ADI frame for all 39 wavelengths and 48 rotations
    print("Computing the ADI frame...")
    fr_adi = vip_hci.medsub.median_sub(cube, -angs, scale_list=wl, mode='fullfr') # 2D ADI frame
    print("Done!")
    ### Plot the frame
    if adi_plot == True:
        plot_frames(fr_adi, vmin=adi_min_scale, vmax=adi_max_scale)
    ### Compute the detection of sources
    if detection==True:
        detect = vip_hci.metrics.detection(fr_adi, fwhm=fwhm[0], psf=psf_norm[0], debug=False, mode='log', snr_thresh=detect_sigma,bkg_sigma=detect_sigma,matched_filter=True,vmin=adi_min_scale,vmax=adi_max_scale,verbose=False) # Sigma limit provided by user
        print("Detected sources : " , "\n", detect)
        detect_pos = np.array(detect) # Converted to array in order to be used later
        sys.exit("Sources detected. To continue, please input the target coordinates in the script and skip this process.")

## SNR maps
if snr_maps == True:
    snrmap = vip_hci.metrics.snrmap(vip_hci.pca.pca(cube, -angs, scale_list=wl, ncomp=ncomp_pca, verbose=True), fwhm[0], nproc=ncores, plot=True)
    vip_hci.fits.write_fits(snr_map_file,snrmap) # Write SNR maps to file
    sys.exit("SNR maps created. To continue, please input follow from the beginning process.")

# Stellar photometry of the companion

## Collapse the images for better photometry measurement
# cube_wl_coll = np.zeros_like(cube[:,0,:,:])
# for i in range(len(wl)):
#         cube_wl_coll[i] = vip_hci.hci_postproc.median_sub(cube[i],-angs,fwhm=fwhm[i],verbose=False) # Rotate & collapse along the rotation axis - 3D image
cube_derot = vip_hci.preproc.cube_derotate(cube,angs) # Rotate the images to the same north
cube_wl_coll = vip_hci.preproc.cube_collapse(cube_derot,wl_cube=True) # Collapse along the rotation axis - 3D image
#cube_coll = vip_hci.preproc.cube_collapse(cube_derot,wl_cube=False) # Collapse along the wavelength axis - 2D image

## Check the collapsed data cubes
if see_collapsed_cube == True:
    ds9 = vip_hci.Ds9Window()
    ds9.display(cube_wl_coll[0])#,cube_coll) # cube_wl_coll on the left and cube_coll on the right

## Aperture photometry of companions and PSF

### Define photometry
noise_phot = np.zeros_like(wl) #Noise photometry
psf_final_sum = np.zeros_like(wl) #PSF photometry
final_sum_K1 = np.zeros_like(radial_dist) #Companion photometry in the K1 band
final_sum_K2 = np.zeros_like(radial_dist) #Companion photometry in the K2 band

### Apertures
aper_noise_psf = photutils.CircularAnnulus(psf_pos,noise_aperture_pos_psf[0],noise_aperture_pos_psf[1])

### Aperture photometry - PSF
for i in range(0,len(wl)):
    ### Aperture
    aper_psf = photutils.CircularAperture(psf_pos, 1./2*fwhm[i])
    ### Flux
    phot_psf = photutils.aperture_photometry(psf[i], aper_psf)
    phot_psf_noise = photutils.aperture_photometry(psf[i], aper_noise_psf)
    psf_bkg_mean = phot_psf_noise['aperture_sum'] / aper_noise_psf.area()
    psf_bkg_sum = psf_bkg_mean * aper_psf.area()
    psf_final_sum[i] = phot_psf['aperture_sum'] - psf_bkg_sum

### Aperture photometry - Companions
for i in range(0,len(radial_dist)):
    ### Apertures dependent on companions
    aper_noise_comp = photutils.CircularAnnulus((512,512),radial_dist[i]-5,radial_dist[i]+5)
    aper_comp_K1 = photutils.CircularAperture((comp_pos[i][0], comp_pos[i][1]),1./2*fwhm[0])
    aper_comp_K2 = photutils.CircularAperture((comp_pos[i][0], comp_pos[i][1]),1./2*fwhm[1])
    ### Flux
    phot_noise_K1 = photutils.aperture_photometry(cube_wl_coll[0], aper_noise_comp)
    phot_noise_K2 = photutils.aperture_photometry(cube_wl_coll[1], aper_noise_comp)
    phot_K1 = photutils.aperture_photometry(cube_wl_coll[0], aper_comp_K1)
    phot_K2 = photutils.aperture_photometry(cube_wl_coll[1], aper_comp_K2)
    bkg_mean_K1 = (phot_noise_K1['aperture_sum']-phot_K1['aperture_sum']) / (aper_noise_comp.area()-aper_comp_K1.area())
    bkg_mean_K2 = (phot_noise_K2['aperture_sum']-phot_K2['aperture_sum']) / (aper_noise_comp.area()-aper_comp_K2.area())
    bkg_sum_K1 = bkg_mean_K1 * aper_comp_K1.area()
    bkg_sum_K2 = bkg_mean_K2 * aper_comp_K2.area()
    final_sum_K1[i] = phot_K1['aperture_sum'] - bkg_sum_K1
    final_sum_K2[i] = phot_K2['aperture_sum'] - bkg_sum_K2

### Scaling the PSF for normalisation -- SHOULD I JUST TAKE PSF_NORM INSTEAD?
psf_scaled = np.zeros_like(psf)
for i in range (0,len(psf)):
    psf_scaled[i] = psf[i]/psf_final_sum[i]

## Contrast curve
if contrast_curves == True:
    cube_negfc = vip_hci.metrics.cube_inject_companions(cube,psf_norm,-angs,flevel=-105,plsc=pxscale,rad_dists=[radial_dist],theta=PA) # Remove companion using NEGFC technique
    print("Companion removed")
    print("Computing contrast curve...")
    contrcurve = vip_hci.metrics.contrast_curve(cube_negfc,-angs,psf,np.average(fwhm),pxscale,psf_final_sum,vip_hci.pca.pca,nbranch=n_branches,
              dpi=300, student=False, debug=True ,plot=True, verbose=True, full_output=True, ncomp=ncomp_pca, scale_list=wl)


## PSF error calculation
if psf_errors == True:
    psferr = vip_hci.fits.open_fits(psf_filepath) # Open the raw PSFs
    stddev = np.zeros_like(wl) # Create an array for the stored standard deviation
    for i in range(len(wl)): # Loop over the wavelengths
        psferr_med = vip_hci.preproc.cosmetics.cube_crop_frames(psferr[i], size_psf, xy=(32, 32), verbose=True, force=True) # Resize the PSF
        psf_norm_err, maxflux_err, fwhm_err = vip_hci.metrics.normalize_psf(psferr_med, fwhm='fit',size=None, threshold=None,mask_core=None, model='gauss',imlib='opencv',interpolation='lanczos4',force_odd=True,full_output=True,verbose=False) # Measure the maximum flux for each PSF
        stddev[i] = np.std(maxflux_err,ddof=1) # Calculate the standard deviation for the PSFs
    if psf_errors_save: # Save the error
        np.savetxt(psf_errors_file, stddev, delimiter='   ') # Saves to file

# Spectrum extraction with NM
if extract_spec == True:

    ## Define some parameters
    f_guess_pl = max(final_sum_K1) # Flux first guess
    f_range_K1 = np.zeros((len(final_sum_K1),200))
    f_range_K2 = np.zeros((len(final_sum_K2),200))
    for i in range(0,len(final_sum_K1)):
        f_range_K1[i] = np.linspace(0.2*np.abs(final_sum_K1[i]), 30*np.abs(final_sum_K1[i]),200)
        f_range_K2[i] = np.linspace(0.2*np.abs(final_sum_K2[i]), 30*np.abs(final_sum_K2[i]),200)
    p_in = np.array([radial_dist,PA]) # Regroup companion positions
    simplex_options = {'xtol':1e-2, 'maxiter':500, 'maxfev':1000} # Set the simplex options
    simplex_guess_K1 = np.zeros((len(radial_dist),3)) # Set the simplex variable: r, PA, flux for every companion - K1
    simplex_guess_K2 = np.zeros((len(radial_dist),3)) # Set the simplex variable: r, PA, flux for every companion - K2
    ## Start Simplex
    for i in range(0,1):
        print("Companion index: ", i + 1) # Companions for IRDIS
        ### K1
        comp_xycoord = [[comp_pos[i][0],comp_pos[i][1]]] # Companion coords
        simplex_K1 = vip_hci.negfc.firstguess(cube[0],-angs,psf_norm[0],ncomp=ncomp_pca,plsc=pxscale,planets_xy_coord=comp_xycoord,fwhm=fwhm[0],annulus_width=ann_width,aperture_radius=aper_radius,simplex_options=simplex_options,f_range=f_range_K1[i],simplex=True,fmerit='sum',collapse='median',svd_mode='lapack',scaling=None,verbose=True,plot=False,save=False) # This takes some time
        x_simplex = simplex_K1[0] * np.cos(np.deg2rad(simplex_K1[1])) + centx
        y_simplex = simplex_K1[0] * np.sin(np.deg2rad(simplex_K1[1])) + centy
        comp_xycoord = [[x_simplex[0],y_simplex[0]]] # Companion coords
        simplex_guess_K1[i] = vip_hci.negfc.firstguess(cube[0],-angs,psf_norm[0],ncomp=ncomp_pca,plsc=pxscale,planets_xy_coord=comp_xycoord,fwhm=fwhm[0],annulus_width=ann_width,aperture_radius=aper_radius,simplex_options=simplex_options,f_range=f_range_K1[i],simplex=True,fmerit='sum',collapse='median',svd_mode='lapack',scaling=None,verbose=True,plot=False,save=False) # This takes some time
        ### K2
        comp_xycoord = [[comp_pos[i][0],comp_pos[i][1]]] # Companion coords
        simplex_K2 = vip_hci.negfc.firstguess(cube[1],-angs,psf_norm[1],ncomp=ncomp_pca,plsc=pxscale,planets_xy_coord=comp_xycoord,fwhm=fwhm[1],annulus_width=ann_width,aperture_radius=aper_radius,simplex_options=simplex_options,f_range=f_range_K2[i],simplex=True,fmerit='sum',collapse='median',svd_mode='lapack',scaling=None,verbose=False,plot=False,save=False) # This takes some time
        x_simplex = simplex_K2[0] * np.cos(np.deg2rad(simplex_K2[1])) + centx
        y_simplex = simplex_K2[0] * np.sin(np.deg2rad(simplex_K2[1])) + centy
        comp_xycoord = [[x_simplex[0],y_simplex[0]]] # Companion coords
        simplex_guess_K2[i] = vip_hci.negfc.firstguess(cube[1],-angs,psf_norm[1],ncomp=ncomp_pca,plsc=pxscale,planets_xy_coord=comp_xycoord,fwhm=fwhm[1],annulus_width=ann_width,aperture_radius=aper_radius,simplex_options=simplex_options,f_range=f_range_K2[i],simplex=True,fmerit='sum',collapse='median',svd_mode='lapack',scaling=None,verbose=False,plot=False,save=False) # This takes some time
        print("K1: ", simplex_guess_K1[i])
        print("K2: ", simplex_guess_K2[i])

## Save the spectrum
if save_spec == True:
    np.savetxt(sspec_file_K1, simplex_guess_K1, delimiter='   ') # Saves to file
    np.savetxt(sspec_file_K2, simplex_guess_K2, delimiter='   ')

# Spectrum extraction with MCMC
if extract_mcmc == True:
    instru= 'IRDIS36059' # Define instrument parameters
    ann_width=annulus_width # Annulus width of MCMC
    aperture_radius=aperture_width # Aperture radius
    fig_merit='sum' # Summation figure of merit
    outpath = mcmc_path.format(source) # Path to save MCMC files

    print("########## MCMC Sampling starting... ##########")

    nwalkers, itermin, itermax = (100,200,500) # as recommended by Oli W
    for i in range(len(final_sum)): # For each wavelength channel
        initialState = simplex_guess[i] # Take r, PA and flux from simplex
        bounds=[[0.75*initialState[0],1.25*initialState[0]],[0.75*initialState[1],1.25*initialState[1]],[0.75*initialState[2],1.30*initialState[2]]] # Initiate bounds
        output_file = source+'_IRDIS_companions_{}'.format(i) # Save to output file

        chain_40 = vip.negfc.mcmc_negfc_sampling(cube[i], -angs,  psf_scaled[i], ncomp_pca, pxscale, initialState, ann_width,
                                                 aperture_radius, cube_ref=None, svd_mode='lapack', nwalkers=nwalkers,
                                                 bounds=bounds, niteration_min=itermin,
                                                 niteration_limit=itermax, check_maxgap=50, nproc= ncores,
                                                 output_file=output_file, display=True,verbose=True, save=True,
                                                 rhat_threshold=1.01, niteration_supp=0, fmerit=fig_merit) # MCMC run per channel
    print("########## MCMC Sampling done! ##########")

## Read MCMC files
if read_mcmc == True:
    import pickle # Import important MCMC libraries
    from pickle import Pickler
    pickler={}
    mcmc_result={}
    outpath = mcmc_path.format(source) # Path to save MCMC files
    for i in range(0,3): # Read all channels and store them to variables
        with open(outpath+source+'K1_IRDIS_companion_S{}/MCMC_results'.format(i),'rb') as fi:
                pickler["myPickler{}".format(i)] = pickle.Unpickler(fi)
                mcmc_result["mcmc_result{}".format(i)] = pickler["myPickler{}".format(i)].load()

    ## Create variable to store MCMC results
    final_pos = []
    final_PA = []
    final_contr = []
    final_pos_gauss = []
    final_PA_gauss = []
    final_contr_gauss = []

    ## Obtain r, PA, flux and error values
    for i in range(0,len(wl)):
        mcmc = mcmc_result["mcmc_result{}".format(i)]
        chain_40 = mcmc['chain']
        index = np.where(mcmc['AR']>0.4)
        print('Companion number: ', i)

        burnin = 0.3
        chain_40_g = chain_40[index]

        isamples_flat = chain_40_g[:,int(chain_40_g.shape[1]//(1/burnin)):,:].reshape((-1,3))
        mu,sigma = vip_hci.negfc.mcmc_sampling.confidence(isamples_flat,
                                                                        cfd = 68,
                                                                        gaussian_fit = True,
                                                                        verbose=True,
                                                                        save=False,
                                                                        full_output=False,
                                                                        title=source,
                                                                        edgecolor = 'm',facecolor = 'b',range=())
    #
    #     pKey = ['r','theta','f']
    #     PA = val_max[pKey[1]]-90
    #     if PA < -180: PA = val_max[pKey[1]]+360
    #
    #     final_pos.append([val_max[pKey[0]],conf[pKey[0]][0],conf[pKey[0]][1]])
    #     final_PA.append([PA,conf[pKey[1]][0],conf[pKey[1]][1]])
    #     final_contr.append([val_max[pKey[2]],conf[pKey[2]][0],conf[pKey[2]][1]])
    #     final_pos_gauss.append([mu[0],sigma[0]])
    #     final_PA_gauss.append([mu[1],sigma[1]])
    #     final_contr_gauss.append([mu[2],sigma[2]])
    #
    # # Get values into usuable arrays
    # spectra_gauss = np.array([item[0] for item in final_contr_gauss])
    # spectra_gauss_err = np.array([item[1] for item in final_contr_gauss])

# Read FASTWIND calibrated spectra
if fastwind == True:
    ## Open file
    f = open(fastwind_path + 'FLUXCONT', 'r')

    ## Read and ignore header lines
    header1 = f.readline()

    # Define the wavelength and fnue
    fast_wavel = 0
    fast_flux = 0
    fast_wavel = np.array([])
    fast_flux = np.array([])

    ## Loop over lines and extract variables of interest within the wavelength range of IRDIS
    for line in f:
        line = line.strip()
        columns = line.split()
        if float(columns[1]) > 9400 and float(columns[1]) < 16400:
            fast_wavel = np.append(fast_wavel,float(columns[1]))
            fast_flux = np.append(fast_flux,float(columns[2]))
        if float(columns[1]) < 9300:
            break
    f.close()
    # The flux is measured the opposite way as for IRDIS
    fast_wavel = fast_wavel[::-1]
    fast_flux = fast_flux[::-1]

    ## Adjust the model spectra to the same wavelengths as IRDIS
    model_spectra = np.interp(wl*1e4,fast_wavel,fast_flux)

    ## Define some parameters
    model_flux = np.zeros_like(wl)

    ## Put the flux from frequency to wavelength space
    for i in range(len(wl)):
        model_flux[i] = (c*1e10) / (wl[i]*1e4)**2 * 10**(model_spectra[i])

    ## Scale the flux to a distance of Ro
    flux = model_flux/(dist_fast)**2 * (120*rad_fast)**2 #Use 100 Ro for distance to flux measurement

# Companion spectrum calibration
if calib_spec == True:
    ## If spectrum file is specified, then read and store it
    if sspec_file is not None:
        simplex_guess = np.zeros((39,3)) # Set the simplex variable: r, PA, flux
        f = open(sspec_file,'r')
        for line in f:
            line = line.strip()
            columns = line.split()
            simplex_guess[:][0] = float(columns[0])
            simplex_guess[:][1] = float(columns[1])
            simplex_guess[:][2] = float(columns[2])
        f.close()
    ## Read calibrated spectrum and store it
    if calib_star_spec_path is not None:
        calib_spectrum = np.zeros_like(wl)
        f = open(calib_star_spec_path,'r')
        for line in f:
            line = line.strip()
            columns = line.split()
            calib_spectrum = float(columns)
        f.close()
    ## Calculate the contrast spectrum of the companion
    contr_spectra = np.zeros_like(wl)
    contr_spectra_err = np.zeros_like(wl)
    for i in range(0,len(wl)):
        contr_spectra[i] = simplex_guess[i][2]/psf_final_sum[i] # Spectrum
        contr_spectra_err[i] = (spectra_gauss_err[i]/simplex_guess[i][2])*(simplex_guess[i][2]/psf_final_sum[i]) # Error on contrast spectrum !! Requires prior MCMC reading !!

    ## Calculate the calibrated companion spectrum
    calib_comp_spec = contr_spectra * calib_spectrum
    calib_comp_spec_err = calib_comp_spec * contr_spectra_err/contr_spectra

    ## Save the calibrated spectrum with errors if requested
    if save_calib_spec == True:
        comp_spec = np.zeros((39,2))
        for i in range(0,39):
            comp_spec[i][0] = calib_comp_spec
            comp_spec[i][1] = calib_comp_spec_err
        np.savetxt(cspec_file, comp_spec, delimiter='   ')

# Magnitude contrasts
if mag_contr == True:

    ## Define the Y, J and H bands
    Y = wl[0:14]
    J = wl[14:25]
    H = wl[25:39]

    ## Calculate their respective contrast spectra
    contr_spectra_Y = contr_spectra[0:14]
    contr_spectra_J = contr_spectra[14:25]
    contr_spectra_H = contr_spectra[25:39]
    contr_err_Y = contr_err[0:14]
    contr_err_J = contr_err[14:25]
    contr_err_H = contr_err[25:39]

    ## Define some stored parameters
    dmag_Y = np.zeros_like(contr_spec_Y)
    dmag_J = np.zeros_like(contr_spec_J)
    dmag_H = np.zeros_like(contr_spec_H)
    dmag_Y_err = np.zeros_like(contr_spec_Y)
    dmag_J_err = np.zeros_like(contr_spec_J)
    dmag_H_err = np.zeros_like(contr_spec_H)

    for i in range(len(contr_spec_Y)):
        dmag_Y[i] = -2.5*mh.log10(contr_spec_Y[i])
        dmag_H[i] = -2.5*mh.log10(contr_spec_H[i])
        dmag_Y_err[i] = np.abs((1/mh.log(10)) * -2.5 * (contr_err_Y[i]/contr_spec_Y[i]))
        dmag_H_err[i] = np.abs((1/mh.log(10)) * -2.5 * (contr_err_H[i]/contr_spec_H[i]))

    for i in range(len(contr_spec_J)):
        dmag_J[i] = -2.5*mh.log10(contr_spec_J[i])
        dmag_J_err[i] = np.abs((1/mh.log(10)) * -2.5 * (contr_err_J[i]/contr_spec_J[i]))

    if print_mag_contr == True:
        print("Y contrast magnitude: " + dmag_Y + " +/- " + dmag_Y_err)
        print("J contrast magnitude: " + dmag_J + " +/- " + dmag_J_err)
        print("H contrast magnitude: " + dmag_H + " +/- " + dmag_H_err)

# Absolute magnitudes
if abs_mag == True:

    ## Zero-point flux m=0
    star_flux_Y = 2026*10**(star_mag_Y/-2.5)
    star_flux_J = 1600*10**(star_mag_J/-2.5)
    star_flux_H = 1080*10**(star_mag_H/-2.5)

    ## Companion flux
    FluxJy_Y = contr_spec_Y * star_flux_Y
    FluxJy_J = contr_spec_J * star_flux_J
    FluxJy_H = contr_spec_H * star_flux_H

    ## Companion magnitude
    mag_comp = np.zeros([3])
    mag_comp[0] = -2.5 * mh.log10(FluxJy_Y/1.)#Yband zero point flux)
    mag_comp[1] = -2.5 * mh.log10(FluxJy_J/1.)#Yband zero point flux)
    mag_comp[2] = -2.5 * mh.log10(FluxJy_H/1.)#Yband zero point flux)
    mag_comp_err = np.zeros([3])
    mag_comp_err[0] = -2.5 * mh.log10(FluxJy_Y/1.)#Yband zero point flux)
    mag_comp_err[1] = -2.5 * mh.log10(FluxJy_J/1.)#Yband zero point flux)
    mag_comp_err[2] = -2.5 * mh.log10(FluxJy_H/1.)#Yband zero point flux)

    ## Extinction
    BV_obs = 6.37-6.24
    BV_O = -0.26
    Rv = 3.1
    E_BV = BV_obs - BV_O
    Av = Rv * E_BV
    Ay = Av * 1.# Missing
    Aj = Av * 0.282
    Ah = Av * 0.175

    ## Absolute magnitude
    M_comp = np.zeros_like(mag_comp)
    M_comp[0] = mag_comp[0] - Ay - 5*(mh.log10(star_dist)-1)
    M_comp[1] = mag_comp[1] - Aj - 5*(mh.log10(star_dist)-1)
    M_comp[2] = mag_comp[2] - Ah - 5*(mh.log10(star_dist)-1)
    M_comp_err = np.zeros_like(mag_comp)
    M_comp_err[0] = mag_comp[0] - Ay - 5*(mh.log10(star_dist)-1)
    M_comp_err[1] = mag_comp[1] - Aj - 5*(mh.log10(star_dist)-1)
    M_comp_err[2] = mag_comp[2] - Ah - 5*(mh.log10(star_dist)-1)

    ## Print absolute magnitudes
    if print_abs_mag == True:
        print("Y apparent magnitude: " + mag_comp[0] + " +/- " + mag_comp_err[0])
        print("J apparent magnitude: " + mag_comp[1] + " +/- " + mag_comp_err[1])
        print("H apparent magnitude: " + mag_comp[2] + " +/- " + mag_comp_err[2])
        print("Y absolute magnitude: " + M_comp[0] + " +/- " + M_comp_err[0])
        print("J absolute magnitude: " + M_comp[1] + " +/- " + M_comp_err[1])
        print("H absolute magnitude: " + M_comp[2] + " +/- " + M_comp_err[2])
