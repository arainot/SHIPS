############################
# Date: 07/08/2019
# Title: Running script for SHIPS for IFS data
# Description: Use this script to run SHIPS for IFS data. In this script you'll find all the necessary parameters to run SHIPS. ONLY SPHERE-DC DATA FOR NOW. VIP is used.
# VIP version: 0.9.11 (Rainot edit.)
# Python version: 3 ONLY
############################

# Set up your parameters

## Define images to analyse
cube_filepath = '/Users/alan/Documents/PhD/Data/SPHERE/IFS/QZCardone/ifs_sortframes_dc-IFS_SCIENCE_REDUCED_SPECTRAL_MASTER_CUBE_SORTED-center_im_sorted.fits'
wavelength_filepath = '/Users/alan/Documents/PhD/Data/SPHERE/IFS/QZCardone/ifs_sortframes_dc-IFS_SCIENCE_LAMBDA_INFO-lam.fits'
angles_filepath = '/Users/alan/Documents/PhD/Data/SPHERE/IFS/QZCardone/ifs_sortframes_dc-IFS_SCIENCE_PARA_ROTATION_CUBE_SORTED-rotnth_sorted.fits'
psf_filepath = '/Users/alan/Documents/PhD/Data/SPHERE/IFS/QZCardone/psf_corrected_final.fits'
# wavelength_filepath = '/Users/alan/Documents/PhD/Data/SPHERE/IFS/HD93403/ifs_sortframes_dc-IFS_SCIENCE_LAMBDA_INFO-lam.fits'
# cube_filepath = '/Users/alan/Documents/PhD/Data/SPHERE/IFS/HD93403/ifs_sortframes_dc-IFS_SCIENCE_REDUCED_SPECTRAL_MASTER_CUBE_SORTED-center_im_sorted.fits'
# angles_filepath = '/Users/alan/Documents/PhD/Data/SPHERE/IFS/HD93403/ifs_sortframes_dc-IFS_SCIENCE_PARA_ROTATION_CUBE_SORTED-rotnth_sorted.fits'
# psf_filepath = '/Users/alan/Documents/PhD/Data/SPHERE/IFS/HD93403/ifs_sortframes_dc-IFS_SCIENCE_PSF_MASTER_CUBE-median_unsat.fits'

## Photometry

# HD93403
# comp_pos = (128.,168.) # Companion position in pixels from the center of the frame (X,Y)
# psf_pos = (32, 32) # PSF position in pixels (X,Y)
# frame_cent = (145,145) # Center of the frame
# radial_dist = 28.6 # Radial distance of companion in pixels
# position_angle = 159.  # Position angle of companion in degrees
# noise_aperture_pos_comp = (23,33) # Position in pixels of the circular annulus aperture for noise measurement in the case of the companion
# noise_aperture_pos_psf = (12,22) # Position in pixels of the circular annulus aperture for noise measurement in the case of the PSF
# size_psf = 31 # What size PSF would you like to use? ODD VALUE ONLY!!

comp_pos = (110.,54.) # Companion position in pixels from the center of the frame (X,Y)
psf_pos = (32, 32) # PSF position in pixels (X,Y)
frame_cent = (145,145) # Center of the frame
radial_dist = 98 # Radial distance of companion in pixels
position_angle = 159.  # Position angle of companion in degrees
noise_aperture_pos_comp = (92,104) # Position in pixels of the circular annulus aperture for noise measurement in the case of the companion
noise_aperture_pos_psf = (12,22) # Position in pixels of the circular annulus aperture for noise measurement in the case of the PSF
size_psf = 31 # What size PSF would you like to use? ODD VALUE ONLY!!

## Computing power
ncores = 4 # Number of cores you are willing to share for the computation

## Do you want to see the image?
see_cube = False # Original cube
see_collapsed_cube = False # Collapsed cube
see_psf_norm = False # Normalised PSF
see_cube_centre = False # Check if the image is centered correctly

## PCA
ncomp_pca = 1 # Number of principal components for PCA
opti_pca = False # Optimise the number of PCA components?
source_pca = (28.,159.) # Source where to optimise the PCA

## SNR maps
snr_maps = False # Would you like to make and save an SNR map to disk?
snr_map_file = '/home/alan/data/Backup_macbook/SPHERE/IFS/QZCardone/SNRmap_VIP11.fits' # Finish the file with .fits

## Detection
adi_frame = False # Would you like to apply ADI on the frame?
adi_plot = False # Would you like to see the resulting plot?
adi_min_scale = -1 # Minimum colour scale for the ADI plot
adi_max_scale = 1 # Maximum colour scale for the ADI plot
detection = False # Would you like the algorithm to detect sources for you? !! WARNING: this is a simple detection !!

detect_sigma = 5 # What sigma limit would you like for the detection?

## Contrast curves
contrast_curves = False # True or False !! computationally intensive !!
n_branches = 1 # Number of branches for contrast curves

## Photometric errors of PSF
psf_errors = True # Compute the photometric errors of the central star's PSF
psf_errors_save = True # Save the errors to a file?
psf_errors_file = "/Users/alan/Documents/PhD/Data/SPHERE/IFS/QZCardone/PSF_errors.txt"

## Aperture Photometry
plot_aper = False # Plot the aperture photometry of the detected companion?

## Spectrum extraction with Simplex Nelder-Mead optimisation
extract_spec = False # Will start the simplex Nelder-Mead optimisation for spectrum extraction
ann_width = 3 # Annulus width of Simplex
aper_radius = 3 # Aperture Radius of PCA
save_spec = False # Save the spectrum to ascii file
sspec_file = '/Users/alan/Documents/PhD/Data/SPHERE/IFS/HD93403/VIP_simplex.txt' # Filepath to save the Simplex spectrum
plot_sspec = False # Plot the resulting spectrum?

## Spectrum extraction with MCMC
extract_mcmc = False # Will compute the MCMC for all 39 wavelengths !! This takes ~1,5h per wavelength and is very computer intensive !!
source = 'HD93403' # Give name for your source
mcmc_path = '/home/alan/data/Backup_macbook/SPHERE/IFS/HD93403/mcmc/' # Directory where MCMC results will be stored
plot_mcmc = False # Plot the mcmc errors with simplex?

## Reading MCMC results
read_mcmc = False # Do you wish to read the MCMC results?
source = 'QZCar' # Give name for your source
mcmc_path = '/Users/alan/Documents/PhD/Data/SPHERE/IFS/QZCardone/spectra/' # Directory where MCMC results are stored

## Load calibrated FASTWIND models of the central star
fastwind = False # Use FASTWIND model spectra for the star
fastwind_path = '/home/alan/data/Backup_macbook/SPHERE/IFS/HD93403/fastwind_model.txt' # Directory where the FASTWIND flux are
rad_fast = 22.1 # Radius of model star
dist_fast = 100. # Distance to consider for the flux of the calibrated spectrum in Ro

## Compute calibrated spectrum of companion
calib_spec = False # Do you wish to calibrate the spectrum of the companion?
save_calib_spec = False # Would you like to save the calibrated spectrum & associated error?
calib_star_spec_path = '/home/alan/data/Backup_macbook/SPHERE/IFS/HD93403/fastwind_model.txt' # Path to calibrated spectrum of central star
# sspec_file = '/Users/alan/Documents/PhD/Data/SPHERE/IFS/HD93403/VIP_simplex.txt' # Path to spectrum file
cspec_file = '/home/alan/data/Backup_macbook/SPHERE/IFS/HD93403/VIP_calib_spectra.txt' # Path to calibrated spectrum

## Magnitude contrasts
mag_contr = False # Do you to calculate the magnitude contrasts for your companion?
print_mag_contr = False # Do you wish to print the magnitude contrasts to the screen?

## Absolute magnitude !! Work in Progress !!
abs_mag = False # Would you like to calculate the absolute magnitudes of your companion?
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
import math as mh
from scipy.integrate import quad, dblquad
### Error libraries
from uncertainties import ufloat
from uncertainties.umath import *
from uncertainties import unumpy

## Define constants
c = 299792458. # Speed of light
Ro = 6.957e8 # Solar Radius
sr2pc = 44334448.0068964 # Convert steraradians to parsec
pxscale = 0.0074 # IFS pixel scale in arcsec/pixel
PA = position_angle + 90 # Correct for VIP unconventional rotation axis

## Open image files
cube = vip_hci.fits.open_fits(cube_filepath)
wl = vip_hci.fits.open_fits(wavelength_filepath)
angs = vip_hci.fits.open_fits(angles_filepath)
psf = vip_hci.fits.open_fits(psf_filepath)

## Define some Parameters
#psf = np.median(psf, axis=1) # Take the median of all PSFs
psf_scaled = np.zeros_like(psf) # The psf will need to be scaled
flevel = np.zeros_like(cube[:,0,0,0]) # Flux level for the companion
flevel = np.array(flevel) # Redefinition - why?

## Check the RAW data cubes
if see_cube == True:
    ds9 = vip_hci.Ds9Window()
    ds9.display(cube[0,0])

## Get FWHM of images & normalised PSF
psf_med = vip_hci.preproc.cosmetics.cube_crop_frames(psf, size_psf, xy=psf_pos, verbose=True, force=True) # Resize the PSF
psf_norm, maxflux, fwhm = vip_hci.metrics.normalize_psf(psf_med, fwhm='fit',size=None, threshold=None,mask_core=None, model='gauss',imlib='opencv',interpolation='lanczos4',force_odd=True,full_output=True,verbose=False) # maxflux is a dummy variable

### Plot it
if see_psf_norm == True:
    plot_frames(psf_norm[0], grid=True, size_factor=10)

## Check if the cube is centred correctly by plotting
if see_cube_centre == True:
    plot_frames(vip_hci.preproc.frame_crop(cube[0,0], 50), grid=True, size_factor=10)

## Optimise the number of PCA components
if opti_pca == True:
    vip_hci.pca.pca(cube[0], angs, fwhm=fwhm[0], source_xy=(129,169),mask_center_px=None, ncomp=(1, 41, 2))
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
        detect = vip_hci.metrics.detection(fr_adi, fwhm=fwhm[0], psf=psf_norm[0], debug=False, mode='log', snr_thresh=detect_sigma,bkg_sigma=detect_sigma,matched_filter=True,vmin=adi_min_scale,vmax=adi_max_scale,verbose=False,plot=True) # Sigma limit provided by user
        print("Detected sources : " , "\n", detect)
        detect_pos = np.array(detect) # Converted to array in order to be used later
        sys.exit("Sources detected. To continue, please input the target coordinates in the script and skip this process.")

# Stellar photometry of the companion

## Collapse the images for better photometry measurement
cube_derot = vip_hci.preproc.cube_derotate(cube,angs) # Rotate the images to the same north
cube_wl_coll = vip_hci.preproc.cube_collapse(cube_derot,wl_cube=True) # Collapse along the rotation axis - 3D image
# cube_coll = vip_hci.preproc.cube_collapse(cube_derot,wl_cube=False) # Collapse along the wavelength axis - 2D image
# cube_wl_coll = np.zeros_like(cube[:,0,:,:])
# for i in range(len(wl)):
#         cube_wl_coll[i] = vip_hci.hci_postproc.median_sub(cube[i],-angs,fwhm=fwhm[i],verbose=False) # Rotate & collapse along the rotation axis - 3D image
# cube_derot = vip_hci.preproc.cube_derotate(cube,angs) # Rotate the images to the same north

## Check the collapsed data cubes
if see_collapsed_cube == True:
    ds9 = vip_hci.Ds9Window()
    ds9.display(cube_wl_coll[0]) # cube_wl_coll on the left and cube_coll on the right

## Aperture photometry of companions and PSF

### Define photometry
noise_phot = np.zeros_like(wl) #Noise photometry
psf_final_sum = np.zeros_like(wl) #PSF photometry
final_sum = np.zeros_like(wl) #Companion photometry

### Apertures
aper_noise_comp = photutils.CircularAnnulus(frame_cent,noise_aperture_pos_comp[0],noise_aperture_pos_comp[1])
aper_noise_psf = photutils.CircularAnnulus(psf_pos,noise_aperture_pos_psf[0],noise_aperture_pos_psf[1])

### Aperture photometry
for i in range(0,wl.shape[0]):
    ### Apertures dependent on channel
    aper_comp = photutils.CircularAperture(comp_pos, 1./2*fwhm[i])
    aper_psf = photutils.CircularAperture(psf_pos, 1./2*fwhm[i])
    ### Noise
    phot_noise = photutils.aperture_photometry(cube_wl_coll[i], aper_noise_comp)
    noise_phot[i] = np.array(phot_noise['aperture_sum'])
    ### PSF
    phot_psf = photutils.aperture_photometry(psf[i], aper_psf)
    phot_psf_noise = photutils.aperture_photometry(psf[i], aper_noise_psf)
    psf_bkg_mean = phot_psf_noise['aperture_sum'] / aper_noise_psf.area()
    psf_bkg_sum = psf_bkg_mean * aper_psf.area()
    # psf_final_sum[i] = phot_psf['aperture_sum'] - psf_bkg_sum
    psf_final_sum[i] = maxflux[i] - psf_bkg_sum
    ### Companion
    phot = photutils.aperture_photometry(cube_wl_coll[i], aper_comp)
    bkg_mean = (phot_noise['aperture_sum']-phot['aperture_sum']) / (aper_noise_comp.area()-aper_comp.area())
    bkg_sum = bkg_mean * aper_comp.area()
    final_sum[i] = phot['aperture_sum'] - bkg_sum

### Scaling the PSF for normalisation -- SHOULD I JUST TAKE PSF_NORM INSTEAD?
psf_scaled = np.zeros_like(psf)
for i in range (0,len(psf)):
    psf_scaled[i] = psf[i]/psf_final_sum[i]

## SNR maps
if snr_maps == True:
    snrmap = vip_hci.metrics.snrmap(vip_hci.pca.pca(cube, -angs, scale_list=wl, ncomp=ncomp_pca, verbose=True), fwhm[0], nproc=ncores, plot=True)
    vip_hci.fits.write_fits(snr_map_file,snrmap) # Write SNR maps to file

## Contrast curve
if contrast_curves == True:
    cube_negfc = vip_hci.metrics.cube_inject_companions(cube,psf_norm,-angs,flevel=-final_sum,plsc=pxscale,rad_dists=[radial_dist],theta=PA) # Remove companion using NEGFC technique
    print("Companion removed")
    print("Computing contrast curve...")
    contrcurve = vip_hci.metrics.contrast_curve(cube_negfc,-angs,psf_norm,np.average(fwhm),pxscale,psf_final_sum,vip_hci.pca.pca,nbranch=n_branches,
              dpi=300, student=False, debug=True ,plot=True, verbose=True, full_output=True, ncomp=ncomp_pca, scale_list=wl)

## PSF error calculation
if psf_errors == True:
    psferr = vip_hci.fits.open_fits(psf_filepath) # Open the raw PSFs
    stddev_psf = np.zeros_like(wl) # Create an array for the stored standard deviation
    for i in range(len(wl)): # Loop over the wavelengths
        psferr_med = vip_hci.preproc.cosmetics.cube_crop_frames(psferr[i], size_psf, xy=psf_pos, verbose=True, force=True) # Resize the PSF
        psf_norm_err, maxflux_err, fwhm_err = vip_hci.metrics.normalize_psf(psferr_med, fwhm='fit',size=None, threshold=None,mask_core=None, model='gauss',imlib='opencv',interpolation='lanczos4',force_odd=True,full_output=True,verbose=False) # Measure the maximum flux for each PSF
        stddev_psf[i] = np.std(maxflux_err,ddof=1) # Calculate the standard deviation for the PSFs
    if psf_errors_save: # Save the error
        np.savetxt(psf_errors_file,stddev_psf,delimiter='   ') # Saves to file

# Spectrum extraction with NM
if extract_spec == True:

    ## Define some parameters
    comp_xycoord = [(comp_pos[0],comp_pos[1])] # Companion coords
    f_guess_pl = max(final_sum) # Flux first guess as the maximum value of the flux
    f_range = np.linspace(0.*f_guess_pl,10*f_guess_pl,100)
    p_in = np.array([radial_dist,PA]) # Regroup companion positions
    simplex_options = {'xtol':1e-2, 'maxiter':500, 'maxfev':1000} # Set the simplex options
    simplex_guess = np.zeros((len(wl),3)) # Set the simplex variable: r, PA, flux

    ## Start Simplex
    for i in range(0,len(wl)):
        print("Wavelength index: ", i + 1) # 39 wavelengths for IFS
        simplex_guess[i] = vip_hci.negfc.firstguess(cube[i],-angs,psf_norm[i],ncomp=ncomp_pca,plsc=pxscale,planets_xy_coord=comp_xycoord,fwhm=fwhm[i],annulus_width=ann_width,aperture_radius=aper_radius,simplex_options=simplex_options,f_range=f_range,simplex=True,fmerit='sum',collapse='median',svd_mode='lapack',scaling=None,verbose=False,plot=False,save=False)
        #simplex_guess[i] = vip_hci.negfc.simplex_optim.firstguess(cube[i], -angs, psf_norm[i], ncomp=1, plsc=0.0074,fwhm=fwhm[i], annulus_width=3, aperture_radius=2, planets_xy_coord=comp_xycoord, cube_ref=None,svd_mode='lapack',f_range=f_range, simplex=True,fmerit='sum',scaling=None, simplex_options=simplex_options,collapse='median',verbose=False)

    ## Save the spectrum
    if save_spec == True:
        np.savetxt(sspec_file, simplex_guess, delimiter='   ') # Saves to file
        print("Spectrum saved successfully!")

# Spectrum extraction with MCMC
if extract_mcmc == True:

    ## If the spectrum extraction was initiated previously, we do not need to read the resulting file. In case it was not, we read the file specified
    if extract_spec == False:
        ## Read the file containing the simplex minimization values and save the values to an array
        simplex_guess = np.zeros((len(wl),3)) # Simplex array
        j = 0
        with open(sspec_file) as f:
            for line in f:
                line = line.strip()
                columns = line.split()
                simplex_guess[j,0] = float(columns[0]) # Radius
                simplex_guess[j,1] = float(columns[1]) # PA
                simplex_guess[j,2] = float(columns[2]) # Flux
                j+=1
        f.close

    instru= 'IFS36059' # Define instrument parameters
    fig_merit='sum' # Summation figure of merit
    outpath = mcmc_path.format(source) # Path to save MCMC files

    print("########## MCMC Sampling starting... ##########")

    nwalkers, itermin, itermax = (100,200,500) # as recommended by Oli W
    for i in range(len(final_sum)): # For each wavelength channel
        initialState = simplex_guess[i] # Take r, PA and flux from simplex
        bounds=[[0.75*initialState[0],1.25*initialState[0]],[0.75*initialState[1],1.25*initialState[1]],[0.75*initialState[2],1.30*initialState[2]]] # Initiate bounds
        output_file = source+'_IFS_wavelength_{}'.format(i) # Save to output file

        chain_40 = vip_hci.negfc.mcmc_negfc_sampling(cube[i], -angs,  psf_norm[i], ncomp_pca, pxscale, initialState, ann_width,
                                                 aper_radius, cube_ref=None, svd_mode='lapack', nwalkers=nwalkers,
                                                 bounds=bounds, niteration_min=itermin,
                                                 niteration_limit=itermax, check_maxgap=50, nproc= ncores,
                                                 output_file=output_file, display=False,verbosity=1, save=True,
                                                 rhat_threshold=1.01, niteration_supp=0, fmerit=fig_merit) # MCMC run per channel
    print("########## MCMC Sampling done! ##########")

## Read MCMC files
if read_mcmc == True:
    import pickle # Import important MCMC libraries
    from pickle import Pickler
    pickler={}
    mcmc_result={}
    outpath = mcmc_path.format(source) # Path to save MCMC files
    for i in range(0,len(wl)): # Read all channels and store them to variables
        with open(outpath+source+'_IFS_wavelength_{}/MCMC_results'.format(i),'rb') as fi:
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
        print('Wavelength channel: ', i)

        burnin = 0.8
        chain_40_g = chain_40[index]

        isamples_flat = chain_40_g[:,int(chain_40_g.shape[1]//(1/burnin)):,:].reshape((-1,3))
        val_max,conf,mu,sigma = vip.negfc.mcmc_sampling.confidence(isamples_flat,
                                                                        cfd = 68,
                                                                        gaussianFit = True,
                                                                        verbose=False,
                                                                        save=False,
                                                                        full_output=True,
                                                                        title=source,
                                                                        edgecolor = 'm',
                                                                        facecolor = 'b',range=())

        pKey = ['r','theta','f']
        PA = val_max[pKey[1]]-90
        if PA < -180: PA = val_max[pKey[1]]+360

        final_pos.append([val_max[pKey[0]],conf[pKey[0]][0],conf[pKey[0]][1]])
        final_PA.append([PA,conf[pKey[1]][0],conf[pKey[1]][1]])
        final_contr.append([val_max[pKey[2]],conf[pKey[2]][0],conf[pKey[2]][1]])
        final_pos_gauss.append([mu[0],sigma[0]])
        final_PA_gauss.append([mu[1],sigma[1]])
        final_contr_gauss.append([mu[2],sigma[2]])

    # Get values into usuable arrays
    spectra_gauss = np.array([item[0] for item in final_contr_gauss])
    spectra_gauss_err = np.array([item[1] for item in final_contr_gauss])

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

    ## Loop over lines and extract variables of interest within the wavelength range of IFS
    for line in f:
        line = line.strip()
        columns = line.split()
        if float(columns[1]) > 9400 and float(columns[1]) < 16400:
            fast_wavel = np.append(fast_wavel,float(columns[1]))
            fast_flux = np.append(fast_flux,float(columns[2]))
        if float(columns[1]) < 9300:
            break
    f.close()
    # The flux is measured the opposite way as for IFS
    fast_wavel = fast_wavel[::-1]
    fast_flux = fast_flux[::-1]

    ## Adjust the model spectra to the same wavelengths as IFS
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

    ## If the spectrum extraction was initiated previously, we do not need to read the resulting file. In case it was not, we read the file specified
    if extract_spec == False:
        ## Read the file containing the simplex minimization values and save the values to an array
        simplex_flux = np.zeros_like(wl) # Simplex array
        j = 0
        with open(sspec_file) as f:
            for line in f:
                line = line.strip()
                columns = line.split()
                simplex_flux[j] = float(columns[2])
                j+=1
        f.close
    ## Otherwise, carry on
    else:
        for i in range(len(wl)):
            simplex_flux[i] = simplex_guess[i][2]

    ## Compute the contrast spectrum
    contr_spectra = simplex_flux / psf_final_sum
    #contr_spectra_err

    ## Define the Y, J and H bands
    Y = wl[0:14]
    J = wl[14:25]
    H = wl[25:39]

    ## Calculate their respective contrast spectra
    contr_spectra_Y = contr_spectra[0:14]
    contr_spectra_J = contr_spectra[14:25]
    contr_spectra_H = contr_spectra[25:39]
    # contr_err_Y = contr_err[0:14]
    # contr_err_J = contr_err[14:25]
    # contr_err_H = contr_err[25:39]

    ## Define some stored parameters
    dmag_Y = np.zeros_like(contr_spectra_Y)
    dmag_J = np.zeros_like(contr_spectra_J)
    dmag_H = np.zeros_like(contr_spectra_H)
    # dmag_Y_err = np.zeros_like(contr_spectra_Y)
    # dmag_J_err = np.zeros_like(contr_spectra_J)
    # dmag_H_err = np.zeros_like(contr_spectra_H)

    for i in range(len(contr_spectra_Y)):
        dmag_Y[i] = -2.5*mh.log10(contr_spectra_Y[i])
        dmag_H[i] = -2.5*mh.log10(contr_spectra_H[i])
        # dmag_Y_err[i] = np.abs((1/mh.log(10)) * -2.5 * (contr_err_Y[i]/contr_spectra_Y[i]))
        # dmag_H_err[i] = np.abs((1/mh.log(10)) * -2.5 * (contr_err_H[i]/contr_spectra_H[i]))

    for i in range(len(contr_spectra_J)):
        dmag_J[i] = -2.5*mh.log10(contr_spectra_J[i])
        # dmag_J_err[i] = np.abs((1/mh.log(10)) * -2.5 * (contr_err_J[i]/contr_spectra_J[i]))

    if print_mag_contr == True:
        print("Y contrast magnitude: ", dmag_Y)# + " +/- " + dmag_Y_err)
        print("J contrast magnitude: ", dmag_J)# + " +/- " + dmag_J_err)
        print("H contrast magnitude: ", dmag_H)# + " +/- " + dmag_H_err)

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

# Plotting
## Aperture Photometry
if plot_aper == True:
    plt.figure(figsize=(12, 9))
    #plt.title('Aperture Photometry IFS')
    #plt.legend()
    # ax.get_xaxis().tick_bottom()
    # ax.get_yaxis().tick_left()
    plt.ylim(0, 1.1*max(final_sum))
    plt.grid(True, 'major', 'y', ls='--', lw=.5, c='k', alpha=.5)
    #plt.grid(True, 'major', 'x', ls='--', lw=.5, c='k', alpha=.3)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.ylabel("Flux [ADU/s]", fontsize=20)
    plt.xlabel('Wavelength [$\mathring{A}$]', fontsize=20)
    plt.plot(wl*1e4, final_sum,lw=2.8)
    plt.show()

## Simplex Optim
if plot_sspec == True:
    simplex_flux = np.zeros_like(wl)
    for i in range(len(wl)):
        simplex_flux[i] = simplex_guess[i][2]
    plt.figure(figsize=(12, 9))
    #plt.title('Aperture Photometry IFS')
    #plt.legend()
    # ax.get_xaxis().tick_bottom()
    # ax.get_yaxis().tick_left()
    plt.ylim(0, 1.1*max(simplex_flux))
    plt.grid(True, 'major', 'y', ls='--', lw=.5, c='k', alpha=.5)
    #plt.grid(True, 'major', 'x', ls='--', lw=.5, c='k', alpha=.3)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.ylabel("Simplex flux [ADU/s]", fontsize=20)
    plt.xlabel('Wavelength [$\mathring{A}$]', fontsize=20)
    plt.plot(wl*1e4, simplex_flux,lw=2.8)
    plt.show()
