############################
# Date: 06/08/2019
# Title: Running script for SHIPS for IFS data
# Description: Use this script to run SHIPS for IFS data. In this script you'll find all the necessary parameters to run SHIPS. ONLY SPHERE-DC DATA FOR NOW. VIP and pyKLIP are used.
# VIP version: 0.9.9 (Rainot edit.)
# pyKLIP version: 1.1 NOT IMPLEMENTED YET
############################

# Set up your parameters

## Define images to analyse
cube_filepath = '/Users/alan/Documents/PhD/Data/SPHERE/IFS/QZCardone/ifs_sortframes_dc-IFS_SCIENCE_REDUCED_SPECTRAL_MASTER_CUBE_SORTED-center_im_sorted.fits'
wavelength_filepath = '/Users/alan/Documents/PhD/Data/SPHERE/IFS/QZCardone/ifs_sortframes_dc-IFS_SCIENCE_LAMBDA_INFO-lam.fits'
angles_filepath = '/Users/alan/Documents/PhD/Data/SPHERE/IFS/QZCardone/ifs_sortframes_dc-IFS_SCIENCE_PARA_ROTATION_CUBE_SORTED-rotnth_sorted.fits'
psf_filepath = '/Users/alan/Documents/PhD/Data/SPHERE/IFS/QZCardone/corrected_psf.fits'

## Photometry
comp_pos = (112,55) # Companion position in pixels (X,Y)
psf_pos = (32, 33) # PSF position in pixels (X,Y)
radial_dist = 97.684904197311155 # Radial distance of companion in pixels
position_angle = 159.37210826761003  # Position angle of companion in degrees
noise_aperture_pos_comp = (92,102) # Position in pixels of the circular annulus aperture for noise measurement in the case of the companion
noise_aperture_pos_psf = (12,22) # Position in pixels of the circular annulus aperture for noise measurement in the case of the PSF

## Computing power
ncores = 4 # Number of cores you are willing to share for the computation

## PCA
ncomp_pca = 1 # Number of principal components for PCA

## Do you want to see the image?
see_cube = False # Original cube
see_collapsed_cube = False # Collapsed cube
see_psf_norm = False # Normalised PSF
see_cube_centre = False # Check if the image is centered correctly

## SNR maps
snr_maps = False # Would you like to make and save an SNR map to disk?
snr_map_file = '/Users/alan/Documents/PhD/Data/SPHERE/IFS/QZCardone/SNRmap_VIP.fits' # Finish the file with .fits

## Contrast curves
contrast_curves = False # True or False !! computationally intensive !!
n_branches = 1 # Number of branches for contrast curves

## Spectrum extraction with Simplex Nelder-Mead optimisation
extract_spec = True # Will start the simplex Nelder-Mead optimisation for spectrum extraction
save_spec = True # Save the spectrum to ascii file
sspec_file = '/Users/alan/Documents/PhD/Data/SPHERE/IFS/QZCardone/VIP_simplex.txt' # Filepath to save the Simplex spectrum

## Spectrum extraction with MCMC
extract_mcmc = False # Will compute the MCMC for all 39 wavelengths !! This takes ~1,5h per wavelength and is very computer intensive !!
source = 'QZCar' # Give name for your source
mcmc_path = '/Users/alan/Documents/PhD/Data/SPHERE/IFS/QZCardone/spectra/' # Directory where MCMC results will be stored

## Reading MCMC results
read_mcmc = False # Do you wish to read the MCMC results?
source = 'QZCar' # Give name for your source
mcmc_path = '/Users/alan/Documents/PhD/Data/SPHERE/IFS/QZCardone/spectra/' # Directory where MCMC results are stored

## Load calibrated FASTWIND models of the central star
fastwind = False # Use FASTWIND model spectra for the star
fastwind_path = '/Users/alan/Nextcloud/PhD/Thesis/SPHERE/spectra/fastwind/qzcarAa1/' # Directory where the FASTWIND flux are
rad_fast = 22.1 # Radius of model star
dist_fast = 100. # Distance to consider for the flux of the calibrated spectrum in Ro

## Compute calibrated spectrum of companion
calib_spec = False # Do you wish to calibrate the spectrum of the companion?
calib_star_spec_path = '/Users/alan/Nextcloud/PhD/Thesis/SPHERE/spectra/fastwind/qzcar_fastwind_spec.txt' # Path to calibrated spectrum of central star
sspec_file = '/Users/alan/Documents/PhD/Data/SPHERE/IFS/QZCardone/VIP_simplex.txt' # Path to spectrum file
# ---------------------------------------------------------------------------

# Running script (DO NOT MODIFY)

# Some definitions

## Load libraries
import __init__
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
pxscale = 0.0074 # IFS pixel scale in arcsec/pixel
PA = position_angle + 90 # Correct for VIP unconventional rotation axis

## Open image files
cube = vip_hci.fits.open_fits(cube_filepath)
wl = vip_hci.fits.open_fits(wavelength_filepath)
angs = vip_hci.fits.open_fits(angles_filepath)
psf = vip_hci.fits.open_fits(psf_filepath)

## Define some Parameters
psf_scaled = np.zeros_like(psf) # The psf will need to be scaled
flevel = np.zeros_like(cube[:,0,0,0]) # Flux level for the companion
flevel = np.array(flevel) # Redefinition - why?

## Check the RAW data cubes
if see_cube == True:
    ds9 = vip_hci.Ds9Window()
    ds9.display(cube[0,0])

## Get FWHM of images & normalised PSF
psf_norm, maxflux, fwhm = vip_hci.metrics.normalize_psf(psf, fwhm='fit', size=int(13), verbose=False,full_output=True) # maxflux is a dummy variable
### Plot it
if see_psf_norm == True:
    plot_frames(psf_norm[0], grid=True, size_factor=4)

## Check if the cube is centred correctly by plotting
if see_cube_centre == True:
    plot_frames(vip_hci.preproc.frame_crop(cube[0,0], 50), grid=True, size_factor=4)


# Stellar photometry of the companion

## Collapse the images for better photometry measurement
cube_derot = vip_hci.preproc.cube_derotate(cube,angs) # Rotate the images to the same north
cube_wl_coll = vip_hci.preproc.cube_collapse(cube_derot,wl_cube=True) # Collapse along the rotation axis - 3D image
cube_coll = vip_hci.preproc.cube_collapse(cube_derot,wl_cube=False) # Collapse along the wavelength axis - 2D image

## Check the collapsed data cubes
if see_collapsed_cube == True:
    ds9 = vip_hci.Ds9Window()
    ds9.display(cube_wl_coll[0],cube_coll) # cube_wl_coll on the left and cube_coll on the right

## Aperture photometry of companions and PSF

### Define photometry
noise_phot = np.zeros_like(wl) #Noise photometry
psf_final_sum = np.zeros_like(wl) #PSF photometry
final_sum = np.zeros_like(wl) #Companion photometry

### Apertures
aper_noise_comp = photutils.CircularAnnulus((145,145),noise_aperture_pos_comp[0],noise_aperture_pos_comp[1])
aper_noise_psf = photutils.CircularAnnulus(psf_pos,noise_aperture_pos_psf[0],noise_aperture_pos_psf[1])

### Aperture photometry
for i in range(0,wl.shape[0]):
    ### Apertures dependent on channel
    aper_comp = photutils.CircularAperture(comp_pos, 1./2*fwhm[i])
    aper_psf = photutils.CircularAperture(psf_pos, 1./2*fwhm[i])
    #Noise
    phot_noise = photutils.aperture_photometry(cube_wl_coll[i], aper_noise_comp)
    noise_phot[i] = np.array(phot_noise['aperture_sum'])
    #PSF
    phot_psf = photutils.aperture_photometry(psf[i], aper_psf)
    phot_psf_noise = photutils.aperture_photometry(psf[i], aper_noise_psf)
    psf_bkg_mean = phot_psf_noise['aperture_sum'] / aper_noise_psf.area()
    psf_bkg_sum = psf_bkg_mean * aper_psf.area()
    psf_final_sum[i] = phot_psf['aperture_sum'] - psf_bkg_sum
    #Companion
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
    cube_negfc = vip_hci.metrics.cube_inject_companions(cube,psf_norm,-angs,flevel=-105,plsc=pxscale,rad_dists=[radial_dist],theta=PA) # Remove companion using NEGFC technique
    print("Companion removed")
    print("Computing contrast curve...")
    contrcurve = vip_hci.metrics.contrast_curve(cube_negfc,-angs,psf,np.average(fwhm),pxscale,psf_final_sum,vip_hci.pca.pca,nbranch=n_branches,
              dpi=300, student=False, debug=True ,plot=True, verbose=True, full_output=True, ncomp=ncomp_pca, scale_list=wl)

elif contrast_curves == False:
    print("No contrast curve")


# Spectrum extraction with NM
if extract_spec == True:

    ## Define some parameters
    comp_xycoord = [[comp_pos[0],comp_pos[1]]] # Companion coords
    f_guess_pl = 200. # Flux first guess
    f_range = np.linspace(0.*f_guess_pl,10 *f_guess_pl,400)
    p_in = np.array([radial_dist,PA]) # Regroup companion positions
    simplex_options = {'xtol':1e-2, 'maxiter':500, 'maxfev':1000} # Set the simplex options
    simplex_guess = np.zeros((39,3)) # Set the simplex variable: r, PA, flux

    ## Start Simplex
    for i in range(0,len(wl)):
        print("Wavelength index: ", i + 1) # 39 wavelengths for IFS
        simplex_guess[i] = vip_hci.negfc.firstguess(cube[i],-angs,psf_scaled[i],ncomp_pca,pxscale,comp_xycoord,simplex_options=simplex_options,f_range=f_range,p_ini=p_in,verbose=False) # This takes some time
        print(simplex_guess[i])
## Save the spectrum
if save_spec == True:
    np.savetxt(sspec_file, simplex_guess, delimiter='   ') # Saves to file

# Spectrum extraction with MCMC
if extract_mcmc == True:
    instru= 'IFS36059' # Define instrument parameters
    ann_width=annulus_width # Annulus width of MCMC
    aperture_radius=aperture_width # Aperture radius
    fig_merit='sum' # Summation figure of merit
    outpath = mcmc_path.format(source) # Path to save MCMC files

    print("########## MCMC Sampling starting... ##########")

    nwalkers, itermin, itermax = (100,200,500) # as recommended by Oli W
    for i in range(len(final_sum)): # For each wavelength channel
        initialState = simplex_guess[i] # Take r, PA and flux from simplex
        bounds=[[0.75*initialState[0],1.25*initialState[0]],[0.75*initialState[1],1.25*initialState[1]],[0.75*initialState[2],1.30*initialState[2]]] # Initiate bounds
        output_file = source+'_IFS_wavelength_{}'.format(i) # Save to output file

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
    fast_wavel_a1 = 0
    fast_flux_a1 = 0
    fast_wavel_a1 = np.array([])
    fast_flux_a1 = np.array([])

    ## Loop over lines and extract variables of interest within the wavelength range of IFS
    for line in f:
        line = line.strip()
        columns = line.split()
        if float(columns[1]) > 9400 and float(columns[1]) < 16400:
            fast_wavel_a1 = np.append(fast_wavel_a1,float(columns[1]))
            fast_flux_a1 = np.append(fast_flux_a1,float(columns[2]))
        if float(columns[1]) < 9300:
            break
    f.close()
    # The flux is measured the opposite way as for IFS
    fast_wavel_a1 = fast_wavel_a1[::-1]
    fast_flux_a1 = fast_flux_a1[::-1]

    ## Adjust the model spectra to the same wavelengths as IFS
    model_spectra_a1 = np.interp(wl*1e4,fast_wavel_a1,fast_flux_a1)

    ## Define some parameters
    model_flux_a1 = np.zeros_like(wl)
    model_FLUX_a1 = np.array([])
    model_WL_a1 = np.array([])

    ## Put the flux from frequency to wavelength space
    for i in range(len(wl)):
        model_flux_a1[i] = (c*1e10) / (wl[i]*1e4)**2 * 10**(model_spectra_a1[i])

    ## Scale the flux to a distance of Ro
    flux_a1 = model_flux_a1/(dist_fast)**2 * (120*rad_fast)**2 #Use 100 Ro for distance to flux measurement

    f = open("/Users/alan/Nextcloud/PhD/Thesis/SPHERE/spectra/qzcar_simplex_flux_final.txt",'r')

# Companion spectrum calibration
if calib_spec == True:
    ## Define parameters
    fcomp = np.array([])
    for line in f:
        line = line.strip()
        columns = line.split()
        fcomp = np.append(fcomp,float(columns[1]))
