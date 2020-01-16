psf_norm, maxflux, fwhm = vip_hci.metrics.normalize_psf(psf, fwhm='fit', size=int(12), verbose=False,full_output=True)


x = vip_hci.metrics.cube_inject_companions(cube,psf_norm,-angs,flevel=200,plsc=pxscale,rad_dists=100,theta=100)
cube_derot = vip_hci.preproc.cube_derotate(x,angs)
ds9.display(cube_derot)

f_guess_pl = 200. # Flux first guess
f_range = np.linspace(0.*f_guess_pl,5 *f_guess_pl,25)

for i in range(0,len(wl)):
   print("Wavelength index: ", i + 1) # 39 wavelengths for IFS
   simplex_guess[i] = vip_hci.negfc.firstguess(x[i],-angs,psf_norm[i],ncomp=ncomp_pca,plsc=pxscale ,planets_xy_coord=[(128,244)],fwhm=fwhm[i],annulus_width=3,aperture_radius=3,simplex_options=simplex_options,f_range=f_range,verbose=True,simplex=True,plot=False)
   print(simplex_guess[i])

firstguess(cube, angs, psfn, ncomp, plsc, planets_xy_coord, fwhm=4,
               annulus_width=3, aperture_radius=4, cube_ref=None,
               svd_mode='lapack', scaling=None, fmerit='sum', imlib='opencv',
               interpolation='lanczos4', collapse='median', p_ini=None,
               f_range=None, simplex=True, simplex_options=None, plot=False,
               verbose=True, save=False)

nc = np.linspace(1,21,num=11)

vip_hci.pca.pca(cube, -angs,scale_list=wl, fwhm=fwhm[0],mask_center_px=None,adimsdi='double',ncomp=(1, 41, 2),nproc=4)

vip_hci.pca.pca(cube[0], angs, fwhm=fwhm[0], source_xy=(501,526),mask_center_px=None, ncomp=(1, 21, 2))

frame = vip_hci.pca.pca(cube[0],-angs,scale_list=wl,adimsdi='double',ncomp=1)
frame = vip_hci.pca.pca(cube_wl_coll,scale_list=wl,ncomp=1)

vip.pca.pca(psf_med_orig,scale_list=wl,ncomp=5)

plot_frames(vip_hci.pca.pca(cube[0],-angs,fwhm=fwhm[0],ncomp=1), grid=False, size_factor=10,vmin=-1,vmax=1)

plot_frames(vip_hci.pca.pca(x,-angs,fwhm=fwhm,scale_list=wl,ncomp=1), grid=False, size_factor=10,vmin=-1,vmax=1)
plot_frames(f,grid=False, size_factor=10,vmin=-1,vmax=1)

for i in range(0,10):
    plot_frames(f[i],grid=False, size_factor=10,vmin=-1,vmax=1)
    plot_frames(f[1],grid=False, size_factor=10,vmin=-1,vmax=1)

plot_cubes(psf_norm, grid=True, size_factor=10)


mh.sqrt(10**(5.41)*3.828*10**26/(4*mh.pi*(40100)**4*(5.67*10**-8)))/(6.957*10**8)

g = 6.674*10**(-11) * 37.3*2e30/((10.5*6.955e8)**2)
4.3e-3*3.1e16*10**6 * 68.5/((22.1*6.955e8)**2)

    contrcurve = vip_hci.metrics.contrast_curve(cube,-angs,psf,np.average(fwhm),pxscale,psf_final_sum,vip_hci.pca.pca,nbranch=n_branches,
              dpi=300, student=True, debug=True ,plot=True, verbose=True, full_output=True, ncomp=ncomp_pca, scale_list=wl)
psf_norm, maxflux, fwhm = vip_hci.metrics.normalize_psf(psf, fwhm='fit', size=None,verbose=False,full_output=True)

import os
rootdir = '/Users/alan/Documents/PhD/Data/SPHERE/IFS'
for subdir, dirs, files in os.walk(rootdir):
    if subdir.startswith("/Users/alan/Documents/PhD/Data/SPHERE/IFS/Q") or subdir.startswith("/Users/alan/Documents/PhD/Data/SPHERE/IFS/HD93403") or subdir.startswith("/Users/alan/Documents/PhD/Data/SPHERE/IFS/P"):
        continue
    if subdir.startswith("/Users/alan/Documents/PhD/Data/SPHERE/IFS/H") or subdir.startswith("/Users/alan/Documents/PhD/Data/SPHERE/IFS/V"):
        for file in files:
            if file.endswith(".fits"):
                print(os.path.join(subdir, file))

import pathlib
path = pathlib.PurePath('/Users/alan/Documents/PhD/Data/SPHERE/IFS/HD93403')
path.name


------- SIMPLEX -------

psf=open_fits('./corrected_psf.fits', n=0, ignore_missing_end=True)
size=31
psf_med=vip_hci.preproc.cosmetics.cube_crop_frames(psf, size, xy=(32, 32), verbose=True, force=True)
#cube_cropped=vip_hci.preproc.cosmetics.cube_crop_frames(cube, 289, xy=(145, 145), verbose=True, force=True)
psfn = vip_hci.hci_dataset.normalize_psf(psf_med, fwhm='fit', size=None, threshold=None, mask_core=None, model='gauss',
                                         imlib='opencv', interpolation='lanczos4', force_odd=True, full_output=True,
                                         verbose=False)
#plots(psfn[0])
star_flux=psfn[1]
fwhm_full=psfn[2]
print('flux =  ', star_flux[0])

psf_norm=psfn[0]

# Estimating Flux
# Initial guess of the position by examining a flux frame or SNR map
planet_xycoord = np.array([[109,54]])
# Naive minimization of the chi^2 by trying this grid of values for the flux
f_range = np.linspace(10,600,50)

simplex_options = {'xtol':1e-1, 'maxiter':500, 'maxfev':400}

# plot with the behaviour of the chi^2 with the naive minimization
figure_options = {'color':'b','marker':'o',
                  'xlim': [f_range[0]-100,f_range[-1]+100],
                  'title':r'$\chi^2_{r}$ vs flux'}

spectra_1_600 = np.zeros_like(wl)
for i in range(0,len(wl)):
    print("Wavelength index: ", i + 1)
    fake_comp = vip_hci.negfc.simplex_optim.firstguess(cube[i], -angs, psf_norm[i], ncomp=1, plsc=0.0074,
                                                           fwhm=fwhm_full[i], annulus_width=3, aperture_radius=2,
                                                           planets_xy_coord=planet_xycoord, cube_ref=None,
                                                           svd_mode='lapack',f_range=f_range, simplex=True,
                                                           fmerit='sum',scaling=None, simplex_options=simplex_options,
                                                           collapse='median',verbose=False)
    print(fake_comp[2])
    spectra_1_600[i] = fake_comp[2]

plt.figure(figsize=(12, 9))
#plt.title('Aperture Photometry IFS')
#plt.legend()
# ax.get_xaxis().tick_bottom()
# ax.get_yaxis().tick_left()
plt.ylim(0, 1.1*max(spectra_1_600))
plt.grid(True, 'major', 'y', ls='--', lw=.5, c='k', alpha=.5)
#plt.grid(True, 'major', 'x', ls='--', lw=.5, c='k', alpha=.3)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.ylabel("Simplex flux [ADU/s]", fontsize=20)
plt.xlabel('Wavelength [$\mathring{A}$]', fontsize=20)
plt.plot(wl*1e4, spectra_1_600,lw=2.8)
plt.show()


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

contrcurve = vip_hci.metrics.contrast_curve(cube,-angs,psf_norm,np.average(fwhm),pxscale,psf_final_sum[0],vip_hci.pca.pca,nbranch=n_branches,
       dpi=300, student=False, debug=True ,plot=True, verbose=True, full_output=True, ncomp=ncomp_pca, scale_list=wl)

cubeM = cube.as_matrix().astype(float)

vip_hci.pca.pca(cube, -angs, cube_ref=None, scale_list=wl, ncomp=1,
        svd_mode='lapack', scaling="temp-standard", mask_center_px=None, source_xy=None,
        delta_rot=1, fwhm=fwhm, adimsdi='single', crop_ifs=True, imlib='opencv',
        interpolation='lanczos4', collapse='median', check_memory=True,
        batch=None, nproc=4, full_output=False, verbose=True)

betapic = vip_hci.hci_dataset.Dataset(cube=cube, angles=-angs, psf=psf)
scaler = preprocessing.StandardScaler()
Xc = scaler.fit_transform(cube)

vip_hci.metrics.throughput(cube, -angs, psf, fwhm, pxscale,
                        nbranch=1, theta=0, inner_rad=1,
                        wedge=(0,360), fc_snr=100, full_output=True,
                        algo=vip_hci.pca.pca, imlib='opencv', verbose=True)

f = open('/Users/alan/Documents/PhD/Data/SPHERE/IRDIS/QZCar/VIP_simplex_K2.txt','r')

simplex_guess = np.zeros((19,3)) # Set the simplex variable: r, PA, flux
contr_spectra = np.zeros((19))
#simplex_guess = np.array([]) # Set the simplex variable: r, PA, flux
j=0
for line in f:
    line = line.strip()
    columns = line.split()
    simplex_guess[j][0] = float(columns[0])
    simplex_guess[j][1] = float(columns[1])
    simplex_guess[j][2] = float(columns[2])
    j+=1
#print(simplex_guess)
f.close()

## Define some stored parameters
dmag = np.zeros(19)

for i in range(len(dmag)-2):
    contr_spectra = simplex_guess[i][2]/psf_final_sum[1]
    dmag[i] = -2.5*mh.log10(contr_spectra)
    print("contrast magnitude: ", dmag[i])


fk1 = 2155.994034403564
fk2 = 1482.5920012322508
r=211.1749682-1.3
t=33.85788191-0.02
cube_negfc = np.zeros_like(cube)
for i in range(0,2):
    if i ==0:
        cube_negfc[i] = vip_hci.metrics.cube_inject_companions(cube[i],psf_norm[i],-angs,flevel=-fk1-530,plsc=pxscale,rad_dists=r,theta=t)
    if i ==1:
        cube_negfc[i] = vip_hci.metrics.cube_inject_companions(cube[i],psf_norm[i],-angs,flevel=-fk2,plsc=pxscale,rad_dists=r,theta=t)
plot_frames(cube_negfc[0,0], size_factor=10, vmin=-10, vmax=10)
vip_hci.fits.write_fits("/home/alan/Desktop/cube_negfc.fits",cube_negfc)
c = np.zeros_like(cube)
for i in range(0,2):
    if i ==0:
        c[i] = vip_hci.negfc.cube_planet_free([(r,t,fk1+510)],cube[i],-angs,psf_norm[i],plsc=pxscale)
    if i ==1:
        c[i] = vip_hci.negfc.cube_planet_free([(r,t,fk2)],cube[i],-angs,psf_norm[i],plsc=pxscale)
plot_frames(c[0,0], size_factor=10, vmin=-10, vmax=10)

cube_negfc = vip_hci.metrics.cube_inject_companions(cube[0],psf_norm[0],-angs,flevel=800,plsc=pxscale,rad_dists=250,theta=100)
cube_wl_coll = vip_hci.hci_postproc.median_sub(cube_negfc,-angs,fwhm=fwhm[0],verbose=False)
ds9 = vip_hci.Ds9Window()
ds9.display(cube_negfc[0])
f_range_K1 = np.linspace(0.2*25, 30*25,200)
vip_hci.negfc.firstguess(cube_negfc,-angs,psf_norm[0],ncomp=ncomp_pca,plsc=pxscale,planets_xy_coord=comp_xycoord,fwhm=fwhm[0],annulus_width=ann_width,aperture_radius=aper_radius,simplex_options=simplex_options,f_range=f_range_K1,simplex=True,fmerit='sum',collapse='median',svd_mode='lapack',scaling=None,verbose=True,plot=False,save=False)

c = np.zeros_like(cube)
res = np.zeros_like(cube)
res_der = np.zeros_like(cube)
for i in range(0,2):
    c[i] = vip_hci.pca.utils_pca.pca_annulus(cube[i], -angs, 1, 20, 384)
plot_frames(c[0,0], size_factor=10, vmin=-1, vmax=1)

r = np.array([])
theta = np.array([])
f = np.array([])
flux = np.array([])
c = np.zeros_like(cube)
with open('/Users/alan/Documents/PhD/Data/SPHERE/IFS/QZCardone/VIP_simplex.txt') as f:
    for line in f:
        line = line.strip()
        columns = line.split()
        r = np.append(r,float(columns[0]))
        theta = np.append(theta,float(columns[1]))
        f = np.append(f,float(columns[2]))

for i in range(len(f)):
    flux=np.append(flux,f[i+1])

for i in range(len(r)):
    c[i]=vip_hci.negfc.cube_planet_free([(r[i]-1,theta[i]-0.6,flux[i])],cube[i],-angs,psf_norm[i],plsc=0.0075, imlib='opencv', interpolation='lanczos4')

cube_filepath = '/Users/alan/Documents/PhD/Data/SPHERE/IRDIS/QZCar/cube_free_IRDIS.fits'
cube_mask = vip_hci.fits.open_fits(cube_filepath)
d = vip_hci.hci_postproc.median_sub(cube_mask[0],-angs,fwhm=fwhm[0],verbose=False)
ds9.display(d)

r_K1 = np.array([])
theta_K1 = np.array([])
f_K1 = np.array([])
r_K2 = np.array([])
theta_K2 = np.array([])
f_K2 = np.array([])
with open('/Users/alan/Nextcloud/PhD/Thesis/SPHERE/spectra/summary_companions.txt') as f:
    header1 = f.readline()
    for line in f:
        line = line.strip()
        columns = line.split()
        r_K1 = np.append(r_K1,float(columns[2]))
        theta_K1 = np.append(theta_K1,float(columns[3])+90)
        f_K1 = np.append(f_K1,float(columns[4]))
f.close
with open('/Users/alan/Documents/PhD/Data/SPHERE/IRDIS/QZCar/VIP_simplex_K2.txt') as f:
    header1 = f.readline()
    header2 = f.readline()
    header3 = f.readline()
    for line in f:
        line = line.strip()
        columns = line.split()
        r_K2 = np.append(r_K2,float(columns[0]))
        theta_K2 = np.append(theta_K2,float(columns[1]))
        f_K2 = np.append(f_K2,float(columns[2]))
f.close
for i in range(len(r_K1)-2):
    if i==0:
        c[0]=vip_hci.negfc.cube_planet_free([(r_K1[i],theta_K1[i],f_K1[i])],cube[0],-angs,psf_norm[0],plsc=0.01225, imlib='opencv', interpolation='lanczos4')
    else:
        c[0]=vip_hci.negfc.cube_planet_free([(r_K1[i],theta_K1[i],f_K1[i])],c[0],-angs,psf_norm[0],plsc=0.01225, imlib='opencv', interpolation='lanczos4')
for i in range(len(r_K2)-2):
    if i==0:
        c[1]=vip_hci.negfc.cube_planet_free([(r_K2[i],theta_K2[i],f_K2[i])],cube[1],-angs,psf_norm[1],plsc=0.01225, imlib='opencv', interpolation='lanczos4')
    else:
        c[1]=vip_hci.negfc.cube_planet_free([(r_K2[i],theta_K2[i],f_K2[i])],c[1],-angs,psf_norm[1],plsc=0.01225, imlib='opencv', interpolation='lanczos4')
vip_hci.fits.write_fits("/home/alan/Desktop/test2.fits",c)
vip_hci.hci_postproc.median_sub(c[0],-angs,fwhm=fwhm[0],verbose=False)

c = np.zeros_like(cube)
c[0]=vip_hci.negfc.cube_planet_free([(2.100067681452761690e+02,3.378271377905575434e+01,2.605832969178423355e+03)],cube[0],-angs,psf_norm[0],plsc=0.01225, imlib='opencv', interpolation='lanczos4')
c[1]=vip_hci.negfc.cube_planet_free([(2.099282631149717986e+02,3.381485450865740461e+01,1.651860225089150845e+03)],cube[1],-angs,psf_norm[1],plsc=0.01225, imlib='opencv', interpolation='lanczos4')
vip_hci.fits.write_fits("/Users/alan/Desktop/cube_free_Ad.fits",c)
t=vip_hci.hci_postproc.median_sub(c[1],-angs,fwhm=fwhm[1],verbose=False)
ds9.display(t)

c_E = np.zeros_like(c)
r_K1 = [2.1109824803199433063e+02]
theta_K1 = [3.387762229578338804e+01]
f_K1 = 2.781458963416393090e+03
r_K2 = 2.099286977127382556e+02
theta_K2 = 3.381478951572226777e+01
f_K2 = 1.651526386738957626e+03
c_E[0]=vip_hci.negfc.cube_planet_free([(r_K1,theta_K1,f_K1)],c[0],-angs,psf_norm[0],plsc=0.01225, imlib='opencv', interpolation='lanczos4')
c_E[1]=vip_hci.negfc.cube_planet_free([(r_K2,theta_K2,f_K2)],c[1],-angs,psf_norm[1],plsc=0.01225, imlib='opencv', interpolation='lanczos4')
t=vip_hci.hci_postproc.median_sub(c_E[1],-angs,fwhm=fwhm[1],verbose=False)
ds9.display(t)

contrcurve = vip_hci.metrics.contrast_curve(cube,-angs,psf_norm,np.average(fwhm),pxscale,psf_final_sum,vip_hci.pca.pca,nbranch=n_branches,
          dpi=300, student=False, debug=True ,plot=True, verbose=True, full_output=True, ncomp=ncomp_pca, scale_list=wl)
c = vip_hci.hci_postproc.median_sub(cube[:],-angs,fwhm=fwhm[:])

mask = vip_hci.metrics.mask_source_centers(cube[0],fwhm[0],440,476)

vip_hci.metrics.frame_inject_companion
frame = array.copy()
frame = mask_circle(frame, radius=2*fwhm)


vip_hci.metrics.contrcurve.contrast_curve(cube, -angs, psf_norm, fwhm=np.average(fwhm),
                                          pxscale=0.01225, starphot=psf_final_sum, algo=vip_hci.pca.pca, sigma=5, nbranch=1,
                                          theta=0, inner_rad=1, wedge=(0, 360), fc_snr=100,
                                          student=True, transmission=None, smooth=True,
                                          interp_order=2, plot=True, dpi=100, imlib='opencv',
                                          debug=True, verbose=True, full_output=False, save_plot=None,
                                          object_name=None, frame_size=None, figsize=(8, 4), ncomp=1,adimsdi='double',
                                          scale_list=wl/wl[1])

qzcar_comp = np.zeros((9,2))
for i in range(0,9):
    qzcar_comp[i][0] = contr['mag_contr']
    qzcar_comp[i][1] = contr['mag_contr']
np.savetxt('/home/alan/Desktop/dist_arcsec_IFS.txt', contr['distance_arcsec'], delimiter='   ')
np.savetxt('/home/alan/Desktop/mag_contr_IFS.txt', contr['mag_contr'], delimiter='   ')
np.savetxt('/home/alan/Desktop/mag_contr_K1.txt', contr['mag_contr'], delimiter='   ')

help(vip_hci.metrics.contrcurve.throughput)
vip_hci.metrics.contrcurve.contrast_curve(cube[0], -angs, psf_norm[0], fwhm=fwhm[0],
                                          pxscale=0.01225, starphot=1793119.1, algo=vip_hci.pca.pca, sigma=5, nbranch=1,
                                          theta=70, inner_rad=1, wedge=(0, 360), fc_snr=100,
                                          student=True, transmission=None, smooth=True,
                                          interp_order=2, plot=True, dpi=100, imlib='opencv',
                                          debug=True, verbose=True, full_output=False, save_plot=None,
                                          object_name=None, frame_size=None, figsize=(8, 4), ncomp=1,adimsdi='single',
                                          scale_list=wl[0])

K1 = np.array([])
K2 = np.array([])
with open('/Users/alan/Documents/PhD/Data/SPHERE/IRDIS/QZCar/Simplex_errors/Simplex_K2_flux_S3.txt') as f:
    for line in f:
        line = line.strip()
        columns = line.split()
        K2 = np.append(K2,float(columns[0]))
        #K1K2 = np.append(K1K2,float(columns[1]))

comp_xycoord = [[comp_pos[4][0],comp_pos[4][1]]] # Companion coords
simplex_guess_K1[4] = vip_hci.negfc.firstguess(cube[0],-angs,psf_norm[0],ncomp=ncomp_pca,plsc=pxscale,planets_xy_coord=comp_xycoord,fwhm=fwhm[0],annulus_width=ann_width,aperture_radius=aper_radius,simplex_options=simplex_options,f_range=f_range_K1[4],simplex=True,fmerit='sum',collapse='median',svd_mode='lapack',scaling=None,verbose=False,plot=False,save=False)

centy, centx = vip_hci.var.frame_center(cube[0,0])
posy = r_0 * np.sin(np.deg2rad(theta_0)) + centy
posx = r_0 * np.cos(np.deg2rad(theta_0)) + centx
posy = r_0 * np.sin(np.deg2rad(theta_0)) + 145
posx = r_0 * np.cos(np.deg2rad(theta_0)) + 145
np.deg2rad(theta_0) = mh.acos((60-145)/97)
print(posx, posy)

## Error estimation
r_0 = 1.971334773278342709e+02
theta_0 = 7.326447481557413255e+01
flux_0 = 3.223110891696441627e+01
centy, centx = vip_hci.var.frame_center(cube[0,0])
n_samples=7

theta_inj= np.zeros(n_samples)
r_inj=np.zeros(n_samples)
for j in range(0,n_samples):
   r_inj[j]=r_0
   theta_inj[j]=(j+1)*(360./(n_samples+1.))+theta_0-360.


spectra_mes_K1 = np.zeros(len(r_inj))
theta_mes_K1= np.zeros_like(theta_inj)
r_mes_K1= np.zeros_like(r_inj)
spectra_mes_K2 = np.zeros(len(r_inj))
theta_mes_K2= np.zeros_like(theta_inj)
r_mes_K2= np.zeros_like(r_inj)

for j in range(0,n_samples):
    print('Companion injection: ', j)
    posy = r_inj[j] * np.sin(np.deg2rad(theta_inj[j])) + centy
    posx = r_inj[j] * np.cos(np.deg2rad(theta_inj[j])) + centx

    cube_inj=vip_hci.metrics.cube_inject_companions(cube, psf_norm, -angs, flevel=flux_0, plsc=0.0074,
                                                   rad_dists=r_inj[j], n_branches=1, theta=theta_inj[j],
                                                   imlib='opencv', interpolation='lanczos4', verbose=False)
    planet_xycoord = np.array([[posx,posy]])
    print("success")
    fake_comp_K1 = vip_hci.negfc.simplex_optim.firstguess(cube_inj[0], -angs, psf_norm[0], ncomp=1, plsc=0.01225,
                                                      fwhm=fwhm[0], annulus_width=3, aperture_radius=3,
                                                      planets_xy_coord=planet_xycoord, cube_ref=None,
                                                      svd_mode='lapack',f_range=f_range_K1[11], simplex=True,
                                                      fmerit='sum',scaling=None, simplex_options=simplex_options,
                                                      collapse='median',verbose=False)

    fake_comp_K2 = vip_hci.negfc.simplex_optim.firstguess(cube_inj[1], -angs, psf_norm[1], ncomp=1, plsc=0.01225,
                                                      fwhm=fwhm[1], annulus_width=3, aperture_radius=3,
                                                      planets_xy_coord=planet_xycoord, cube_ref=None,
                                                      svd_mode='lapack',f_range=f_range_K2[11], simplex=True,
                                                      fmerit='sum',scaling=None, simplex_options=simplex_options,
                                                      collapse='median',verbose=False)

    r_mes_K1[j] = fake_comp_K1[0]
    theta_mes_K1[j] = fake_comp_K1[1]
    spectra_mes_K1[j] = fake_comp_K1[2]
    r_mes_K2[j] = fake_comp_K2[0]
    theta_mes_K2[j] = fake_comp_K2[1]
    spectra_mes_K2[j] = fake_comp_K2[2]

flux_K1 = 3.223110891696441627e+01
flux_K2 = 1.230475984428790959e+01
std_flux_K1=np.zeros_like(wl)
std_flux_corr_K1=np.zeros_like(wl)
std_r_K1=np.sqrt(np.sum(abs(r_inj - r_mes_K1)**2)/(n_samples-1))
std_theta_K1=np.sqrt(np.sum(abs(theta_inj - theta_mes_K1)**2)/(n_samples-1))
std_flux_K1=np.sqrt(np.sum(abs(spectra_mes_K1 - flux_K1)**2)/(n_samples-1))
std_flux_K2=np.zeros_like(wl)
std_flux_corr_K2=np.zeros_like(wl)
std_r_K2=np.sqrt(np.sum(abs(r_inj - r_mes_K2)**2)/(n_samples-1))
std_theta_K2=np.sqrt(np.sum(abs(theta_inj - theta_mes_K2)**2)/(n_samples-1))
std_flux_K2=np.sqrt(np.sum(abs(spectra_mes_K2 - flux_K2)**2)/(n_samples-1))

print(std_r_K1)

std_flux_K2=np.sqrt(np.sum(abs(K2 - 3.196914362045397695e+01)**2)/(13-1))


## Read MCMC
import pickle # Import important MCMC libraries
from pickle import Pickler
pickler={}
mcmc_result={}
for i in range(0,2): # Read all channels and store them to variables
    with open('/Users/alan/Documents/PhD/Data/SPHERE/IRDIS/QZCar/MCMC/S9S10/QZCarK1_IRDIS_companions_S{}/MCMC_results'.format(i),'rb') as fi:
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
for i in range(0,2):
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
R2 = 22**2 # Radius of the star considered
mask = np.zeros_like(cube[:,0,:,:])
mask_value = np.zeros((2,44,44))
shift = -angs+angs[0] # Angular shift
for h in range(len(angs)):
    x=0
    y=0
    for i in range(300,344): # Iterate through all pixels in X/Y positions
        for j in range(394,438):
            r2 = (i-322)**2+(j-416)**2# Calculate the radius of the pixel from the center pixel
            if r2<=R2: # If the pixels are in the circle of radius R
                mask_value[:,y,x] = cube[:,h,j,i]
            y+=1
        y=0
        x+=1
    x=0
    y=0
    for i in range(326,370): # Iterate through all pixels in X/Y positions
        for j in range(356,400):
            r2 = (i-348)**2+(j-378)**2# Calculate the radius of the pixel from the center pixel
            if r2<=R2: # If the pixels are in the circle of radius R
                cube_mask[:,h,j,i] = mask_value[:,y,x]
            y+=1
        y=0
        x+=1
