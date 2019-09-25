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
