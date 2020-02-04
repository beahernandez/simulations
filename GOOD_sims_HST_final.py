
# coding: utf-8

# In[1]:

import sys
import os
import math
import numpy
import galsim
import numpy as np
from astropy.io import ascii
from astropy.io import fits
from scipy import stats

# In[2]:


def generate_rand_from_pdf(pdf, x_grid,n):
    cdf = np.cumsum(pdf)
    cdf = cdf / cdf[-1]
    values = np.random.rand(n)
    value_bins = np.searchsorted(cdf, values)
    random_from_cdf = x_grid[value_bins]
    return random_from_cdf

def kde(x, x_grid, bandwidth=0.2, **kwargs):
    """Kernel Density Estimation with Scipy"""
    kde = stats.gaussian_kde(x, bw_method=bandwidth / x.std(ddof=1), **kwargs)
    return kde.evaluate(x_grid)

sky_level = 40. #30. #5. #0.8   #euclid
pixel_scale = 0.05               # arcsec / pixel
noise_variance = (sky_level*pixel_scale**2)**2            # ADU^2
gain=14 #e-/ADU


gal_flux_min =2 # 3. #1.e4     # Range for galaxy flux
gal_flux_max =200.# 3.e1 #1.e6  
gal_hlr_min = 0.6       # arcsec
gal_hlr_max = 1.3       # arcsec
gal_e_min = 0.          # Range for ellipticity
gal_e_max = 0.8

psf_fwhm = 0.1 #0.65         # arcsec
psf_beta = 5 #Moffat profile

nx_tiles =100                   #
ny_tiles =100                   #
stamp_xsize = 100                #
stamp_ysize = 100                #

gal_ellip_max = 0.8 
gal_ellip_rms = 0.3 

disk_n = 1.5
disk_r0 = 0.85
#2.3

#random_seed = 6424512
#gal_resolution = 0.98 
gal_re = 0.1 #0.85

gal_g1=0.05
gal_g2=0.07

shift_radius_sq = 1.0 
print 'done'

#SB Error: fourierDraw() requires an FFT that is too large, 6144
#If you can handle the large FFT, you may update gsparams.maximum_fft_size.
# In[ ]:

random_seed = 6424532 #12
initial_n=0
gal_g1=-0.4 + initial_n*0.016
gal_g2=-0.45+ initial_n*0.018
#cat_file_name ='galsim_default_input.asc'
#cat = galsim.Catalog(cat_file_name)
shears1=[]
shears2=[]
numbers=[]
ellipticities=[]
size=[]
flux=[]
mag=[]
sersic_indices=[]
print 'hey'
indices=[]

gara=galsim.GaussianDeviate(mean=1.,sigma=1) #random_seed,

############################gasn=galsim.GaussianDeviate(random_seed,mean=25.5,sigma=0.8) #gasn=galsim.GaussianDeviate(random_seed,mean=26,sigma=1.) #Currently draws for magnitudes. #0.5

ref_data = ascii.read('/vol/aibn148/data1/bhernandez/PhD/PSF_HST/tmp.regions_VI3')
goodgal = np.array([ref_data[i][2] for i in range(len(ref_data))])
goodregion = np.array([ref_data[i][1] for i in range(len(ref_data))])
VI = np.array([ref_data[i][5] for i in range(len(ref_data))])
#VI = VI[(goodgal>0)&(goodregion>0)]

#ref_rad_tmp = [ref_data[i][4] for i in range(len(ref_data))]
#ref_rad_tmp = np.array(ref_rad_tmp)
#ref_rad = ref_rad_tmp[(goodgal>0.)&(goodregion>0.)&(VI<0.3)]
ref_sn_tmp = [ref_data[i][3] for i in range(len(ref_data))]
ref_sn_tmp = np.array(ref_sn_tmp)
ref_sn = ref_sn_tmp[(goodgal>0.)&(goodregion>0.)&(VI<0.3)]
ref_mag_tmp = [ref_data[i][0] for i in range(len(ref_data))]
ref_mag_tmp = np.array(ref_mag_tmp)
ref_mag = ref_mag_tmp[(goodgal>0.)&(goodregion>0.)&(VI<0.3) & (ref_mag_tmp>-50)]
    
x_grid_mag=np.linspace(min(ref_mag), max(ref_mag), 1000)
MAG = kde(ref_mag, x_grid_mag, bandwidth=0.1) #stats.gaussian_kde(ref_mag, bw_method=0.1 / ref_mag.std(ddof=1))
random_mag= generate_rand_from_pdf(MAG, x_grid_mag,1)


#gara=galsim.GaussianDeviate(mean=3,sigma=3.)
#gasn=galsim.GaussianDeviate(mean=1.,sigma=2.)
#From fitting:
#A=2.02261391351e-08
#Zp=25.3513858337
A=2.52722419076 
Zp=27.1007978301
#####Magnitude calibration
####Zp=26.7 #For F814W
####t_exp=100.#10. #In seconds
####gain=14. #HST like e-/ADU

######I_ref=t_exp/gain*10**(-0.4*(0.)) #In ADU because of gain
######m_ref=Zp
######I_sky=sky_level*pixel_scale**2
######m_sky = 5. / 2. * np.log10(I_ref/I_sky ) + m_ref
######m_sky=22.35 #mag arcsec−2

######sky_level= A*10**(-0.4*(m_sky-Zp)) #Per arcsec−2 #/pixel_scale**2
######noise_variance = (sky_level*pixel_scale**2)**2            # In real units?^2

#print m_sky

#####cat_file_name = '/vol/euclid1/euclid1_3/bhernandez/galsims/COSMOS_23.5_training_sample/real_galaxy_catalog_23.5.fits'
######'/vol/euclid1/euclid1_3/bhernandez/galsims/test/real_galaxy_catalog_23.5_example_sel.fits'
#####real_galaxy_catalog = galsim.RealGalaxyCatalog(cat_file_name)

#####cat=fits.getdata(cat_file_name)
#####id_list = cat['IDENT']

big_fft_params = galsim.GSParams(maximum_fft_size=30000)


ud = galsim.UniformDeviate() #random_seed
gd = galsim.GaussianDeviate(ud, sigma=gal_ellip_rms)

gd_sersic = galsim.GaussianDeviate( mean=1., sigma=1) #random_seed,


#psf_model=fits.getdata('/export/data1/bhernandez/PhD/PSF_HST/Tiny_Tim/result00_psf.fits')
psf= galsim.InterpolatedImage('/vol/aibn148/data1/bhernandez/PhD/PSF_HST/test_bg3_cdgstridegood.fits', gsparams=big_fft_params, scale=0.0165)#focus/zero_cdgstridegood.fits', gsparams=big_fft_params, scale=0.0165) #1.*0.0231) #Pixel scale from header of image
#psf=psf.rotate(90*galsim.degrees)
print psf.calculateFWHM(scale=pixel_scale)
#psf= galsim.InterpolatedImage('/export/data1/bhernandez/PhD/PSF_HST/Tiny_Tim/result00.fits', gsparams=big_fft_params, scale=0.0231) #1.*0.0231) #Pixel scale from header of image
#print psf.calculateFWHM(scale=pixel_scale)
psf=psf.dilate((psf_fwhm)/psf.calculateFWHM(scale=pixel_scale)) #/pixel_scale
#psf = galsim.Gaussian(fwhm = psf_fwhm, gsparams=big_fft_params)
#psf = galsim.Moffat(beta=psf_beta, flux=1., fwhm=psf_fwhm)

#bp_file = os.path.join(galsim.meta_data.share_dir, 'wfc_F814W.dat.gz')
#bandpass = galsim.Bandpass(bp_file, wave_type='ang').thin().withZeropoint(25.94)
cosmos_cat = galsim.COSMOSCatalog()#dir=/vol/aibn148/data1/bhernandez/share/galsim/COSMOS_25.2_training_sample)

for n in range(initial_n,16): #50):# #30):#43):
    x=[]
    y=[]
    ellipticities=[]
    betas=[]
    size=[]
    flux=[]
    mag=[]
    sersic_indices=[]
    tmp_indices=[]
    images = []
    
    
    indices = np.random.random_integers(0, high=81490, size=ny_tiles*nx_tiles) #np.arange()
    im_size = 64
    real_gal_list = cosmos_cat.makeGalaxy(indices, gal_type='parametric', gsparams=big_fft_params)#, chromatic=True) #gal_type='real', noise_pad_size=im_size*pixel_scale)

    gal_image = galsim.ImageF(stamp_xsize * nx_tiles-1 , stamp_ysize * ny_tiles-1,
                              scale=pixel_scale)
    gal_image_rot = galsim.ImageF(stamp_xsize * nx_tiles-1 , stamp_ysize * ny_tiles-1,
                              scale=pixel_scale)

    psf_image = galsim.ImageF(stamp_xsize * nx_tiles-1 , stamp_ysize * ny_tiles-1,
                              scale=pixel_scale)
    k=1
    

    gal_g1+=0.016
    gal_g2+=0.018
    print gal_g1, gal_g2
    shears1.append(gal_g1)
    shears2.append(gal_g2)
    numbers.append('%03d'%n)
    real_idx=-1
    for iy in range(ny_tiles):
        for ix in range(nx_tiles):
            
            
             #beta =ud() * 2. * math.pi * galsim.radians
            val = gd()
            ellip = math.fabs(val)
            if ellip>gal_ellip_max:
                ellip=gal_ellip_max
            #print ellip
            #ellipticities.append(ellip)
            b = galsim.BoundsI(ix*stamp_xsize +1, (ix+1)*stamp_xsize -1, 
                               iy*stamp_ysize +1, (iy+1)*stamp_ysize-1)
            sub_gal_image = gal_image[b]
            sub_gal_image_rot = gal_image_rot[b]
            sub_psf_image = psf_image[b]
        
            rsq = 2 * shift_radius_sq #Random shift
            while (rsq > shift_radius_sq):
                dx = (2*ud()-1) * shift_radius_sq
                dy = (2*ud()-1) * shift_radius_sq
                rsq = dx**2 + dy**2 
        
            
        
            #psf_re = psf.getHalfLightRadius()
            #gal_re = psf_re * gal_resolution
            
            reu = np.abs(gara())# 0.2+ud()*()
            re=gal_re*reu
            #print reu*disk_r0
        
            m=np.float(generate_rand_from_pdf(MAG, x_grid_mag,1))+2#+0.8 #-0.8 #gasn() #gal_flux_min+ud()*(gal_flux_max-gal_flux_min)
            
            gal_flux=A*10**(-0.4*(m-Zp))
            #gal_flux=t_exp/gain*10**(-0.4*(m-Zp)) #In ADU because of gain
            #######sersic_idx=gd_sersic()#sersic_idx=0.5+ud()*(4.-0.5)
            
            #######while sersic_idx>6.1 or sersic_idx<0.3:
                #######sersic_idx=gd_sersic()
            
            
            #gal = galsim.Sersic(sersic_idx, flux=gal_flux, half_light_radius=re)
            real_idx+=1
            gal=real_gal_list[real_idx]
            #gal = galsim.RealGalaxy(real_galaxy_catalog,flux=gal_flux,id=id_list[56])#id_list[real_idx])
            param=cosmos_cat.getParametricRecord(indices[real_idx])
            hlr=param['hlr'][0]
            #print re,hlr
            sersic_idx=param['sersicfit'][2]
            beta =param['sersicfit'][7]
            ellip = param['sersicfit'][3]
            try:
                #gal = gal.dilate(50.)
                gal = gal.dilate(re/hlr)
                #print real_idx
            except:
                print real_idx, 'Too large'
                
            gal = gal.withFlux(gal_flux)
            half_light_radius=re
            gal_mag=m
            #print gal_flux,gal_mag
            tmp_indices.append(sersic_idx)
            
            #gal = galsim.DeVaucouleurs(flux=gal_flux, half_light_radius=re)
            gal = gal.shear(e=ellip, beta=beta*galsim.radians)
            
            gal_rot=gal
            gal_rot=gal_rot.rotate(90*galsim.degrees)

        
            gal = gal.shear(g1=gal_g1, g2=gal_g2) #Add shear
            #gal = gal.shift(dx,dy)
            gal_rot = gal_rot.shear(g1=gal_g1, g2=gal_g2) #Add shear
    
            final = galsim.Convolve([psf, gal], real_space=False) #Final galaxy and psf
            final = final.shift(dx,dy)
        
            final_rot = galsim.Convolve([psf, gal_rot], real_space=False) #Final galaxy and psf
            final_rot = final_rot.shift(dx,dy)
        
            psf_flux=gal_flux_min*10**(ud()*(np.log10(gal_flux_max/gal_flux_min)))
            #print psf_flux
            star = psf.withFlux(psf_flux)
            star = star.shift(dx,dy)
        
            sub_gal_image += sky_level * pixel_scale**2
            try:
              final.drawImage(image=sub_gal_image, scale=pixel_scale)  #Draws image
        
              final_rot.drawImage(image=sub_gal_image_rot, scale=pixel_scale)  #Draws image
            except:
              k-=1
        
            star.drawImage(sub_psf_image, scale=pixel_scale)
        
            x.append(ix*stamp_xsize+stamp_xsize/2.+dx)
            y.append(iy*stamp_ysize+stamp_ysize/2.+dy)
            ellipticities.append(ellip)         
            betas.append(beta)
            size.append(re)
            sersic_indices.append(sersic_idx)
            flux.append(gal_flux)
            mag.append(gal_mag)
            
   
            k+=1
    
    #t_image=[ellipticities,size,sersic_indices,flux,mag] #, ellipticities]
    t_image=[x,y,ellipticities,betas,size,sersic_indices,flux,mag]
    ascii.write(t_image, 'HST_bg_ellip_10000_sims%03d_galaxies_00_input.txt'%n)
    noise = galsim.GaussianNoise(ud, sigma=math.sqrt(noise_variance))
    noise_psf = galsim.GaussianNoise(ud, sigma=math.sqrt(noise_variance))
    #noise = galsim.CCDNoise(ud, gain=gain, read_noise=read_noise)
    gal_image.addNoise(noise)
    gal_image_rot.addNoise(noise)
    psf_image.addNoise(noise)
    if n==0:
        psf_file_name='HST_bg_ellip_10000_sims_starfield.fits'
        psf_image.write(psf_file_name)
        
    gal_file_name='HST_bg_ellip_10000_sims%03d_galaxies_00.fits'%n
    gal_file_name_rot='HST_bg_ellip_10000_sims%03d_galaxies_90.fits'%n

    
    gal_image_rot.write(gal_file_name_rot)
    gal_image.write(gal_file_name)   
    
    #indices.append(tmp_indices)

    
#ascii.write(indices)


# In[17]:

######tables=[numbers,shears1,shears2] #, ellipticities]
######ascii.write(tables, 'input_shear_HST_realparam_TTsubsampled_finetuningmag_10000.txt')



# In[ ]:


