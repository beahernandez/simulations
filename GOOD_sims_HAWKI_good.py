
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

#Input properties
sky_level = 40. #5. #0.8   #euclid
pixel_scale = 0.106 #HAWKI              # arcsec / pixel
noise_variance = (sky_level*pixel_scale**2)**2            # ADU^2
gain=14 #e-/ADU


gal_flux_min =2 # 3. #1.e4     # Range for galaxy flux
gal_flux_max =200.# 3.e1 #1.e6  
gal_hlr_min = 0.6       # arcsec
gal_hlr_max = 1.3       # arcsec
gal_e_min = 0.          # Range for ellipticity
gal_e_max = 0.8

psf_fwhm = 0.4 #HAWKI #0.65         # arcsec
psf_beta = 5 #Moffat profile

nx_tiles =70                   #
ny_tiles =70                   #
stamp_xsize = 100                #
stamp_ysize = 100                #

gal_ellip_max = 0.8 
gal_ellip_rms = 0.2 #0.3 

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


random_seed = 6424532 #12
initial_n=4
gal_g1=-0.4 + initial_n*0.016
gal_g2=-0.45+ initial_n*0.018
#cat_file_name ='galsim_default_input.asc'
#cat = galsim.Catalog(cat_file_name)

#initialization
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
big_fft_params = galsim.GSParams(maximum_fft_size=30000)

#tmp=ascii.read('/vol/euclid1/euclid1_3/bhernandez/galsims/HAWKI/setup/good_catalogue.cat')
#hawki_cat= tmp[(tmp['col6']>21.0) & (tmp['col6']<24.2) & (tmp['col7']>10.) & (np.sqrt(tmp['col3'] **2 +tmp['col4']**2)<1.0)]
#select_cat=fits.getdata('/vol/euclid1/euclid1_3/bhernandez/galsims/HAWKI/setup/select_colorselected.cat')



#mag_list=hawki_cat['col6']

#Opening the catalogues
input_cat=ascii.read('/vol/euclid1/euclid1_3/bhernandez/galsims/HAWKI/setup/cosmos/cosmos_3dhst.v4.1_f160w.galfit')
index_list="/vol/euclid1/euclid1_3/bhernandez/galsims/HAWKI/setup/color_cut_indices.txt"
f = open(index_list, 'r') # 'r' = read
cat_index = np.genfromtxt(f, delimiter='\n',dtype=int)
f.close()

select_cat=input_cat[cat_index]
select_cat=select_cat[(select_cat['RA']>150.2) & (select_cat['DEC']>2.2)]
input_cat=input_cat[(input_cat['RA']>150.2) & (input_cat['DEC']>2.2)]
print len(input_cat)

A=2.52722419076 
Zp=27.1007978301

RA_max=max(input_cat['RA'])
RA_min=min(input_cat['RA'])
DEC_max=max(input_cat['DEC'])
DEC_min=min(input_cat['DEC'])

candels_size_x=int((RA_max-RA_min)*3600./0.05) #11217
candels_size_y=int((DEC_max-DEC_min)*3600./0.05)
big_fft_params = galsim.GSParams(maximum_fft_size=30000)

#Initialization of the random objects
ud = galsim.UniformDeviate() #random_seed
gd = galsim.GaussianDeviate(ud, sigma=gal_ellip_rms)

gd_sersic = galsim.GaussianDeviate( mean=1., sigma=1) #random_seed,


psf = galsim.Moffat(beta=psf_beta, flux=1., fwhm=psf_fwhm)


#Loop over all images
for n in range(initial_n,50):
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
    flag=[]
    
    
    indices = np.random.random_integers(0, high=81490, size=ny_tiles*nx_tiles) #np.arange()
    im_size = 64
    ##real_gal_list = cosmos_cat.makeGalaxy(indices, gal_type='parametric', gsparams=big_fft_params)#, chromatic=True) #gal_type='real', noise_pad_size=im_size*pixel_scale)
    
    #Creation of the empty images
    gal_image = galsim.ImageF(candels_size_x ,candels_size_y,
                              scale=pixel_scale)
    gal_image_rot = galsim.ImageF(candels_size_x , candels_size_y,
                              scale=pixel_scale)

    psf_image = galsim.ImageF(candels_size_x , candels_size_y,
                              scale=pixel_scale)
    k=1
    
    #Shear
    gal_g1+=0.016
    gal_g2+=0.018
    print gal_g1, gal_g2
    shears1.append(gal_g1)
    shears2.append(gal_g2)
    numbers.append('%03d'%n)
    real_idx=-1
    
    #Loop over the galaxies
    for obj in range(len(input_cat)):
        #print obj
	
	ra=input_cat['RA'][obj]
	dec=input_cat['DEC'][obj]
	
        #Check if the galaxies are in the colour-selected catalogue
	if np.where((ra-0.0005<select_cat['RA']) & (select_cat['RA']<ra+0.0005) & (dec-0.0005<select_cat['DEC'])&(select_cat['DEC']<dec+0.0005))>1: 
            flag.append(1)
        else:
            flag.append(0)

        val = gd()
        ellip = math.fabs(val)
        if ellip>gal_ellip_max:
            ellip=gal_ellip_max

        b = galsim.BoundsI(1, stamp_xsize, 
                        1, stamp_ysize)
        sub_psf_image = psf_image[b]
        
        x_tmp = (ra-RA_min)*3600./0.05#ud()*(image_size-stamp_xsize)
        y_tmp = (dec-DEC_min)*3600./0.05#ud()*(image_size-stamp_xsize)
        image_pos = galsim.PositionD(x_tmp,y_tmp)
        
        #Random shift on the pixel level
        rsq = 2 * shift_radius_sq #Random shift
        while (rsq > shift_radius_sq):
            dx = (2*ud()-1) * shift_radius_sq
            dy = (2*ud()-1) * shift_radius_sq
            rsq = dx**2 + dy**2 

        #re=re_list[idx_cats]

        m=input_cat['mag'][obj]#mag_list[np.int(np.random.random_integers(0, high=max_idx_m,size=1))] #idx_cats]#
        gal_flux=A*10**(-0.4*(m-Zp))
        
        #Sersic index
        re=input_cat['re'][obj]
        sersic_idx=input_cat['n'][obj]
        if sersic_idx>6.2:
            sersic_idx=6.2
        elif sersic_idx<0.3:
            sersic_idx=0.3
            
        #sersic_cat=tmp[(tmp['n']>0.3) & (tmp['n']<6.2)]
        
        gal = galsim.Sersic(sersic_idx, flux=gal_flux, half_light_radius=re)
        real_idx+=1
 
        hlr=input_cat['re'][obj]
        beta=input_cat['pa'][obj]* galsim.radians
        ellip=input_cat['q'][obj]
        if ellip>0.8:
            ellip=0.8
        elif ellip<-0:
            ellip=0
 
            
        gal = gal.withFlux(gal_flux)
        half_light_radius=re
        gal_mag=m
        #print gal_flux,gal_mag
        tmp_indices.append(sersic_idx)
        

        gal = gal.shear(e=ellip, beta=beta) #Introduces the ellipticity
        
        gal_rot=gal #Creates the rotated copy
        gal_rot=gal_rot.rotate(90*galsim.degrees)

    
        gal = gal.shear(g1=gal_g1, g2=gal_g2) #Add shear
        #gal = gal.shift(dx,dy)
        gal_rot = gal_rot.shear(g1=gal_g1, g2=gal_g2) #Add shear

        final = galsim.Convolve([psf, gal], real_space=False, gsparams=big_fft_params) #Final galaxy and psf
        #final = final.shift(dx,dy)
    
        final_rot = galsim.Convolve([psf, gal_rot], real_space=False, gsparams=big_fft_params) #Final galaxy and psf
        #final_rot = final_rot.shift(dx,dy)
    
        psf_flux=gal_flux_min*10**(ud()*(np.log10(gal_flux_max/gal_flux_min)))
        #print psf_flux
        star = psf.withFlux(psf_flux) #Creates the PSF image
        star = star.shift(dx,dy)
    
        #Saving the objects into the images
        try:
        #sub_gal_image += sky_level * pixel_scale**2
            stamp = final.drawImage(scale=pixel_scale) 
            stamp_rot = final_rot.drawImage(scale=pixel_scale)
        

            stamp.setCenter(image_pos.x,image_pos.y)
            stamp_rot.setCenter(image_pos.x,image_pos.y)
        
            bounds = stamp.bounds & gal_image.bounds
            bounds_rot = stamp_rot.bounds & gal_image_rot.bounds

            # Finally, add the stamp to the full image.
            gal_image[bounds] += stamp[bounds]
            gal_image_rot[bounds_rot] += stamp_rot[bounds_rot]
            
            star.drawImage(sub_psf_image, scale=pixel_scale)
        
        except:
            print 'Error in drawing'
    
        #Adds the parameters to the input file
        x.append(image_pos.x)
        y.append(image_pos.y)
        ellipticities.append(ellip)         
        betas.append(beta)
        size.append(re)
        sersic_indices.append(sersic_idx)
        flux.append(gal_flux)
        mag.append(gal_mag)
        
        

        k+=1
    
    #Saves images and files
    #t_image=[ellipticities,size,sersic_indices,flux,mag] #, ellipticities]
    t_image=[x,y,ellipticities,betas,size,sersic_indices,flux,mag,flag]
    ascii.write(t_image, 'HAWKI_realcats_cosmos1_sims%03d_galaxies_00_input.txt'%n)
    noise = galsim.GaussianNoise(ud, sigma=math.sqrt(noise_variance))
    noise_psf = galsim.GaussianNoise(ud, sigma=math.sqrt(noise_variance))
    #noise = galsim.CCDNoise(ud, gain=gain, read_noise=read_noise)
    gal_image.addNoise(noise)
    gal_image_rot.addNoise(noise)
    psf_image.addNoise(noise)
    if n==0:
        psf_file_name='HAWKI_realcats_cosmos1_sims_starfield.fits'
        psf_image.write(psf_file_name)
        
    gal_file_name='HAWKI_realcats_cosmos1_sims%03d_galaxies_00.fits'%n
    gal_file_name_rot='HAWKI_realcats_cosmos1_sims%03d_galaxies_90.fits'%n

    
    gal_image_rot.write(gal_file_name_rot)
    gal_image.write(gal_file_name)   
    



