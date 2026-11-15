### Set of basic python functions to work with Sentinel-1 image using posada
# author: Tarasenko Anastasia with posada code (Solab),
# date: Jul-Aug 2017 


import numpy as np

import scipy as sp
import scipy.signal as sg
import scipy.ndimage as ndimage
from scipy.interpolate import interp1d

import sys, os

# sys.path.append('/home/mag/Documents/repos/solab/POSADA/') #common version on student server
## sys.path.append('/home/tarasenko/codes/posada3/') #was old version (before switching to posada3 on python3)
sys.path.append('/home/tarasenko/prj/posada3/') #last version on my computerin St.P.

from posada.Handlers import Sentinel_1 as SAR
from posada.Models import Swath
from posada.Readers import SAFE

from posada.Handlers.resources import sentinel_koeffs as sk
from posada.Handlers.sentinel_1 import change_thermal_noise
from posada.Handlers import basesar



from cmod.cmod_gpu import rcs2windOpenCl
from cmod.cmod_vect import rcs2windPar



def get_values_without_thermalnoise(url, sigma0, inc):
# created from  Sentinel_1(BaseSAR)._get_values_nice_image_hh_vv(self, slices=None, **kwargs) in posada/Handlers/sentinel-1.py
    tnoise_new, flatline = change_thermal_noise(sk.tnoise, sk.incvec)
    fnoise_new = interp1d(sk.incvec, tnoise_new, kind='slinear')

    sigma0 = np.ma.array(sigma0, mask=np.zeros_like(sigma0))
    sigma0_corr = sigma0 - fnoise_new(inc)

    return tnoise_new, sk.incvec, sigma0_corr




def get_data(url, res_m = 500, pol = 'hh', wind_model = 'interim'):
# Get data from Sentinel1 image
# sar.get_field_names() 
    sarf = SAFE(url=url, polarisation=pol)
    sar = Swath()
    sar.load(sarf)
    sar.add_handler(SAR)
    sar.set_spacing(sar.resolution2spacing(res_m))

    sigma0 = sar.get_values('sigma0')
    inc = sar.get_values('incidence')
    lat = sar.get_values('lat')
    lon = sar.get_values('lon')

    # Denoised sigma_0     
    tnoise_n, incvec, sigma0_corr = get_values_without_thermalnoise(url, sigma0, inc)  
       
    # Nice images from Albert Murtazin
    sigma0_nice = sar.get_values('nice_image')
    
    # Extract reanalysis wind with posada
    if wind_model =='interim':
        wdir_model = sar.get_values('wind_dir_era_interim')
    else:
        wdir_model = sar.get_values('wind_dir_ncep_gfs')

        
    if lat[0,0]>lat[-1,-1]:
         trj = 'desc'
    else:
         trj = 'asc'
# #rotate everything to get good vectors for quiver
#         sigma0 = np.flipud(np.fliplr(sigma0))
#         sigma0_corr = np.flipud(np.fliplr(sigma0_corr))
#         sigma0_nice = np.flipud(np.fliplr(sigma0_nice))
        
#         wdir_model = np.flipud(np.fliplr(wdir_model))
         wdir_model = wdir_model-180

#         lat = np.flipud(np.fliplr(lat))
#         lon = np.flipud(np.fliplr(lon))
    print (trj +' trajectory')           
        
        
    return sigma0, lat, lon, inc, sigma0_corr, sigma0_nice, wdir_model 



def get_wind_data(url, sigma0, inc, wdir, lat, lon):
 
    # this is part of the function 'compute_wind_speed' in basesar.py
    PR = basesar.PR_Mouche(inc, wdir)
    try:
        wind_speed = rcs2windOpenCl(sar=sigma0 * PR, windir=wdir, theta=inc)
    except Exception:
        wind_speed = rcs2windPar(sar=sigma0 * PR, cmdv=5, windir=wdir, theta=inc, nprocs=4)
  
    return wind_speed



def get_ice_mask_draft(sigma_hh, sigma_hv, inc, resolution):
    fn = 'data/tnoise_new/tnoise_min.mat'
    #fn = 'data/tnoise_new/tnoise_min_HStreaks.mat'
    tnoise_min = np.array(sp.io.loadmat(fn, squeeze_me=True)['tnoise'])
    incvec = np.array(sp.io.loadmat(fn, squeeze_me=True)['incvec'])

    fn = 'data/tnoise_new/tnoise_mean.mat'
    tnoise_mean = np.array(sp.io.loadmat(fn, squeeze_me=True)['tnoise'])
    incvec = np.array(sp.io.loadmat(fn, squeeze_me=True)['incvec'])

