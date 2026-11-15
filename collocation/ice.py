
# coding: utf-8

# In[1]:

#! /usr/bin/env python

#==================
# Open sea ice product; create a mask from SSMI data; plot SSMI chart (resolution 12.5 km), 

# Important!
# - date should be in format '2015-10-28'


# authors: Anastasia Tarasenko (SIAM LOPS Ifremer) 
# creation: oct 2016
#==================

import os, sys
sys.path.append('/home/cercache/users/atarasen/codes/python/')
import ncdfView as nv
import numpy as np
import numpy.ma as ma
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap


# In[2]:

def SSMI_ice_open(date):
  path_grid_0='/home/cersat9/cersat/ftp/products/gridded/psi-concentration/data/grid_north_12km.nc.gz'
  path_grid='/home/cercache/users/atarasen/data/SSMI/'+path_grid_0[-21:-3]
  os.system('gunzip -c %s > %s'%(path_grid_0, path_grid))  
  nc_grid=nv.ncdfView(path_grid, Quiet=True)
  lon_ice=nc_grid('longitude')
  lat_ice=nc_grid('latitude')
  
  path_sea_ice_Z='/home/cersat9/cersat/ftp/products/gridded/psi-concentration/data/arctic/daily/netcdf/'+date[0:4]+'/'
  ice_conc=date[0:4]+date[5:7]+date[8:10]+'.nc.Z'
  
  path_sea_ice='/home/cercache/users/atarasen/data/SSMI/'+date[0:4]+'/'
  if os.path.isdir(path_sea_ice)==False:
    os.system('mkdir %s' %(path_sea_ice))
    
  path_sea_ice_nc=path_sea_ice+ice_conc[:-2]
  print path_sea_ice_nc
  if os.path.isfile(path_sea_ice_nc)==False:
    os.system('gunzip -c %s > %s'%(path_sea_ice_Z+ice_conc, path_sea_ice_nc))
    print 'path_sea_ice_nc unzipped'
  
  if os.stat(path_sea_ice_nc).st_size==0:
    ice=[]
  else:
    nc  = nv.ncdfView(path_sea_ice_nc, Quiet=True)
    ice=nc('concentration')
  return ice, lon_ice, lat_ice


# In[3]:

def SSMI_ice_mask(ice, percent=5):
  ice_masked = ma.masked_less(ice, percent) #the outputs are: ice_masked [55, 60, --, ... ] and ice_masked.mask ([true, true, false...])
  return ice_masked


# In[4]:

def SSMI_ice_show(date, ice=[], lon_ice=[], lat_ice=[]):
  path_sea_ice='/home/cercache/users/atarasen/data/SSMI/'+date[0:4]+'/'
  if ice==[]:
    ice, lon_ice, lat_ice=SSMI_ice_open(date)
    
  # setup Basemap
  m = Basemap(projection='npstere',boundinglat=lat_ice[0][0],lon_0=0,resolution='l')
  m.drawcoastlines()
  m.fillcontinents(color='coral',lake_color='aqua')
  # draw parallels and meridians.
  m.drawparallels(np.arange(-80.,81.,10.))
  m.drawmeridians(np.arange(-180.,181.,20.))
  m.drawmapboundary(fill_color='aqua')
  
  #im1 = m.pcolormesh(lon_ice, lat_ice, ice, shading='flat', cmap=plt.cm.jet)
  #plt.figure(1, figsize=(12,8), dpi=300)
  m.imshow(ice)
  plt.show()
  plt.savefig(path_sea_ice+date+'.png', dpi=300)
  
  


# In[6]:


def SSMI_ice_buoy_show(date, ice, lon_ice, lat_ice, lon_buoy, lat_buoy):
 path_sea_ice='/home/cercache/users/atarasen/data/SSMI/'+date[0:4]+'/'

 plt.figure(1, figsize=(12,8), dpi=300)
 plt.imshow(ice)
 
 plt.savefig(path_sea_ice+date+'.png', dpi=300)  


# In[5]:

def osisaf_icemask_open(date):
  path_sea_ice_mask='/home2/taveeg/data/operation/project/myocean/sst-tac/odyssea/v2/glob/mask/'+start_d[0:4]+'/'
  ice_mask=path_sea_ice_mask+start_d[0:4]+start_d[5:7]+start_d[8:10]+'-EUR-L4HRmask-GLOB-v02-fv01-ODYSSEA.nc'
  print  ice_mask
  nc  = nv.ncdfView(ice_mask, Quiet=True)
  lat_ice=nc('lat')
  lon_ice=nc('lon')
  ice_mask=nc('mask')
  ice_mask=np.squeeze(ice_mask) #get rid of non-used dimension
  return ice_mask, lon_ice, lat_ice

