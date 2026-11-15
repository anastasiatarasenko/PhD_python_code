#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 28 01:12:11 2019

@author: anastasiiatarasenko
"""
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import os, sys

import conda
conda_file_dir = conda.__file__
os.environ["PROJ_LIB"] = os.path.join(os.path.join(conda_file_dir.split('lib')[0], 'share'), 'proj')
from mpl_toolkits.basemap import Basemap

from matplotlib.colors import LinearSegmentedColormap
import matplotlib.lines as mlines

sys.path.append('/Users/anastasiiatarasenko/codes/satellite/')
import open_sst_L4

  
bbox = [74, 85, 90, 170]
latmin, latmax, lonmin, lonmax = bbox
reg_lat, reg_lon = np.meshgrid(np.arange(latmin, latmax, 0.01), 
                               np.arange(lonmin, lonmax, 0.05))


  
def discrete_cmap(N, base_cmap=None, base = None):
    """
    Create an N-bin discrete colormap from the specified input cmap
    INPUT
        either base_cmap (string - cmap name) or base (colormap instance) 
        should be provided
    OUTPUT
        descrete colormap
    """
    
# Note that if base_cmap is a string or None, you can simply do
#    return plt.cm.get_cmap(base_cmap, N)
# The following works for string, None, or a colormap instance:
    if base_cmap:
        base = plt.cm.get_cmap(base_cmap)
    
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    return LinearSegmentedColormap.from_list(cmap_name, color_list, N)



def search_bbox(bbox, lat, lon):
    latmin, latmax, lonmin, lonmax = bbox

    ma_lon = np.where(np.logical_or(lon<=lonmin, lon>=lonmax), np.ones(lon.shape), 0)
    ma_lat = np.where(np.logical_or(lat<=latmin, lat>=latmax), np.ones(lat.shape), 0)

    prep_mask = ma_lat+ma_lon
    mask = np.where(prep_mask==0, 0, 1)
    mask = np.array(mask, dtype=bool)
    inds = np.where(mask==True)
    return mask, inds


def get_topo(etopo = '/Users/anastasiiatarasenko/data/etopo1_3_Kara-Laptev-EastSib.nc',
             two_d = False):
     """ Function extracting the topography from etopo 1/2" dataset
     """
     nc_f = Dataset(etopo, 'r')
     etopo_lat = nc_f.variables['lat'][:]
     etopo_lon = nc_f.variables['lon'][:]
     etopo_depths = nc_f.variables['Band1'][:]
     etopo_depths = np.ma.masked_where(etopo_depths>0, etopo_depths)
     
     if two_d == True:
         etopo_lat, etopo_lon = np.meshgrid(etopo_lat, etopo_lon)
         etopo_depths = etopo_depths.T
         
     return etopo_depths, etopo_lat, etopo_lon

#def plot_topo_section(etopo_depths, etopo_lat, etopo_lon, 
#                      zonal = True, lat_L = 78, 
#                      meridional = False, lon_L = 126):
#    """
#    Only one parameter should be chosen at time: zonal or meridional
#    """
#    if zonal==True:
#        cond_lat = np.logical_and(etopo_lat>= lat_L-0.01, etopo_lat<= lat_L+0.01)
#        cond_lon = np.logical_and(etopo_lon>=110, etopo_lat<= 170)
#
#        inds_78 = np.where(np.logical_and(cond_lat, cond_lon))
#        lons_vector = np.arange(74, 84, 0.01)
#        lats_vector = lat_L*np.ones(lons_vector.shape)
#
#        xy_to_interp = np.array([lons_vector, lats_vector]).T 
#        
#        
#    elif meridional == True:              
#        cond_lat = np.logical_and(etopo_lat>= 74, etopo_lat<= 85)
#        cond_lon = np.logical_and(etopo_lon>= lon_L-0.05, etopo_lat<= 84)
#
#        
#    etopo_depths[inds_78]
    

def set_ice_edge(m, reg_lat, reg_lon, bbox, fpath):
    """
    This function can be adupted for other cases, but was written for NABOS case with 2 particular datasets
    """
    fname_sicmin = '20180921120000-DMI-L4_GHRSST-STskin-DMI_OI-ARC_IST-v02.0-fv01.0.nc'
#    sst_sic, lat, lon, date, n_mes, masks = open_sst_L4.open_sst_l4(fpath, 
#                                                                    fname_sicmin, 
#                                                                    bbox, 
#                                                                    ice_mask_cond=True,
#                                                                    interp = True)
    sst_sic, lat, lon, date, n_mes, masks = open_sst_L4.open_sst_l4_interp(fpath, 
                                                                    fname_sicmin, 
                                                                    bbox, 
                                                                    ice_mask_cond=True)


    land_mask, ice_mask_interp_min = masks
    ice_mask_interp_min[np.logical_and(lat<80,lon< 145)]=0 #filter out small pieces of ice in weird places
    
       
#    fname_sicmed = '20180828120000-DMI-L4_GHRSST-STskin-DMI_OI-ARC_IST-v02.0-fv01.0.nc'
#    ice_mask_interp_med = open_sst_L4.get_seaiceedge(fpath, fname_sicmed, bbox)
    
    fname_sicmax = '20180821120000-DMI-L4_GHRSST-STskin-DMI_OI-ARC_IST-v02.0-fv01.0.nc'
#    sst_sic, lat, lon, date, n_mes, masks = open_sst_L4.open_sst_l4(fpath, 
#                                                                    fname_sicmax, 
#                                                                    bbox, 
#                                                                    ice_mask_cond=True,
#                                                                    interp = True)
    sst_sic, lat, lon, date, n_mes, masks = open_sst_L4.open_sst_l4_interp(fpath, 
                                                                    fname_sicmax, 
                                                                    bbox, 
                                                                    ice_mask_cond=True)
    land_mask, ice_mask_interp_max = masks
    ice_mask_interp_max[np.logical_and(lat<77,lon < 152)]=0 #filter out small pieces of ice in weird places
    ice_mask_interp_max[np.logical_and(lat<80.8,lon < 105)]=0 #filter out small pieces of ice in weird places

    x, y = m(reg_lon, reg_lat)

    m.contour(x, y, ice_mask_interp_max, 
                        colors = 'r', 
                        linewidths = 1, 
                        linestyles = 'dotted', 
                        corner_mask = False,
                        alpha= 0.5) 
    
#    m.contour(x, y, ice_mask_interp_med, 
#                        colors = 'w', 
#                        linewidths = 3, 
#                        linestyles = 'dotted', 
#                        corner_mask = False)
    
    m.contour(x, y, ice_mask_interp_min, 
                        colors = 'm', 
                        linewidths = 1, 
                        linestyles = 'dotted', 
                        corner_mask = False,
                        alpha= 0.5) 

    
    
    cntrmax = mlines.Line2D([], [], color='r', linestyle=':',
                            linewidth=3, label='Aug 21')
#    cntrmed = mlines.Line2D([], [], color='w', linestyle=':',
#                            linewidth=3, label='Aug 28')
    cntrmin = mlines.Line2D([], [], color='m', linestyle=':',
                            linewidth=3, label='Sept 21')
    
    legend = plt.legend([cntrmax, cntrmin], 
                        ['Aug 21', 'Sept 21'], loc = 'upper right')
#    legend = plt.legend([cntrmax,cntrmed, cntrmin], 
#                        ['Aug 5','Aug 28', 'Sept 21'])
    for label in legend.get_texts():
        label.set_fontsize('large')
    
    frame = legend.get_frame()
    frame.set_facecolor('0.90')
    return ice_mask_interp_max, ice_mask_interp_min



def set_map(resolution = 'l', data =[], reg_lat=[], reg_lon=[], 
            title = 'STD of SSS August-September 2018', size = (14, 8),
            figax = None, cmap = 'jet', clim_v =[]):
    """
    Function that make a base map for the NABOS expedition
    """
    bbox = [74, 85, 90, 170]
    latmin, latmax, lonmin, lonmax = bbox

    if figax:
        fig_sat, ax  = figax
    else:
        fig_sat = plt.figure(figsize = size)
        ax = fig_sat.add_axes([0.1,0.1,0.8,0.8])
        
    m = Basemap(projection='stere',lon_0=130, lat_0 = 78, 
                llcrnrlat=71,urcrnrlat=80,\
                llcrnrlon=110,urcrnrlon=190,
                resolution ='l')
    m.fillcontinents()
#    m.drawcoastlines()
    m.drawparallels(np.arange(70.,85, 2.),labels=[1,0,0,0],fontsize=14)
    
    m.drawmeridians(np.arange(0.,360.,5.),labels=[0,0,0,1],fontsize=14, latmax = 90)
    plt.title(title)

    if len(data)>0:
        im = m.pcolormesh(reg_lon, reg_lat, data, cmap = cmap, latlon = True)
        cbar = m.colorbar(im, "right", size="5%", pad='5%')
        cbar.ax.tick_params(labelsize=18)
        if len(clim_v)>0:
            plt.clim(vmin= clim_v[0], vmax = clim_v[1])
    return m, ax, fig_sat




def set_map_eastseas(resolution = 'l', data =[], 
                     reg_lat=reg_lat, reg_lon=reg_lon, title =[]):
    fig_sat = plt.figure(figsize=(16, 10))
    ax = fig_sat.add_axes([0.1,0.1,0.8,0.8])
    m = Basemap(projection='stere',
                lon_0=145, lat_0=70,
                llcrnrlat=60,urcrnrlat=78,
                llcrnrlon=95,urcrnrlon=240, 
                resolution ='l')
    
    m.fillcontinents()
    m.drawcoastlines()
    m.drawparallels(np.arange(50.,85, 2.),labels=[1,0,0, 0],fontsize=15)
    m.drawmeridians(np.arange(0.,360.,10.),labels=[0,0,0,1],fontsize=15, latmax = 90)    
    if len(title)>0: plt.title(title, fontsize=15)
    return m, ax, fig_sat



def set_map_arctic(resolution = 'l', title =[]):
    fig_sat = plt.figure(figsize = (14, 14))
    ax = fig_sat.add_axes([0.1,0.1,0.8,0.8])

    m = Basemap(projection='npstere', boundinglat=60, lon_0=140, 
                resolution = resolution)
    
    m.drawcoastlines()
    m.fillcontinents(color='white',lake_color='white')
    m.drawparallels(np.arange(0, 85, 5), labels=[True, False, False, False]) # labels = [left,right,top,bottom]
    m.drawmeridians(np.arange(-180, 180, 10),labels=[False, False, False, True])
    m.drawmapboundary(fill_color='white')
    return m, ax, fig_sat

#
#struct2 = ndimage.generate_binary_structure(2, 2)
#erode2 = ndimage.binary_erosion(sss, struct2)
#edges2 = sss_ma.mask ^ erode2
#plt.contour(reg_lon, reg_lat, edges2)
