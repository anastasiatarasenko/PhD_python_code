### Set of basic python functions to extract the wind speed and direction from Sentinel-1 image
# Should be adapted to use without SOLAB servers, as for wind directions posada needs reanalysis-data from ncep or era-interim

# author: Tarasenko Anastasia with posada code (Solab), Jul 2017 


import numpy as np

import matplotlib
import matplotlib.pyplot as plt

import scipy.signal as sg
import scipy.ndimage as ndimage

import sys, os

# sys.path.append('/home/mag/Documents/repos/solab/POSADA/') #common version on student server
## sys.path.append('/home/tarasenko/codes/posada3/') #was old version (before switching to posada3 on python3)
sys.path.append('/home/tarasenko/prj/posada3/') #last version on my computer

from posada.Handlers import Sentinel_1 as SAR
from posada.Models import Swath
from posada.Readers import SAFE


sys.path.append('/home/tarasenko/codes/python/visu_stuff/')
import nice_visu as visu



def s1_trj(lat, lon, ws, wdir_model):
#Model wind direction is given in mathemathical coordinate system (cos - righward, sin - ipward, anticlockwise). At the same time,  model wind direction ncep and interim) coordinate systems oriented to the north (sin is going northwards). This creates problems when the trajectory is descending
    
    if lat[0,0]-lat[-1,-1]>0:
        trj = 'desc'
        ws = np.flipud(np.fliplr(ws))
        wdir_model = np.flipud(np.fliplr(wdir_model-180))
        if lat != [] and lon != []:
            lat = np.flipud(np.fliplr(lat))
            lon = np.flipud(np.fliplr(lon))
    else:
        trj = 'asc'
    print (trj +' trajectory') 
    return trj, ws, wdir_model, lat, lon


def wind_s1(pn, fn, res_m = 500, pol = 'hh'):
# Get data from Sentinel1 image
# sar.get_field_names() 
    sarf = SAFE(url=pn+fn, polarisation=pol)
    sar = Swath()
    sar.load(sarf)
    sar.add_handler(SAR)
    sar.set_spacing(sar.resolution2spacing(res_m))

    s0 = sar.get_values('sigma0')
    lat = sar.get_values( 'lat')
    lon = sar.get_values( 'lon')
    
   
    # Nice images from Albert Murtazin
    s0_corr = sar.get_values('nice_image')

    # Wind values
    wind_speed = sar.get_values('wind_speed')
    
    # Wind directions from NCEP
    wdir_ncep = sar.get_values('wind_dir_ncep_gfs')-180
    
    # Wind directions from ERA Interim
    wdir_interim = sar.get_values('wind_dir_era_interim')-180
    
    
    if lat[0,0]-lat[-1,-1]>0:
        trj = 'desc'
        
        #rotate everything to get good vectors for quiver
        s0 = np.flipud(np.fliplr(s0))
        s0_corr = np.flipud(np.fliplr(s0_corr))

        wind_speed = np.flipud(np.fliplr(wind_speed))
        wdir_ncep = np.flipud(np.fliplr(wdir_ncep))
        wdir_interim = np.flipud(np.fliplr(wdir_interim))

        lat = np.flipud(np.fliplr(lat))
        lon = np.flipud(np.fliplr(lon))
    else:
        trj = 'asc'
    print (trj +' trajectory') 
    
    return s0, s0_corr, lat, lon, wind_speed, wdir_ncep, wdir_interim, trj





def wind_s1_thermn_corr(pn, fn, res_m = 500, pol = 'hh'):
# Get data from Sentinel1 image
# sar.get_field_names() 
    sarf = SAFE(url=pn+fn, polarisation=pol)
    sar = Swath()
    sar.load(sarf)
    sar.add_handler(SAR)
    sar.set_spacing(sar.resolution2spacing(res_m))

    s0 = sar.get_values('sigma0')
    lat = sar.get_values( 'lat')
    lon = sar.get_values( 'lon')
    

    s0_therm_noise_corr = sar._get_values_nice_image_hh_vv
    
    
    # Nice images from Albert Murtazin
    s0_corr = sar.get_values('nice_image')

    # Wind values
    wind_speed = sar.get_values('wind_speed')
    
    # Wind directions from NCEP
    wdir_ncep = sar.get_values('wind_dir_ncep_gfs')-180
    
    # Wind directions from ERA Interim
    wdir_interim = sar.get_values('wind_dir_era_interim')-180
    
    
    if lat[0,0]-lat[-1,-1]>0:
        trj = 'desc'
        
        #rotate everything to get good vectors for quiver
        s0 = np.flipud(np.fliplr(s0))
        s0_corr = np.flipud(np.fliplr(s0_corr))

        wind_speed = np.flipud(np.fliplr(wind_speed))
        wdir_ncep = np.flipud(np.fliplr(wdir_ncep))
        wdir_interim = np.flipud(np.fliplr(wdir_interim))

        lat = np.flipud(np.fliplr(lat))
        lon = np.flipud(np.fliplr(lon))
    else:
        trj = 'asc'
    print (trj +' trajectory') 
    
    return s0, s0_corr, lat, lon, wind_speed, wdir_ncep, wdir_interim, trj





    
    
    
    
def uv_wdir(wind_speed, wdir):
    #wdir is given in geographical gedrees (0 at North, 80 at East, etc)
    u_wind_1d = wind_speed.squeeze() * np.cos(np.deg2rad(wdir).squeeze())
    v_wind_1d = wind_speed.squeeze() * np.sin(np.deg2rad(wdir).squeeze())
    
    u_wind = u_wind_1d.reshape(wind_speed.shape)
    v_wind = v_wind_1d.reshape(wind_speed.shape)

    return u_wind, v_wind



def wind_mean(wind_speed, wdir):
# Mean wind speed and direction:

    ws_mean = np.mean(wind_speed.squeeze()[wind_speed<=18])
    wdir_mean = np.mean(wdir.squeeze())

    u_mean, v_mean = uv_wdir(ws_mean, wdir_mean)
    
    print('Mean wind direction in deg (from North):')
    print(wdir_mean)
    
    print('Mean wind speed')
    print(ws_mean)
    
    return ws_mean, wdir_mean, u_mean, v_mean



def section_extract(data, x0, y0, wdir_mean, distance, resolution):
#     wdir should be in degrees (as in reanalysis products, distance in meters)

    delta_x = distance*np.cos(np.deg2rad(wdir_mean))/resolution
    delta_y = distance*np.sin(np.deg2rad(wdir_mean))/resolution
    
    #print('delta_x, delta_y', delta_x, delta_y)
    x1, y1 = x0+delta_x, y0+delta_y
    section_coords = [x0, x1, y0, y1]

    if x1>data.shape[0] or y1 > data.shape[1]:
        print("Section length is beyond the image size")
        data_section, distance_vector = [], []    
    else:
        distance_vector = np.arange(0, distance, resolution)
        num = np.round(distance/resolution)
        x, y = np.linspace(x0, x1, num), np.linspace(y0, y1, num)

# Extract the values along the line, using cubic interpolation
        data_section = ndimage.map_coordinates(data, np.vstack((x,y)))
    return data_section, distance_vector, section_coords




def wind_plot_all(fn, s0, s0_corr, ws, wdir_ncep, wdir_interim, lon, lat, every_x = 50):   
    # Plot it   
    
    fig_4basics = plt.figure(figsize=(20, 18))
    plt.subplot(2, 2, 1)
    plt.pcolormesh(10*np.log10(s0), cmap='jet')
    plt.colorbar()
    plt.clim(-40, -5)
    plt.title(fn)

    
    plt.subplot(2, 2, 2)
    plt.pcolormesh(ws, cmap='jet')
    plt.clim(0, 25)
    plt.colorbar()
    
    u_ncep, v_ncep = uv_wdir(ws, wdir_ncep)
    u_interim, v_interim = uv_wdir(ws, wdir_interim)

    
    coords_0, coords_1 = np.meshgrid(np.arange(ws.shape[1]), np.arange(ws.shape[0]))
        ###x, y = np.meshgrid(np.arange(ws.shape[0]), np.arange(ws.shape[1]))


        
    plt.quiver(coords_0[::every_x, ::every_x], coords_1[::every_x, ::every_x],\
               v_ncep[::every_x, ::every_x], u_ncep[::every_x, ::every_x],\
               scale=30, scale_units='inches', color = 'white')
    
    plt.quiver(coords_0[::every_x, ::every_x], coords_1[::every_x, ::every_x],\
               v_interim[::every_x, ::every_x], u_interim[::every_x, ::every_x],\
               scale=30, scale_units='inches', color = 'red')    

    plt.title('wind speed from S1 with CMOD5\n wind dir: NCEP (white) / ERA-Interim (red)')


#nice images from Albi
    plt.subplot(2, 2, 3)
    plt.pcolormesh(s0_corr)
    plt.clim(-5, 5)
    plt.colorbar()
    plt.title('Nice image')


#Difference between the models
    plt.subplot(2, 2, 4)
    plt.pcolormesh(wdir_ncep-wdir_interim, cmap='jet');
    plt.clim([-40, 40])
    plt.colorbar()
    plt.title('Difference in wind directions [ncep - interim]')

    plt.show()
    return fig_4basics



def wind_plot(fn, ws, wdir, lat, lon, every_x = 50, case = 'ERA INTERIM', plot_case = [],  section_coords=[], sect_cols=[]):
    #fig_wind = plt.figure(figsize=(12, 12))

    
    if plot_case == 'geo':
        plt.pcolormesh(lon, lat, ws, cmap='jet',)
        plt.xlabel('longitude')
        plt.ylabel('latitude')
    else:
        plt.pcolormesh(ws, cmap='jet',)

    cbar = plt.colorbar()
    cbar.set_label('m/s')
    plt.clim(0, 25)
    
    u, v = uv_wdir(ws, wdir)
    _, _, u_mean, v_mean = wind_mean(ws, wdir)
    
    mean_p_0 = int(ws.shape[0]/2)
    mean_p_1 = int(ws.shape[1]/2)
    
    
    if plot_case == 'geo':
        coords_0, coords_1 = lon, lat
    else:
        coords_v, coords_h = np.meshgrid(np.arange(ws.shape[1]), np.arange(ws.shape[0]))
        ###x, y = np.meshgrid(np.arange(ws.shape[0]), np.arange(ws.shape[1]))
        #plt.xticks(range(0, 850, 200), np.arange(0, 450, 100))
        #plt.yticks(range(0, 850, 200), np.arange(0, 450, 100))
        #plt.xlabel('distance, km')


        
    plt.quiver(coords_v[::every_x, ::every_x], coords_h[::every_x, ::every_x],\
               v[::every_x, ::every_x], u[::every_x, ::every_x],\
               scale=30, scale_units='inches', color = 'red')

    plt.quiver(coords_v[mean_p_0, mean_p_1], coords_h[mean_p_0, mean_p_1],\
               v_mean, u_mean,\
               scale=30, scale_units='inches', color ='black')  
       

    plt.title('Wind speed from\n'+fn+'\n wind direction from '+case)
    plt.tight_layout()
    #plt.show()


# Section on the map
    if section_coords!=[]:
        if type(sect_cols) != str:               
            if sect_cols==[]:
                sorted_names = visu.color_list()
                sect_cols = sorted_names[::4]
            
            for n in np.arange(len(section_coords['x0'])):   
                plt.plot([section_coords['y0'][n], section_coords['y1'][n]], \
                         [section_coords['x0'][n], section_coords['x1'][n]], \
                         color = sect_cols[n], marker = 'o', linewidth = 4)
        else:
            for n in np.arange(len(section_coords['x0'])):   
                plt.plot([section_coords['y0'][n], section_coords['y1'][n]], \
                             [section_coords['x0'][n], section_coords['x1'][n]], \
                             color = sect_cols, marker = 'o', linewidth = 4)
    #return fig_wind

    
def section_plot_1d(distance_vector, data_section):
    distance_km = distance_vector/1000
    plt.plot(distance_km, data_section)
    plt.xlabel('distance, km')
    plt.ylabel('wind speed, m/s')



