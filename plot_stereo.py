#!/usr/bin/env python
# Plotting Sentinel-1 data in stereo projections:
# Arctic - whole Arctic region
# region - depending on the chosen case plot over the whole area (Greenland, Kara, Kara strait, ...)
# zoom - zoom one particular image according to its boundaries, but keep it in polar stereo projection

import numpy as np
import os, sys
import datetime
from geopy.distance import vincenty


sys.path.append('/home/cercache/users/atarasen/codes/python/')
import nice_functions as nf
#import ncdfView as nv

sys.path.append('/home3/homedir7/perso/atarasen/codes/python/collocation')
import buoy as buoy_module

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.basemap import Basemap, maskoceans
from matplotlib.path import Path

import scipy.io
sys.path.append('/home3/homedir7/perso/atarasen/codes/python/posada/')
sys.path.append('/home3/homedir7/perso/atarasen/codes/python/cmod/')

from posada.Handlers.sar import SAR
from posada.Models.swath import Swath
from posada.Readers.safe import SAFEGeoTiffFile

#     Set up colorbar settings (for all images)

path_out_def='/home3/homedir7/perso/atarasen/Figures/'
col_min={'hh':0, 'hv':0, 'hh_log':-30, 'hv_log': -40, 'hh_sobel':0, 'hv_sobel':0}
col_max={'hh':0.1, 'hv':0.01, 'hh_log':0, 'hv_log': -15,'hh_sobel':50, 'hv_sobel':30}
col_bar_label={'hh': r'sigma_0 [dN]', 'hv': r'sigma_0 [dN]', 'hh_log': r'sigma_0 [dB]', \
                   'hv_log': r'sigma_0 [dB]','hh_sobel':'','hv_sobel':''}
col_map_var={'hh':'bone', 'hv':'bone', 'hh_log':'bone', 'hv_log':'bone', 'hh_sobel':'gist_heat', 'hv_sobel':'gist_heat'}
col_buoy={'hh':'r-', 'hv':'r-', 'hh_log':'r-', 'hv_log':'r-', 'hh_sobel':'b-', 'hv_sobel':'b-'}
#col_buoy={'hh':'r*', 'hv':'r*', 'hh_log':'r*', 'hv_log':'r*', 'hh_sobel':'b*'}


width_dic={'Greenland': 3000000, 'Fram strait': 1200000,  'Kara': 1800000, 'Kara strait': 600000}
height_dic={'Greenland': 3000000, 'Fram strait': 1200000, 'Kara': 1700000, 'Kara strait': 600000}
lat_ts_dic= {'Greenland': 70, 'Fram strait': 78, 'Kara': 70, 'Kara strait': 70}
lat_0_dic = {'Greenland': 75, 'Fram strait': 78, 'Kara': 77, 'Kara strait': 71.5}
lon_0_dic = {'Greenland': -10, 'Fram strait': -1, 'Kara': 73, 'Kara strait': 59.}


def latlon_minmax(lat_ice, lon_ice):
    lat_min=np.nanmin([lat_ice[0, 0], lat_ice[0, -1], lat_ice[-1, -1], lat_ice[-1, 0]])
    lat_max=np.nanmax([lat_ice[0, 0], lat_ice[0, -1], lat_ice[-1, -1], lat_ice[-1, 0]])
    lon_min=np.nanmin([lon_ice[0, 0], lon_ice[0, -1], lon_ice[-1, -1], lon_ice[-1, 0]])
    lon_max=np.nanmax([lon_ice[0, 0], lon_ice[0, -1], lon_ice[-1, -1], lon_ice[-1, 0]])
    
    #nest expressions check for nans:
    if lat_min!=lat_min: 
        lat_min=np.nanmin(lat_ice.ravel())
        lat_max=np.nanmax(lat_ice.ravel())
        
    if lon_min!=lon_min:
        lon_min=np.nanmin(lon_ice.ravel())
        lon_max=np.nanmax(lon_ice.ravel())
    
    #lon_min=np.min([lon_ice[0, 0], lon_ice[-1, -1] ])
    #lon_max=np.max([lon_ice[0, 0], lon_ice[-1, -1] ])
    return lat_min, lat_max, lon_min, lon_max



# Trying to solve the eternal problem of passing 0 meridian:
def s1_lon360(lon_s1):
    l_360=[]
    for l in lon_s1.flatten(order='C'): #C means to flatten in row-major (C-style) order
        if l<0:
            l_360.append(360+l)
        else:
            l_360.append(l)
    lon_360=np.reshape(l_360, lon_s1.shape, order='C')    
    return lon_360

def lon_center_360(lon_max, lon_min):
    if lon_max>300 and lon_min<60:
        print 'image crossing 0 meridian'
        lon_center= np.mean((lon_max-360, lon_min),)
        if lon_center<0:
            lon_center=lon_center+360  
    else:
        lon_center=np.mean((lon_max,lon_min),)        
    return lon_center

def lon_center_180(lon_max, lon_min):
#    if both values have different sign around 180 meridian
    if (lon_min<-150) and (lon_max>150):
        print 'image crossing 180 meridian'
        l_diff_min=180+lon_min
        l_diff_max=180-lon_max
        if l_diff_min>l_diff_max:
            lon_center = (l_diff_min + l_diff_max)/2-180
        else:
            lon_center = 180-(l_diff_min + l_diff_max)/2
    else:
        lon_center=np.mean((lon_max,lon_min),)
    return lon_center
   

# Import buoy positions and select (either mean daily or all points)    
def buoy_traj(S1_name, case, buoy_case, daily=False):
    day_x = datetime.datetime(int(S1_name[17:21]), int(S1_name[21:23]), int(S1_name[23:25]), 12)

    fpath='/home/cercache/users/atarasen/data/Buoys/Arctic-buoy-data/'+case+'/' 
    years, days, lons, lats = buoy_module.importtext(fpath, buoy_case+'.dat')
    if daily==True:
        lats, lons, buoy_days = buoy_module.buoy_daily(years, days, lons, lats)
    else:
        buoy_days=[]
        for ii, day in enumerate(days):
            buoy_days.append(buoy_module.yearday2struct_time(years[ii], day))
            
    ind=np.where(np.array(buoy_days)<=day_x)
    
    buoy_lats=np.array(lats)[ind] 
    buoy_lons=np.array(lons)[ind]
    
    return buoy_lats, buoy_lons, np.array(buoy_days)[ind]
   
def buoy_traj_plot(buoy_lats, buoy_lons, m):
    x_b, y_b = m(buoy_lons, buoy_lats)
    m.plot(x_b, y_b, col_buoy[pol])    

    
#real plotting
def simple(data, pol=[]):
    col_map = plt.get_cmap(col_map_var[pol])
    im = plt.pcolormesh(data, cmap=col_map)
    plt.show()
    #plt.clim(col_min[pol],col_max[pol])
    

    #cbar = plt.colorbar(im, cmap=col_map)
    #cbar.set_label(col_bar_label[pol])



def Arctic(ice_data, lat_ice, lon_ice, S1_name, pol, case, buoy_case=None, path_out=path_out_def, fig_title='', col_bar=True):
    fig_name=path_out+S1_name+'.png'

    proj='npstere'
    lat_min = -80.;  lat_max=81.;  lat_period=20.
    lon_min = -180.; lon_max=181.; lon_period=20.
    m = Basemap(projection=proj,boundinglat=60,lon_0=0,resolution='l')

    x, y = m(lon_ice, lat_ice)
    m.drawcoastlines()
    m.fillcontinents(color='white',lake_color='white')
    m.drawparallels(np.arange(0, 85, 5), labels=[True, False, False, False]) # labels = [left,right,top,bottom]
    m.drawmeridians(np.arange(-180, 180, 10),labels=[False, False, False, True])
    m.drawmapboundary(fill_color='white')
    
    col_map = plt.get_cmap(col_map_var[pol])
    im = m.pcolormesh(x, y, ice_data, cmap=col_map)

    plt.clim(col_min[pol],col_max[pol])
    if col_bar==True:
        cbar = m.colorbar(im,location='right',pad="5%", cmap=col_map)
        cbar.set_label(col_bar_label[pol])
#   open buoy trajectory and plot it over: 
    if buoy_case!=None:
        buoy_lats, buoy_lons, buoy_days = buoy_traj(S1_name, case, buoy_case, daily=False)
        x_b, y_b = m(buoy_lons, buoy_lats)
        m.plot(x_b, y_b,col_buoy[pol])    

    if fig_title=='' and buoy_case!=None:
        fig_title=S1_name[:32]+"..."+S1_name[-4:]+"\n pol = "+pol+"\n Buoy "+buoy_case
    elif fig_title=='':
        fig_title = S1_name[:32]+"..."+S1_name[-4:]+"\n pol = "+pol
    plt.title(fig_title)
        
    
    
def region(ice_data, lat_ice, lon_ice, S1_name, pol, case, buoy_case, area_case='', path_out=path_out_def, fig_title='', col_bar=True, mask_coord=[]):    
# lat_ts is latitude of true scale; lon_0,lat_0 is central point.        
    if area_case=='':
        area_case=np.copy(case)
    m = Basemap(width=width_dic[area_case],height=height_dic[area_case],resolution='l',projection='stere',\
        lat_ts=lat_ts_dic[area_case], lat_0=lat_0_dic[area_case], lon_0=lon_0_dic[area_case])

    x, y = m(lon_ice, lat_ice)
    print x.shape
    
    
    m.drawcoastlines()
    m.fillcontinents(color='white',lake_color='white')    
    m.drawparallels(np.arange(0, 85, 5), labels=[True, False, False, False]) # labels = [left,right,top,bottom]
    m.drawmeridians(np.arange(-180, 180, 10),labels=[False, False, False, True])
    m.drawmapboundary(fill_color='white')
    col_map = plt.get_cmap(col_map_var[pol])

#    data_mask=np.zeros([x.shape[0], y.shape[1]])
#    if mask_coord!=[]:
#       coord_lon=mask_coord[1]; coord_lat=mask_coord[0]
#        print coord_lon
        
#        m_poly_x, m_poly_y = m(sum(coord_lon, []), sum(coord_lat, []))
        
#        intersection_vert = np.zeros([len(m_poly_x), 2])
#        for i_p, in range(len(m_poly_x)):
#            intersection_vert[i_p, 0] = m_poly_x[i_p]
#            intersection_vert[i_p, 1] = m_poly_y[i_p]
                                     
#        m_polygon = Path(intersection_vert)
#        print m_polygon
#        good_points=[]
#         for (x_p, y_p) in [x, y]:
               #if m_polygon.contains_point(x_p, y_p):
    #               good_points.append([x_p, y_p]) 


#        m.scatter(m_poly_x, m_poly_y, cmap='afmhot')
    
#     im = m.pcolormesh(x[ind], y[ind], ice_data[ind], cmap=col_map)
    im = m.pcolormesh(x, y, ice_data, cmap=col_map)

    


    plt.clim(col_min[pol],col_max[pol])
    if col_bar==True:
        cbar = m.colorbar(im,location='right',pad="5%", cmap=col_map)
        cbar.set_label(col_bar_label[pol])

#   open buoy trajectory and plot it over: 
    if buoy_case!=None:    
        buoy_lats, buoy_lons, buoy_days = buoy_traj(S1_name, case, buoy_case, daily=True)
        buoy_traj_plot(buoy_lats, buoy_lons, m)

    if fig_title=='' and buoy_case!=None:
        fig_title=S1_name[:32]+"..."+S1_name[-4:]+"\n pol = "+pol+"\n Buoy "+buoy_case
    elif fig_title=='':
        fig_title = S1_name[:32]+"..."+S1_name[-4:]+"\n pol = "+pol
    plt.title(fig_title)
        

        
def zoom(ice_data, lat_ice, lon_ice, S1_name,\
         case, buoy_case=None, path_out=path_out_def, fig_title='',\
         pol= '', col_bar=True, col_map_name='', clim=[], x=[], y=[], x_b=[], y_b=[]):
# Convert lats into 360:
#    lon_ice=s1_lon360(lon_ice)
     
    lat_min, lat_max, lon_min, lon_max = latlon_minmax(lat_ice, lon_ice)
    lon_center = lon_center_360(lon_max,lon_min)
    lat_center = np.mean((lat_max,lat_min),)
    print lat_min, lat_max, lon_min, lon_max

    height_lats = int(vincenty((lat_min, lon_min),(lat_max, lon_min)).meters)
    width_lons = int(vincenty((lat_min, lon_min),(lat_min, lon_max)).meters)

    #print height_lats,  width_lons       
    #m = Basemap(width=600000,height=600000,
    
    m = Basemap(width=width_lons,height=height_lats,
    resolution='h',projection='stere',\
    lat_ts=lat_center,lat_0=lat_center,lon_0=lon_center)
            
    x, y = m(lon_ice, lat_ice)
    m.drawcoastlines()
    #m.fillcontinents(color='white',lake_color='white')
    m.drawparallels(np.arange(0, 85, 5), labels=[True, False, False, False])
    m.drawmeridians(np.arange(-180, 180, 5), labels=[False, False, False, True])    
    m.drawmapboundary(fill_color='white')
    
    
    if col_map_name=='':
        col_map = plt.get_cmap(col_map_var[pol])
    else:
        col_map = plt.get_cmap(col_map_name)
    im = m.pcolormesh(x, y, ice_data, cmap=col_map)
    
    
    if clim==[] and pol!='': 
        plt.clim(col_min[pol],col_max[pol])       
    elif clim!=[]:
        plt.clim(clim[0],clim[1])

        
    
    if pol!='':
        if col_bar==True:
            cbar = m.colorbar(im,location='right',pad="5%", cmap=col_map)
            cbar.set_label(col_bar_label[pol])
    else:
        if col_bar==True:
            cbar = m.colorbar(im,location='right',pad="5%")
    
       
#   open buoy trajectory and plot it over: 
    if buoy_case!=None:
        if x_b==[]:
            buoy_lats, buoy_lons, buoy_days = buoy_traj(S1_name, case, buoy_case, daily=False)
            buoy_traj_plot(buoy_lats, buoy_lons, m)
 

    if fig_title=='' and buoy_case!=None:
        fig_title=S1_name[:32]+"..."+S1_name[-4:]+"\n pol = "+pol+"\n Buoy "+buoy_case
        fig_name=path_out+buoy_case+'_'+S1_name+'_zoom.png'

    elif fig_title=='':
        fig_title=S1[:32]+"..."+S1_name[-4:]+"\n pol = "+pol
        fig_name=path_out+S1_name+'_zoom.png'

    plt.title(fig_title)
    
    return x, y, x_b, y_b, lat_center, lon_center, height_lats, width_lons



def zoom_area(data, lat, lon, S1_name,\
              case, width_lons, height_lats, lat_center, lon_center,\
              path_out=path_out_def, fig_title='',\
              pol= '', col_bar=True, col_map_name='', clim=[]):

    
    m = Basemap(width=width_lons,height=height_lats,
    resolution='h',projection='stere',\
    lat_ts=lat_center,lat_0=lat_center,lon_0=lon_center)
            
    x, y = m(lon, lat)
    m.drawcoastlines()
    m.fillcontinents(color='white',lake_color='white')
    m.drawparallels(np.arange(0, 85, 5), labels=[True, False, False, False])
    m.drawmeridians(np.arange(-180, 180, 5), labels=[False, False, False, True])    
    m.drawmapboundary(fill_color='white')
    
    
    if col_map_name=='':
        col_map = plt.get_cmap(col_map_var[pol])
    else:
        col_map = plt.get_cmap(col_map_name)
    im = m.pcolormesh(x, y, data, cmap=col_map)
    
    
    if clim==[] and pol!='': 
        plt.clim(col_min[pol],col_max[pol])       
    elif clim!=[]:
        plt.clim(clim[0],clim[1])
       
    
    if pol!='':
        if col_bar==True:
            cbar = m.colorbar(im,location='right',pad="5%", cmap=col_map)
            cbar.set_label(col_bar_label[pol])
    else:
        if col_bar==True:
            cbar = m.colorbar(im,location='right',pad="5%")
    

    if fig_title=='':
        fig_title=S1[:32]+"..."+S1_name[-4:]+"\n pol = "+pol
        fig_name=path_out+S1_name+'_zoom.png'

    plt.title(fig_title)
    
    return lat_center, lon_center, height_lats, width_lons

    
    
    

#     plt.savefig(fig_name, dpi=300)
#     plt.show()    
#     nf.pause()
          
    
    
#    ax = plt.gca()
#    divider = make_axes_locatable(ax)
#    cax = divider.append_axes("right", size="5%", pad=0.05)
#     plt.colorbar(im, cax=cax, cmap=plt.cm.inferno)