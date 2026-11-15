#! /usr/bin/env python


#==================
# Script to download Sentinel data (1, 2 or 3) in the region of interest defined by buoys positions (using buoy_date_region module) and sea ice edge from MW products (using ice_open module)

# authors: Anastasia Tarasenko (& peps team for S1 downloading)
# creation: oct 2016
#==================
import numpy as np
import numpy.ma as ma
import os, sys
from operator import add
import copy
from scipy import interpolate

import shapely
from shapely.geometry import Point, Polygon, MultiPolygon
from shapely.ops import cascaded_union, unary_union
from descartes.patch import PolygonPatch

import pickle

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.basemap import Basemap

import datetime

# sys.path.append('/home/cercache/users/atarasen/codes/python/collocation')
sys.path.append('/home/cercache/users/atarasen/codes/python/')
sys.path.append('/home3/homedir7/perso/atarasen/codes/python/collocation')

import nice_functions as nf

import buoy as buoy_module
reload(buoy_module)

import ice as ice_module

import S1_open_process_and_visu as S1_open
reload(S1_open)

#import ncdfView as nv

def date_list(t1, t2):
    numdays = (t2-t1).days
    dt_list = [t1 + datetime.timedelta(days=1) for x in range(0, numdays)]
    return dt_list	



def interpolate_ice(lon_ice, lat_ice, lat_x, lon_x, area_lon=0.5, area_lat=0.5):
    ind_area0, ind_area1 = np.where(np.logical_and(np.abs(lon_ice-lon_x)<area_lon, np.abs(lat_ice-lat_x)<area_lat))

    # Enlarge the area of search, if no data was found  
    area_n=area_lat*2
    while (not ind_area0.any()):
        ind_area0, ind_area1 = np.where(np.logical_and(np.abs(lon_ice-lon_x)<(area_n), np.abs(lat_ice-lat_x)<(area_lat*area_n)))
        area_n=area_n+area_lat
        
    
    l=[]
    for i_count, ind0 in enumerate(ind_area0):
        ind1 = ind_area1[i_count]
          
        x=lon_ice[ind0, ind1]-lon_x
        y=lat_ice[ind0, ind1]-lat_x
        l.append(np.sqrt(x**2 +y**2))

    ind_ind = l.index(np.min(l))  
    true_ind = [ind_area0[ind_ind], ind_area1[ind_ind]] 
    return true_ind    


#Define if the files are sequential based on their names
def mosaic_condition(s1_name1, s1_name2, dt=5):
    d_1 = S1_open.s1nametodatetime(s1_name1)
    d_2 = S1_open.s1nametodatetime(s1_name2)
    delta_t = datetime.timedelta(hours=(dt/60.)) #dt is given in minutes
    
    #if abs(d_2-d_1)<delta_t and d_2!=d_1:
    if abs(d_2-d_1)<delta_t:
        mosaic_ans=True
    else:
        mosaic_ans=False
    return mosaic_ans



def full_latlon(lats, lons, step=0.1):
##    function that recreates the full 2D lat-lon matrix based on 4 corner values(p1(lat, lon), p2(lat, lon)...)
    b_lat = sorted([np.round(lats[0], 1), np.round(lats[3], 1)])
    b_lon = sorted([np.round(lons[0], 1), np.round(lons[1], 1)])
    
    lats_step_base=np.ndarray.tolist(np.arange(b_lat[0], b_lat[1], step))
    lons_step_base=np.ndarray.tolist(np.arange(b_lon[0], b_lon[1], step))    
    
# create 2D mesdgrid for the lat[0]:lat[3] points, lon[0]:0.1:lon[1]
    new_lat, new_lon = np.meshgrid(lats_step_base, lons_step_base)   
    
# make last row: lat[1]:n-of-points-in-1st-row:lat[2]
    len_i = (np.array(new_lat)).shape[0];
    d_i = (lats[1]-lats[0])/len_i    
    len_j = (np.array(new_lat)).shape[1]
    d_j = (lons[3]-lons[0])/len_j
    
    full_lat=copy.deepcopy(new_lat); full_lon = copy.deepcopy(new_lon)
    for i in range(len_i):
        for j in range(len_j):
            full_lat[i,j]=np.round(new_lat[i, j]+d_i*i, 1)   
            full_lon[i,j]=np.round(new_lon[i, j]+d_j*j, 1)                         
    full_lon =  np.rot90(full_lon, 2)  
    #print 'shape lat', full_lat.shape
    #print 'shape_lon', full_lon.shape
    return full_lat, full_lon



### INTERSECTION FUNCTIONS
# Additional small functions helping to plot the overlapped section
# function 1
def v_color(ob):
    COLOR = {
    True:  '#6699cc',
    False: '#ff3333'
    }
    return COLOR[ob.is_valid]

#Function 2
def plot_coords(ax, ob, colors='r'):
    x, y = ob.xy
    ax.plot(x, y, 'o', color=colors, zorder=1)
    #ax.plot(x, y, 'o', color='#999999', zorder=1)
    
    
# Main function that check overlapping
def intersection_condition(corners1, corners2, plot_cond=False):
    polygons_to_compare=[]  

    for i_c, corners in enumerate([corners1, corners2]):
# If the provided collection of 4-corners-points coordinates consists of description of several files 
# (thus the form of the list is [ [ [.., ..], [], [], [] ], [....] ]:    
        if isinstance(corners[0][0], list)==True:
            polygons_full = []
            for ppp in range(len(corners)):
                polygons_full.append(Polygon(np.array(corners)[ppp]))
            merged=cascaded_union(polygons_full)
        else:
            print 'one image'
            merged=Polygon(np.array(corners))            
        polygons_to_compare.append(merged)
       
#        if plot_cond==True:
##             if merged.type=='MultiPolygon':
#                 for polygon in merged:
#                     plot_coords(ax, polygon.exterior, col[i_c])
#                     patch_small = PolygonPatch(polygon, facecolor=v_color(polygon), edgecolor=col[i_c], alpha=0.5, zorder=2)
#                     ax.add_patch(patch_small)
#             else:
#                 plot_coords(ax, merged.exterior, col[i_c])
#                 patch_small = PolygonPatch(merged, facecolor=v_color(merged), edgecolor=col[i_c], alpha=0.5, zorder=2)
#                 ax.add_patch(patch_small)
                           
    c = polygons_to_compare[0].intersection(polygons_to_compare[1])
    print 'area of intersection', c.area
    critique_area = 0.5
    cond_area = (c.area > critique_area)
    print 'cond_area', cond_area
    
    if c.is_empty==False and cond_area==True:
        ans = True
        print 'Images overlap'
        
        if plot_cond==True:
            plt.close()
            fig = plt.figure(1, figsize=(8, 8), dpi=90)
            ax = fig.add_subplot(111)
            col = ['k', 'b', 'r', 'm', 'w', 'y', 'g', 'c']
            
            for merged in polygons_to_compare:
                if merged.type=='MultiPolygon':
                    for polygon in merged:
                        plot_coords(ax, polygon.exterior, col[i_c])
                        patch_small = PolygonPatch(polygon, facecolor=v_color(polygon), edgecolor=col[i_c], alpha=0.5, zorder=2)
                        ax.add_patch(patch_small)
                else:
                    plot_coords(ax, merged.exterior, col[i_c])
                    patch_small = PolygonPatch(merged, facecolor=v_color(merged), edgecolor=col[i_c], alpha=0.5, zorder=2)
                    ax.add_patch(patch_small)        
                
        
        coord=[]; coord_lat=[]; coord_lon=[]
        if c.type=='MultiPolygon':
            for i_color, polygon_c in enumerate(c):
                coord.append(list(polygon_c.exterior.coords))
                
                #coord_x_temp, coord_y_temp = polygon_c.exterior.xy
                #coord_lat.append(coord_x_temp.tolist())
                #coord_lon.append(coord_y_temp.tolist())
                if plot_cond==True:
                        plot_coords(ax, polygon_c.exterior)
                        patch = PolygonPatch(polygon_c, facecolor=col[i_color], edgecolor=v_color(c), alpha=0.5, zorder=2)
                        ax.add_patch(patch)
            #coord=[coord_lat, coord_lon]

        else:
            coord =list(c.exterior.coords)
            coord.append(list(c.exterior.coords))
            #coord_x_temp, coord_y_temp = c.exterior.xy
            #coord_lat.append(coord_x_temp.tolist())
            #coord_lon.append(coord_y_temp.tolist()) 
            #coord=[coord_lat, coord_lon]

            if plot_cond==True:
                    plot_coords(ax, c.exterior)
                    patch = PolygonPatch(c, facecolor='r', edgecolor=v_color(c), alpha=0.5, zorder=2)
                    ax.add_patch(patch)   

    else:
        coord =[]
        ans = False       
        
    if plot_cond==True:    
        plt.show()
    
    return ans, coord, c


# Function that extracts the points inside overlapped area and return a vector of data    
def overlap(c, sigma0, lat, lon, data_flat_total=[], lat_flat_total=[], lon_flat_total=[]):
    if c.geom_type =='MultiPolygon':
        for one_poly in c:
            minx, miny, maxx, maxy = one_poly.bounds
            print minx, miny, maxx, maxy
            print np.min(np.min(lat)), np.min(np.min(lon)), np.max(np.max(lat)), np.max(np.max(lon))
            cond1 = np.logical_and(lat>minx, lat<maxx)
            cond2 = np.logical_and(lon>miny, lon<maxy)
            
            ind = np.where(np.logical_and(cond1, cond2))
            lat_flat = lat[zip(ind)].ravel();
            lon_flat = lon[zip(ind)].ravel();
            data_flat = sigma0[zip(ind)].ravel()           
            
            points_collection =[Point(lat_i, lon_i) for lat_i, lon_i in zip(lat_flat, lon_flat)]
            for i_point, point in enumerate(points_collection):
                if one_poly.contains(point)==True:
                    data_flat_total.append(data_flat[i_point])
                    lat_flat_total.append(lat_flat[i_point])
                    lon_flat_total.append(lon_flat[i_point])

    return data_flat_total, lat_flat_total, lon_flat_total



# Create a regular grid with certain resolultion, knowing that:
# Resolution for latitudes: 1 deg N ~ 111 km; 0.1 ~ 11 km; 0.01 ~ 1 km; 0.001 ~ 100 m
# Resolution for longitudes is more complex close to the poles: 
# at 80N: 1deg ~ 19 km; 0.1 ~ 1.9 km; 0.01 ~ 100 m
# at 75N: 1deg ~ 29 km; 0.1 ~ 2.9 km; 0.01 ~ 290 m
# at 70N: 1deg ~ 38 km; 0.1 ~ 3.8 km; 0.01 ~ 380 m
def regular_grid(c, deg_lat=3, deg_lon=2):
# where deg_lat, deg_lon are derges of decimal precision after the dot (3 decimal after the dot for deg_lat, 2 decimals after the dot for lon)
    dlat = 10**(-deg_lat); dlon = 10**(-deg_lon)
    overlap_bounds = c.bounds  # x-y bounding box is a (minx, miny, maxx, maxy) tuple == (min_lat, max_lon, max_lat, max_lon)
    #print overlap_bounds
    
    r_lat_vector =  np.arange(round(overlap_bounds[0], deg_lat),round(overlap_bounds[2], deg_lat), dlat)   
    r_lon_vector = np.arange(round(overlap_bounds[1], deg_lon), round(overlap_bounds[3], deg_lon), dlon)
    regular_lat, regular_lon = np.meshgrid(r_lat_vector, r_lon_vector)
   
    return regular_lat, regular_lon


# Function that interpolates 1D vectors into regular grid 
def regular_interpolate(data_flat_total, lat_flat_total, lon_flat_total, c, deg_lat=3, deg_lon=2):
    regular_lat, regular_lon = regular_grid(c, deg_lat, deg_lon)
    xy_to_interpolate=np.array([lat_flat_total, lon_flat_total]).T
    s1_overlap_flat = np.array(data_flat_total)
    
    # if method is "linear" : you finish with nans where there is no data
    # if "nearest": griddata fulfills the empty space with "stripes" of interpolated data 
    regular_s1 = interpolate.griddata(xy_to_interpolate, s1_overlap_flat, (regular_lat, regular_lon), method='linear')
    regular_s1 = ma.masked_where(np.isnan(regular_s1), regular_s1)

    if regular_s1==[]:
        print 'Attention!!! regular_s1 is empty!'
    return  regular_s1, regular_lat, regular_lon
    
    

def nc_write():
    data = data

  
    
    
    
    
    

if __name__ == '__main__':
    # Choose which mask to apply ('OSISAF', which starts from Sep, 21, 2015  or recompute from 'SSMI')
    case_ice='SSMI'

    # case_list=['EastSiberian_Chukchi']
    # case_list=['Svalbard']
    case_list=['Greenland'] 
    # case_list=['Kara', 'Canadian', 'CentralArctic', 'EastSiberian_Chukchi',]
    # case_list=['Kara','Greenland']

    # time interval for finding the minmax coordinate around buoy, in days

    for case in case_list:
        fpath='/home/cercache/users/atarasen/data/Buoys/Arctic-buoy-data/'+case+'/' 
        fdir=os.listdir(fpath)
        print fdir

        for fname in fdir:
            if fname[-4:]=='.dat':
    #         if fname=='itp76rawlocs.dat':

                print fname

    #       check if final directories exist
                fpathout_fig='/home3/homedir7/perso/atarasen/Figures'
                fpathout_base='/home/cercache/users/atarasen/data/Sentinel-1_analysis/S1/'+case+'/'
                if os.path.isdir(fpathout_base)==False:
                    os.system('mkdir %s' %(fpathout_base))

                fpathout='/home/cercache/users/atarasen/data/Sentinel-1_analysis/S1/'+case+'/'+fname[:-4]
                if os.path.isdir(fpathout)==False:
                    os.system('mkdir %s' %(fpathout))
                    print fpathout, 'is created'
    #                 nf.pause()


    #       import initial time series of the min-max coordinates around buoy + t1 and t2 of these "time windows"
                nf.tic()
                year, day, lon, lat = buoy_module.importtext(fpath, fname)
                lats, lons, days = buoy_module.buoy_daily(year, day, lon, lat)


                true_ind0=[]; true_ind1=[]
                lat_b=[]; lon_b=[]
                count_miz=0

                for n in range(0,len(days)):
    #                 print fname
                    print n, 'from', len(days)
                    day_str=days[n].strftime('%Y-%m-%d')
                    if case_ice=='SSMI':
    #                     print 'case SSMI', day_str

                        try:
    # open appropriate SSMI sea ice concentration mask (ice - concentration, with negative values over the land (-128), sea ice concentration is from 0 to 100, 
                            ice, lon_ice, lat_ice = ice_module.SSMI_ice_open(day_str)
                        except:
                            pass

    #     if the file with ice information is NOT empty
                        if ice!=[]:
                            ice_data=np.squeeze(ice)	#make the masked variable readable 2d array with filled data
    #                         print 'Ice data shape', ice_data.shape                       
    #                         print 'lats[n], lons[n]', lats[n], lons[n]


    # have to skip part of trajectory, cause there is no data in the 'pole hole: for SSMI instrument:
                            if lats[n] < 87.2:
                                true_ind_n = interpolate_ice(lon_ice, lat_ice, lats[n], lons[n], 0.1, 0.1)

                                true_ind0.append(true_ind_n[0])
                                true_ind1.append(true_ind_n[1])

                                lat_b.append(lat_ice[true_ind_n[0]])    
                                lon_b.append(lon_ice[true_ind_n[1]])                     


    # compare with MW data pixel by pixel
                                percent=30
                                area=12      # enlarge area by 4 pixels = 50 km from each side (distance btwn neigbough 
                                            # pixels is 12.5 km)
                                open_water_area=False
                                ice_area=False

                                buoy_area=np.zeros((area+1, area+1)) 
                                buoy_area_lon=np.zeros((area+1, area+1)); buoy_area_lat=np.zeros((area+1, area+1))

                                for i_count, ii in enumerate(np.arange(-1*area/2, area/2+1)):
                                    for j_count, jj in enumerate(np.arange(-1*area/2, area/2+1)):
                                        buoy_area_lon[i_count, j_count]=lon_ice[true_ind_n[0]+ii, true_ind_n[1]+jj]
                                        buoy_area_lat[i_count, j_count]=lat_ice[true_ind_n[0]+ii, true_ind_n[1]+jj]  
                                        buoy_area[i_count, j_count]=ice_data[true_ind_n[0]+ii, true_ind_n[1]+jj]

                                        if (buoy_area[i_count, j_count]<=percent)==True:
    #                                         print ('Open water (C<%s) is found around the bouy!'%(percent))
                                            open_water_area=True
                                        elif np.logical_and(buoy_area[i_count, j_count]>=percent, \
                                                            buoy_area[i_count, j_count]<=100)==True:
    #                                         print ('sea ice is found around the bouy, yuppi!')   
                                            ice_area=True                                     

                                if open_water_area==True and ice_area==True:
                                    print ('we are not far from MIZ!')
                                    count_miz=count_miz+1
                                    lonmin_s=str(int(np.min(buoy_area_lon))); 
                                    lonmax_s=str(int(np.ceil(np.max(buoy_area_lon))))

                                    if lonmax_s=='360':
                                        lonmin_s='0'; lonmax_s='1'

                                    latmin_s=str(int(np.min(buoy_area_lat))); 
                                    latmax_s=str(int(np.ceil(np.max(buoy_area_lat))))

                                    print lonmin_s, lonmax_s, latmin_s, latmax_s



                                    print 'download S1 data with peps'
                                    print 'fpathout',fpathout, 
                                    print 'lonmin_s', lonmin_s,  'lonmax_s',lonmax_s
                                    print 'latmin_s', latmin_s, 'latmax_s', latmax_s, 
                                    print 'day_str', day_str, 'day_str', day_str
                                    os.system('/home/cercache/users/atarasen/codes/python/peps/peps_download.py -w %s -c S1 --lonmin %s --lonmax %s --latmin %s --latmax %s -a /home/cercache/users/atarasen/codes/python/peps/peps.txt -d %s -f %s'%(fpathout, lonmin_s, lonmax_s, latmin_s, latmax_s, day_str, day_str))