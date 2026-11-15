#!/usr/bin/env python
import numpy as np

import os, sys
sys.path.append('/home/cercache/users/atarasen/codes/python/')
sys.path.append('/home3/homedir7/perso/atarasen/codes/python/')
import nice_functions as nf
import nice_functions2 as nf2
#import ncdfView as nv

import plot_stereo
reload(plot_stereo)

sys.path.append('/home3/homedir7/perso/atarasen/codes/python/collocation')
import buoy as buoy_module
reload (buoy_module)

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.basemap import Basemap, maskoceans

import scipy.io
from scipy import ndimage
import cv2


import datetime

sys.path.append('/home3/homedir7/perso/atarasen/codes/python/posada/')
sys.path.append('/home3/homedir7/perso/atarasen/codes/python/cmod/')


from posada.Handlers.sar import SAR
from posada.Models.swath import Swath
from posada.Readers.safe import SAFEGeoTiffFile



############################

def checkdir_and_create(dir_name):
    import os, sys
    if os.path.isdir(dir_name)==False:
        os.system('mkdir %s' %(dir_name))

def s1nametodatetime(S1_name):
    s1_datetime = datetime.datetime.strptime(S1_name[17:21]+' '+S1_name[21:23]+' '+S1_name[23:25]+' '+ \
                                        S1_name[26:28]+' '+S1_name[28:30]+' '+S1_name[30:32],\
                                        '%Y %m %d %H %M %S' ) # Look here for explanations: http://www.tutorialspoint.com/python
    return s1_datetime


def exclude_borders(fname, spacing_m=500):  
    if fname[10]=='H':
        coeff = spacing_m/25 # initial {high} resolution is 25 m
        ind_border = np.round(100/coeff)
    else:        
        coeff = spacing_m/40 # initial {medium} resolution is 40 m
        ind_border = np.round(100/coeff)
    if ind_border==0:
        ind_border_l = 0
        ind_border_r = -1
    else:
        ind_border_l = ind_border
        ind_border_r = -ind_border
        
    return ind_border_l, ind_border_r


def s1_get_corners(fn, pol, excl_bord = True, spacing_m = 500, print_cond=False):
# Find the geografical coordinates of the corner points of your image
    sarf = SAFEGeoTiffFile(url=fn, polarization=pol)
    
    sarf.open()
    sw = Swath()
    sw.load(sarf)
    sar = SAR(sw)
    
    a = sar.get_infos()
    if excl_bord==False: 
        coord1=[a['geolocation_grid']['latitude'][0,0], a['geolocation_grid']['longitude'][0,0]]
        coord2=[a['geolocation_grid']['latitude'][0,-1], a['geolocation_grid']['longitude'][0,-1]]
        coord3=[a['geolocation_grid']['latitude'][-1,-1], a['geolocation_grid']['longitude'][-1,-1]]
        coord4=[a['geolocation_grid']['latitude'][-1, 0], a['geolocation_grid']['longitude'][-1, 0]]
    else:
        ind_border_l, ind_border_r = exclude_borders(fn[-71:], spacing_m)
        coord1=[a['geolocation_grid']['latitude'][0, ind_border_l], a['geolocation_grid']['longitude'][0, ind_border_l]]
        coord2=[a['geolocation_grid']['latitude'][0, ind_border_r], a['geolocation_grid']['longitude'][0, ind_border_r]]
        coord3=[a['geolocation_grid']['latitude'][-1, ind_border_r], a['geolocation_grid']['longitude'][-1,ind_border_r]]
        coord4=[a['geolocation_grid']['latitude'][-1, ind_border_l], a['geolocation_grid']['longitude'][-1, ind_border_l]]
        
    corners=[coord1, coord2, coord3, coord4]
    lats = sum(corners,[])[::2]
    lons = sum(corners,[])[1::2]
    
    if print_cond==True:
        print lats
        print lons
    return corners, lats, lons



def s1(fn, pol='hh', spacing_m=500, excl_bord=False):
    sarf = SAFEGeoTiffFile(url=fn, polarization=pol)
    sarf.open()
    sw = Swath()
    sw.load(sarf)
    sar = SAR(sw)

    spacing = np.int0(sar.meters2pixels(spacing_m))
    #ssr = sar.get_data('roughness', spacing=spacing)
    sigma0 = sar.get_data('sigma0', spacing=spacing)
    
    ele = sar.get_data('elevation', spacing=spacing)
    inc = sar.get_data('incidence', spacing=spacing)
    lat = sar.get_data('lat', spacing=spacing)
    lon = sar.get_data('lon', spacing=spacing)
    
    if excl_bord ==True:
        ind_border_l, ind_border_r = exclude_borders(fn[-71:], spacing_m)
        sigma0 = sigma0[:, ind_border_l:ind_border_r]
        ele = ele[:, ind_border_l:ind_border_r]
        inc = inc[:, ind_border_l:ind_border_r]
        lat = lat[:, ind_border_l:ind_border_r]
        lon = lon[:, ind_border_l:ind_border_r]
           
    return sigma0, lat, lon, inc, ele
    

def s1_simple_processing(ice_data, elev_angle, pol='hh', corr=True):
#    for ix in sigma0[:, 0].shape:
#          sigma_mean=np.mean(sigma0,0) #sigma mean along the azimouth

    if corr==True:
        mat = scipy.io.loadmat('/home/cercache/users/atarasen/codes/matlab/NESZ_water.mat')
        elev_vector=np.ndarray.tolist(mat['elev_vector'][0])
        ele_step=elev_vector[1]-elev_vector[0]
        if pol == 'hh':
            nesz = np.ndarray.tolist(mat['nesz_HH'][0])
        elif pol =='hv':
            nesz = np.ndarray.tolist(mat['nesz_HV'][0])

        data_corr=np.ones(ice_data.shape)
        for i, v in enumerate(elev_vector):
            r, c = np.where(np.logical_and((v-ele_step/2)<elev_angle, elev_angle<=(v+ele_step/2)))
            if not(r.any())==False:
                data_corr[r, c] = ice_data[r, c] - nesz[i]
    else:
        data_corr=ice_data
        
    data_corr_lg = np.multiply(np.log10(data_corr), 10)
    
    return data_corr, data_corr_lg



def sobel_filter(sigma0, pol,  plot=False, S1_name=None):
    sigma=np.multiply(np.log(sigma0), 10)
    sx = ndimage.sobel(sigma, axis=0, mode='constant')
    sy = ndimage.sobel(sigma, axis=1, mode='constant')
    sob_sigma = np.hypot(sx, sy) # (sx**2+sy**2)**0.5
    
    if plot==True:
#         col_map = plt.get_cmap('gist_heat')
        col_map = plt.get_cmap('jet')

        plt.imshow(sob_sigma, cmap=col_map)
        plt.clim(0,50)
#         plt.colorbar()
        plt.title(S1_name[:32]+"..."+S1_name[-4:]+"\n pol = "+pol)
    #     plt.savefig(fig_name, dpi=300)
    #     plt.show()   
    
    
    return sob_sigma

def sobel_filter2(sigma0, pol, x=1, y=1, plot=False, S1_name=None):
    sigma=np.multiply(np.log(sigma0), 10)
    
    sob_sigma = cv2.Sobel(sigma,cv2.CV_64F,x,y,ksize=5)

    
    if plot==True:
#         col_map = plt.get_cmap('gist_heat')
        col_map = plt.get_cmap('jet')

        plt.imshow(sob_sigma, cmap=col_map)
        plt.clim(0,50)
#         plt.colorbar()
        plt.title(S1_name[:32]+"..."+S1_name[-4:]+"\n pol = "+pol)
    #     plt.savefig(fig_name, dpi=300)
    #     plt.show()   
    
    
    return sob_sigma

def laplacian_filter(sigma0, pol,  plot=False, S1_name=None):
    lapl_sigma = cv2.Laplacian(sigma0,cv2.CV_64F, ksize=5)
    
    if plot==True:
        col_map = plt.get_cmap('jet')
        plt.imshow(sob_sigma, cmap=col_map)
        plt.clim(0,50)
#         plt.colorbar()
        plt.title(S1_name[:32]+"..."+S1_name[-4:]+"\n pol = "+pol)
    #     plt.savefig(fig_name, dpi=300)
    #     plt.show()   
    
    return lapl_sigma




def mask_land(data, lat, lon, filled = True, fill_v = 0):
    masked_oceans = maskoceans(lon, lat, data, resolution = 'h')  
    data_masked_land=np.ma.masked_array(data, mask = np.logical_not(masked_oceans.mask))
    if filled == True:
        data_masked_land = np.ma.filled(data_masked_land, fill_value = fill_v)
    return data_masked_land






if __name__ == "__main__":
    case_list=['Greenland'] 
    # buoy_list=['itp76rawlocs', '300234060179710']
    for case in case_list:
        fpath_0='/home/cercache/users/atarasen/data/Sentinel-1_analysis/S1/'+case+'/'  #eg 'Greenland' 
        fpathout_fig='/home3/homedir7/perso/atarasen/Figures/'+case+'/'
        checkdir_and_create(fpathout_fig)

        fdir_case=os.listdir(fpath_0)
        #print 'fdir_case', fdir_case

    #     for fbuoy_dir in buoy_list:    
        for fbuoy_dir in fdir_case:
                print 'fbuoy_dir', fbuoy_dir            
                fdir=os.listdir(fpath_0+fbuoy_dir)
                #print fdir
                for i, fname in enumerate(fdir):
                    if fname[-4:]=='.zip':
                        print i, ' from ', len(fdir)
                        print fpath_0+fbuoy_dir+'/'+fname

                        nf.tic()
                        fig_name = fpathout_fig+fbuoy_dir+'_'+fname[17:31]+'_'+fname[-8:-4]+'.png'
                        if os.path.isfile(fig_name)==False or os.path.isfile(fig_name)==True:
                            fig = plt.figure(figsize=(12,12))
                            cond_hh=True; cond_hv=True

                            for pol in ['hh', 'hv']:
                                    print 'polarization:', pol
                                    if pol=='hh':
                                        try:
                                            sigma0, lat, lon180, inc, ele = s1(fpath_0+fbuoy_dir+'/'+fname, pol)
                                            lon = plot_stereo.s1_lon360(lon180)
                                            data_corr, data_corr_lg = s1_simple_processing(sigma0, ele, pol)

                                            ax = fig.add_subplot(221)
                                            plot_stereo.region(data_corr_lg, lat, lon, fname[:-4], pol+'_log',\
                                                               'Greenland', fbuoy_dir, col_bar=False)

                                            ax1=fig.add_subplot(222)
                                            fig_title1='pol = HH'
                                            x, y, x_b, y_b = plot_stereo.zoom(data_corr_lg, lat, lon, fname[:-4],\
                                                                              pol+'_log', 'Greenland', fbuoy_dir,\
                                                                              fig_title=fig_title1)

                                            ax2=fig.add_subplot(223)
                                            sob_sigma = sobel_filter(sigma0, pol)
                                            fig_title2='HH + Sobel filter (Gaussian & d/dx +d/dy)'
                                            plot_stereo.zoom(sob_sigma, lat, lon, fname[:-4], pol+'_sobel', \
                                                             'Greenland', fbuoy_dir, fig_title=fig_title2,\
                                                             x=x, y=y, x_b=x_b, y_b=y_b)
                                         
                                        except:
                                             cond_hh=False
                                             print 'error in reading', pol, 'polarization'
                                             pass

                                    elif pol=='hv':
                                         try:
                                            sigma0, _, _, _, _ = s1(fpath_0+fbuoy_dir+'/'+fname, pol)
                                            data_corr, data_corr_lg = s1_simple_processing(sigma0, ele, pol)

                                            ax3=fig.add_subplot(224)
                                            fig_title3='pol = HV'
                                            plot_stereo.zoom(data_corr_lg, lat, lon, fname[:-4], pol+'_log', \
                                                             'Greenland', fbuoy_dir, fig_title=fig_title3, \
                                                             x=x, y=y, x_b=x_b, y_b=y_b)
                                         except:
                                             cond_hv=False
                                             print 'error in reading', pol, 'polarization'
                                             pass
                            if cond_hh==True or cond_hv==True:
                                        plt.savefig(fig_name, dpi=300)
                                        print fig_name, ' is created'
                                        #plt.show()
                        plt.close()
                        nf.toc()
