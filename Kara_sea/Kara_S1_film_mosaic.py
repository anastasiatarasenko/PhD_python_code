#!/usr/bin/env python

# The library that does mosaic and collocations of S1 images
# author: Anastasia Tarasenko

# Idea: create a set from [lon and lat] values for both S1 image to check if they intersect (at least partly). 
# The set is not an N-dimentional matrix, it is more close to a list of unique numbers.
# To create one list from 2 coordinates we just add them to each other (the less time-consuming operation).

import numpy as np
import matplotlib.pyplot as plt
import datetime, time

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.basemap import Basemap

import scipy.io
from scipy import ndimage
import cv2


import os, sys
sys.path.append('/home/cercache/users/atarasen/codes/python/')
sys.path.append('/home3/homedir7/perso/atarasen/codes/python/')
sys.path.append('/home3/homedir7/perso/atarasen/codes/python/collocation')

import nice_functions2 as nf
#import nice_functions2 as nf2
#import ncdfView as nv

import plot_stereo
reload(plot_stereo)

import S1_open_process_and_visu as S1_open
reload (S1_open)


import buoy as buoy_module
reload (buoy_module)

def s1nametodatetime(S1_name):
    s1_datetime = datetime.datetime.strptime(S1_name[17:21]+' '+S1_name[21:23]+' '+S1_name[23:25]+' '+ \
                                        S1_name[26:28]+' '+S1_name[28:30]+' '+S1_name[30:32],\
                                        '%Y %m %d %H %M %S' ) # Look here for explanations: http://www.tutorialspoint.com/python/time_strptime.htm
    return s1_datetime


### Plot separately for HH and HV for general case without buoys (just all available data over the Kara sea, eg)



case_list=[ 'Kara']
conditions={'Oct2014': [datetime.datetime(2014,10,1), datetime.datetime(2014, 11, 15)], 'Jul2015':[datetime.datetime(2015, 7, 20), datetime.datetime(2015,7, 25)], 'Dec2015': [datetime.datetime(2015, 12, 16), datetime.datetime(2016, 1, 10)]}
s1_names_dic={}

for case in case_list:
        fpath_0='/home/cercache/users/atarasen/data/Sentinel-1_analysis/'+case+'/'  
        fdir_0=os.listdir(fpath_0)
        
# Loop 1: transform all names of S1 files into datetime structures, sort the dictionary by time of sensing.

        for i, fname in enumerate(fdir_0):
            if fname[-4:]=='.zip':
                s1_names_dic[s1nametodatetime(fname)]=fname
        print 'S1 names sorted'                    

        

#Loop 2: sort the datetime structures 
        S1_time_sorted = sorted(s1_names_dic.keys())
        dt = datetime.timedelta(hours=(10./60)) #time difference of two neighbours should be less than 1-2 minutes

        
# More easy is to create seperaty HH and HV images (anyway the polarized data is opened separately every time)        
        for pol in ['hv']:
            print 'polarization:', pol
            
            fpathout_fig='/home3/homedir7/perso/atarasen/Figures/'+case+'/'+pol+'/'
            S1_open.checkdir_and_create(fpathout_fig)
            
            
#Run over the loop of sorted times 
            for i, s1_datetime in enumerate(S1_time_sorted):
                print s1_names_dic[s1_datetime], i, ' from ', len(S1_time_sorted)
                if i <len(S1_time_sorted)-1 and s1_datetime-S1_time_sorted[i+1]==0:
                    pass
            
                if i == 0 or (i>0 and s1_datetime-S1_time_sorted[i-1]>dt):
                    fig = plt.figure(figsize=(16,8))
                    fig.subplots_adjust(left=0.1,right=0.9,bottom=0.1,top=0.9)

                    spacing=1000
                    fig_name=fpathout_fig+s1_datetime.strftime('%Y-%m-%d_%H-%M-%S')+pol+'_mosaic_'+str(spacing)+'m.png'
                    print "Open a new figure"
                    
                try:
                        sigma0, lat, lon180, inc, ele = S1_open.s1(fpath_0+'/'+s1_names_dic[s1_datetime], pol,\
                                                                   spacing_m=spacing)
                        lon = plot_stereo.s1_lon360(lon180)
                        data_corr, data_corr_lg = S1_open.s1_simple_processing(sigma0, ele, pol)
                        sob_sigma = S1_open.sobel_filter(sigma0, pol)

                        ax1 = fig.add_subplot(121)
                        fig_title1='pol = '+pol
                        plot_stereo.region(data_corr_lg, lat, lon, s1_names_dic[s1_datetime][:-4], pol+'_log', \
                                           case=case, area_case=case, fig_title=fig_title1, buoy_case=None,\
                                           col_bar=True)
                        
                        ax2 = fig.add_subplot(122)  
                        fig_title2= pol + ' Sobel (Gaussian & d/dx +d/dy)'	
                        plot_stereo.region(sob_sigma, lat, lon, s1_names_dic[s1_datetime][:-4], pol+'_sobel', case=case, \
                                           area_case=case, fig_title=fig_title2, buoy_case=None,\
                                           col_bar=True)
                except:
                        pass

                if (i<(len(S1_time_sorted)-1) and S1_time_sorted[i+1]-S1_time_sorted[i]>dt) or i==len(S1_time_sorted)-1:
                    plt.suptitle('S1 '+s1_datetime.strftime('%Y-%m-%d %H:%M:%S')+'\n res: '+str(spacing)+' m',y=.95, size=16)
                    fig.savefig(fig_name, dpi=300,  bbox_inches="tight")
            #         pad_inches = 0)

                    print fig_name, ' is created'
                    #plt.show(fig)

             	    plt.close()

