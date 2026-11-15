
# coding: utf-8

# In[3]:

# Kara sea general case study:
# download all data for the formation period in 2014, eg.

import numpy as np
import matplotlib.pyplot as plt
import datetime, time

import matplotlib.pyplot as plt

import os, sys
sys.path.append('/home/cercache/users/atarasen/codes/python/')
sys.path.append('/home3/homedir7/perso/atarasen/codes/python/')
sys.path.append('/home3/homedir7/perso/atarasen/codes/python/collocation')

import nice_functions as nf
#import nice_functions2 as nf2
#import ncdfView as nv

import plot_stereo
reload(plot_stereo)

import S1_open_process_and_visu as S1_open
reload (S1_open)

import collocation as colloc
reload(colloc)

import buoy as buoy_module
reload (buoy_module)

sys.path.append('/home3/homedir7/perso/atarasen/codes/python/posada/')
sys.path.append('/home3/homedir7/perso/atarasen/codes/python/cmod/')


from posada.Handlers.sar import SAR
from posada.Models.swath import Swath
from posada.Readers.safe import SAFEGeoTiffFile


# In[4]:

fpathout = '/home/cercache/users/atarasen/data/Sentinel-1_analysis/Kara/'


# In[ ]:

# part 1: download all the data for Kara Sea:
lonmin_s='52'; lonmax_s='104';
latmin_s='69'; latmax_s='81'
day_str1 ='2014-10-01'
day_str2 ='2014-11-10'


os.system('/home/cercache/users/atarasen/codes/python/peps/peps_download.py -w %s -c S1 --lonmin %s --lonmax %s --latmin %s --latmax %s -a /home/cercache/users/atarasen/codes/python/peps/peps.txt -d %s -f %s'%(fpathout, lonmin_s, lonmax_s, latmin_s, latmax_s, day_str1, day_str2))


# In[5]:

# part 2: 

# (1) check for the mosaic and collocation condition
fdir_0=sorted(os.listdir(fpathout))
pol= 'hh'

s1_names_dic={}
lat_dic={}
lon_dic={}
corners_dic={}

# Step 1.1: transform all names of S1 files into datetime structures, find all potential mosaics, 
# sort the dictionary by time of sensing.

nf.tic()
for counter, fname in enumerate(fdir_0):
    if fname[-4:]=='.zip':
        print fname, counter, 'of', len(fdir_0)
        try:
            
            corners, lats, lons = S1_open.s1_get_corners(fpathout+'/'+fname, pol)                
            if s1_names_dic=={}: #  Fist special case
                fname_1=fname[:]
                date_1=S1_open.s1nametodatetime(fname_1)            
                s1_names_dic[date_1]=[fname_1]
                
                lat_dic[date_1]=[lats]
                lon_dic[date_1]=[lons]
                corners_dic[date_1]=[corners]

            else:
                if colloc.mosaic_condition(fname_1, fname)==True:
                    s1_names_dic[date_1].append(fname)

                    lat_dic[date_1].append(lats)
                    lon_dic[date_1].append(lons)
                    corners_dic[date_1].append(corners)
                else:
                    fname_1=fname[:]
                    date_1=S1_open.s1nametodatetime(fname_1)
                    s1_names_dic[date_1]=[fname_1]
                    lat_dic[date_1]=[lats]
                    lon_dic[date_1]=[lons]

                    corners_dic[date_1]=[corners]
        except:
            pass

        
nf.toc()
        
dates = sorted(s1_names_dic.keys())         
print 'S1 names sorted'                    


# In[ ]:

datetime.datetime(2014, 10, 3, 3, 54, 3)


# In[ ]:

reload(colloc)
n_of_collocations= 0
pol = 'hh'
spacing=1000

coord_dic={}

# Step 1.2: Check if the images overlap:
for i_date, date in enumerate(dates):
    if i_date < (len(dates)-1):
        ii=1
        
        while ii<(len(dates)-1) and dates[i_date+ii]-dates[i_date]<=datetime.timedelta(hours=72):
            print date, dates[i_date+ii]
            
#             int_ans, coord = colloc.intersection_condition(corners_dic[date], corners_dic[dates[i_date+ii]], \
#                                                            s1_names_dic[date], s1_names_dic[dates[i_date+ii]],\
#                                                            plot_cond=True)
            int_ans, coord = colloc.intersection_condition(corners_dic[date], corners_dic[dates[i_date+ii]],                                                            plot_cond=True)

            if int_ans==True:
                print s1_names_dic[date], s1_names_dic[dates[i_date+ii]]
                print 'Successful overlap'

                coord_dic[date]={'coord':coord, 'collocated_with': dates[i_date+ii]}
                                 
                n_of_collocations = n_of_collocations + 1
                
                
                
                nf.tic()
#                 time to plot the overlapped days
                fig = plt.figure(figsize=(16,8))
                ax1 = fig.add_subplot(111)
                for mosaic in [s1_names_dic[date], s1_names_dic[dates[i_date+ii]]]:
                    for fn in mosaic:
                        print fpathout+fn
                        sigma0, lat, lon180, _, elev_angle = S1_open.s1(fpathout+fn, pol, spacing_m=spacing)
                        _, sigma_corr_lg = s1_simple_processing(sigma0, elev_angle, pol=pol)
                        lon = plot_stereo.s1_lon360(lon180)
                        
                        plot_stereo.region(sigma_corr_lg, lat, lon, fn[:-4], pol+'_log',                                                case='Kara', area_case='Kara', fig_title=fn, buoy_case=None,                                               col_bar=True)
                nf.toc()
                plt.show()
                
                break
            ii=ii+1
    
print n_of_collocations
print len(dates)
        
        


# In[ ]:

date=datetime.datetime(2014, 11, 4, 2, 48, 13)
print s1_names_dic[date]
print corners_dic[date]


# In[ ]:

## Version 2: 
# Open images one by one
# Step 1.1: transform all names of S1 files into datetime structures, find all potential mosaics, 
# sort the dictionary by time of sensing.

fdir_0=sorted(os.listdir(fpathout))
pol= 'hh'
s1_names_dic={}; lat_dic={}; lon_dic={}

nf.tic()

for counter, fname in enumerate(fdir_0):
    if fname[-4:]=='.zip':
        print fname, counter, 'of', len(fdir_0)
        corners = S1_open.s1_get_corners(fpathout+'/'+fname, pol)
        lats = sum(corners,[])[::2]
        lons = sum(corners,[])[1::2]

        date_1=S1_open.s1nametodatetime(fname)

        s1_names_dic[date_1]=fname
        lat_dic[date_1]=lats
        lon_dic[date_1]=lons
nf.toc()
       
dates = sorted(s1_names_dic.keys())         
print 'S1 names sorted'                    


# In[ ]:

## find first overlaps
reload(colloc)
spacing = 10000
case = 'Kara'
k= 0
# Step 1.2: Check if the images overlap:
for i_date, date in enumerate(dates):
    if i_date < (len(dates)-1):
        ii=1
        while ii<(len(dates)-1):
            int_ans = colloc.intersection_condition(lon_dic[date], lat_dic[date], lon_dic[dates[i_date+ii]], lat_dic[dates[i_date+ii]])
            if int_ans==True:
                print 'Successful overlap' 
                k = k+1
             
#                 time to plot the overlapped days

                fig = plt.figure(figsize=(16,8))
                for fn in [s1_names_dic[date], s1_names_dic[dates[i_date+ii]]]:
                    print fn
                    sigma0, lat, lon180, inc, ele = S1_open.s1(fpathout+'/'+fn, pol, spacing_m=spacing)
                    lon = plot_stereo.s1_lon360(lon180)
            
                    plot_stereo.region(sigma0, lat, lon, fn[:-4], pol+'_log',                                            case=case, area_case=case, fig_title=fn, buoy_case=None,                                           col_bar=True)
                plt.show()
                    
            
                break
    
print k
print len(dates)

# and then search for nearest "mosaics" for plotting


# In[ ]:

fname='S1A_EW_GRDH_1SDH_20141008T134704_20141008T134809_002739_003129_6176.zip'


# In[ ]:

n_files=len(S1_time_sorted)
n=0
while n<n_files:
    for pol in ['hh']:
        print 'polarization:', pol

#     fpathout_fig='/home3/homedir7/perso/atarasen/Figures/'+case+'/'+pol+'/'
#     S1_open.checkdir_and_create(fpathout_fig)


#Run over the loop of sorted times 
        for i, s1_datetime in enumerate(S1_time_sorted):
            print s1_names_dic[s1_datetime], i, ' from ', len(S1_time_sorted)
            if i <len(S1_time_sorted)-1 and s1_datetime-S1_time_sorted[i+1]==0:
                pass

            if i == 0 or (i>0 and s1_datetime-S1_time_sorted[i-1]>dt):

