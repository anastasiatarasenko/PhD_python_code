# Kara sea general case study:
# download all data for the formation period in 2014, eg.

import numpy as np
import matplotlib.pyplot as plt
import datetime, time

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.basemap import Basemap

import scipy.io
from scipy import ndimage

import os, sys
sys.path.append('/home/cercache/users/atarasen/codes/python/')
sys.path.append('/home3/homedir7/perso/atarasen/codes/python/')
sys.path.append('/home3/homedir7/perso/atarasen/codes/python/collocation')

import nice_functions2 as nf
reload(nf)
#import nice_functions2 as nf2
#import ncdfView as nv

import plot_stereo
reload(plot_stereo)

import S1_open_process_and_visu as S1_open
reload (S1_open)

import buoy as buoy_module
reload (buoy_module)

import datetime



# part one: download all the data for Kara Sea:
fpathout = '/home/cercache/users/atarasen/data/Sentinel-1_analysis/Kara'

date_case ='Jul2015'

#change cases manually:
if date_case=='Oct2014':
    lonmin_s='52'; lonmax_s='104';
    latmin_s='69'; latmax_s='81'
    day_str1 ='2014-10-01'; str1=datetime.datetime.strptime(day_str1, '%Y-%m-%d')
    n_days = 40
#day_str2 ='2014-10-10'; str2=datetime.datetime.strptime(day_str2; '%Y-%m-%d')


elif date_case =='Jul2015':
    lonmin_s='55.95'; lonmax_s='62.6';
    latmin_s='70.2'; latmax_s='75.3'
    day_str1 ='2015-07-20'; str1=datetime.datetime.strptime(day_str1, '%Y-%m-%d')
    n_days = 5

elif date_case =='Dec2015':
    #finish here for the small eddy in Kara sea
    n_days = 16
    
    
    
date_list = [str1 + datetime.timedelta(days=x) for x in range(0, n_days)]

for i, date in enumerate(date_list):
    if i<len(date_list)-1:
        date1=nf.struct2string(date) 
        date2=nf.struct2string(date_list[i+1])
    else: 
        date1=nf.struct2string(date, 'YYYY-mm-dd HH:MM:SS')
        date2=nf.struct2string(date+datetime.timedelta(hours=23.9))
        print date1
    os.system('/home/cercache/users/atarasen/codes/python/peps/peps_download.py -w %s \
-c S1 --lonmin %s --lonmax %s --latmin %s --latmax %s -a /home/cercache/users/atarasen/codes/python/peps/peps.txt \
-d %s -f %s'%(fpathout, lonmin_s, lonmax_s, latmin_s, latmax_s, date1, date2))

