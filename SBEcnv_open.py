#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 31 23:47:16 2018

@author: anastasiiatarasenko
"""
import os
import numpy as np
import scipy.signal as sg 

#import scipy
import datetime
import time
import pickle
import matplotlib.pyplot as plt


#fpath = '/Volumes/CAP_TARA/NABOS2018/CTD/'

def collect_head(content):
    st_name_CATS = ' '
    st_name_NABOS = ' '
    for i, aline in enumerate(content):
            if content[i][0]=='*':
                tempor = aline.split()
                try:
                    ind_lat = tempor.index('Latitude') + 2
                    lat_str = tempor[ind_lat] +' '+ tempor[ind_lat+1]
                except:
                    pass
                try:
                    ind_long = tempor.index('Longitude') + 2
                    long_str = tempor[ind_long] +' '+ tempor[ind_long+1]
                except:
                    pass
                try:
                    ind_time = tempor.index('(Time)') + 2
                    time_str = (tempor[ind_time] +' '+ tempor[ind_time+1]+' '+
                                   tempor[ind_time+2]+ ' '+ tempor[ind_time+3])
                except:
                    pass
                
                try:
                    tempor.index('Station:')
                    st_name = (tempor[-1])
                    print ('Station name: '+st_name)
                except:
                    pass
                
                try:
                    tempor.index('CATS:')
                    st_name_CATS = (tempor[-1])
                    print ('CATS name: '+st_name_CATS)
                except:    
                    pass
                
                try:
                    tempor.index('NABOS:') or tempor.index('nabos:') 
                    st_name_NABOS = (tempor[-1])
                    print ('NABOS name: '+st_name_NABOS)

                except:
                    pass
                
                
    return st_name, st_name_CATS, st_name_NABOS, lat_str, long_str, time_str

    
def latlon_str2num(lat_str, lon_str):
    lat = int(lat_str[0:2])+float(lat_str[3:5])/60.+.01*float(lat_str[5::])
    lon = int(lon_str[0:3])+float(lon_str[4:5])/60.+.01*float(lon_str[5::])
    return lat, lon


#def filter_jumping(var_, depth):
#    ind_max = np.where(depth == np.nanmax(depth))
#    if np.asarray(ind_max).shape[1]>1:
#        ind_max = ind_max[0][-1]
#        
#    grd = scipy.gradient(depth)
#    newdepth = masked
##    for idpth in enumerate(grd[0:ind_max]):
#         
#
#    fltrd_data = var_[np.where(np.ma.ismasked(new_depth)==False)]
#    
#    return fltrd_data
#    

def get_data(content):
    time_sec =[]; depSM = []
    t1 = []; t2 = []; s1 = []; s2 = []; turb = []; chla = []
    for i, aline in enumerate(content):
        if content[i][0]!='*' and content[i][0]!='#':
            tempor = aline.split()
            try:
                time_sec.append(float(tempor[1]))
                depSM.append(float(tempor[4]))
                t1.append(float(tempor[5]))
                t2.append(float(tempor[6]))
                s1.append(float(tempor[9]))
                s2.append(float(tempor[10]))
                turb.append(float(tempor[14]))
                chla.append(float(tempor[15]))
            except:
                print tempor
    for ivar in [depSM, t1, t2, s1, s2, turb, chla, time_sec]:
        ivar = np.asarray(ivar)  
    return depSM, t1, t2, s1, s2, turb, chla, time_sec




def filter_soak_simple(depths, ncounts_max_soaking=11000, st_name=''):
    fpath_figs = '/Users/anastasiiatarasenko/figures/SBE911_CTD'

    
    if len(depths)<20000:
        ncounts_max_soaking = int(.3*(len(depths)))
        print('The file contains very few measurements, look only the first '+str(ncounts_max_soaking))

    depths_cut = depths[0 : ncounts_max_soaking]
    inds_0, d_min = 0, depths[0]
#    d_min = np.nanmin(depths_cut)
#    print(d_min)
#    inds_0 = np.where(depths_cut == d_min)
#    print(inds_0)
#    if len(inds_0[0])>1:
#        inds_0 = inds_0[0][-1] #this is needed to have only one integer instead of given nested array [[ind]]
#    else:
#        inds_0 = inds_0[0][0]
    
    
    inds_20 = np.where(np.logical_and(depths_cut>20, depths_cut<24))
    fig_depths= plt.figure(figsize =(6, 4),)
    plt.plot(depths_cut)
    plt.plot(inds_20[0], depths_cut[inds_20], 'ro')
    
    inds_5 = np.where(np.logical_and(depths_cut>0, depths_cut<6))
    if len(inds_5[0])==0:
        print('Recording started after 6 m depth...')
#        print('Seems that soaking was cut already')
#        inds_0, d_min = 0, depths[0]
        
    else:
        plt.plot(inds_5[0], depths_cut[inds_5], 'bo')
        
#        For some reason when comparing all values with one,
#        if all the answers are negative, .any() gives True as well
#       - use of np.where is more robust
        ans = np.where(inds_20[0]<inds_5[0][-1])
#        print('np.where(inds_20[0]<inds_5[0][-1]?')
#        print(ans)
        if inds_20[0].any()<inds_5[0][-1] and len(ans[0])>0:
            print('soaking here!')
            d_min = np.nanmin(depths[inds_20[0][0] : ncounts_max_soaking])
            inds_0 = inds_20[0][0] + np.where(depths[inds_20[0][0]:ncounts_max_soaking]==d_min)
            print('all inds_0: ', inds_0)
            print('d_min: ', d_min)
            print('inds_0.shape:', inds_0.shape)
            if inds_0.shape[1]>1:
                inds_0 = inds_0[0][-1] #this is needed to have only one integer instead of given nested array [[ind]]
            else:
                inds_0 = inds_0[0][0]
            print('Final inds_0:', inds_0)
    
            plt.plot(inds_20[0][0], depths[inds_20[0][0]], 'm*')
        else:
            print('Seems that soaking was cut already')
#            inds_0, d_min = 0, depths[0]
    
    plt.plot(inds_0, depths[inds_0], 'y*')
    
    fig_depths.savefig(fpath_figs+'/SBE911_'+st_name+'_depths.png',
                                dpi = 300)
    plt.show()

    return inds_0, d_min 



def up_downcasts(depths, nminutes = 5):
    dmax = np.nanmax(depths)
    ind_bttm = np.where(depths==dmax)
    print(ind_bttm)
    if len(ind_bttm)>1:
        ind_bttm = ind_bttm[-1]
    ind_bttm = ind_bttm[0][0]
    return ind_bttm

#    sg.argrelmax(depth, order = 12*nminutes)
    
    
def get_10m_data(data, time_0, var_str, depth_min=5, depth_max=10, 
                 ind0 = 0, ind_btm = -1, filtered = False, *st_name):
    """ Collect data from depth_min to depth_max (upcast/downcast - everything inside the file)
    * time_0    time of CTD beginning (in python format)    
    * var_str   should be a name of desired variable
    * numbs     dictionary define the position of appropriate column 
         in data tuple extracted with get_data function
    """
    numbs={'t1': 1, 't2': 2, 's1': 3, 's2': 4, 'turb': 5, 'chla': 6}
    
    seconds = [datetime.timedelta(seconds = float(sec_tmp)) for sec_tmp in data[-1]]
    depths = np.asarray(data[0])
    var_array = data[numbs[var_str]]
    if filtered==True:
        var_array = sg.medfilt(var_array, 25)
   
    depths_tmp = depths[ind0:ind_btm]
    cond = np.logical_and(depths_tmp<=depth_max, 
                          depths_tmp>=depth_min)
    ans = np.where(cond)
    if len(ans[0])==0:
        print('Demanded depths: ['+str(depth_min)+', '+str(depth_max)+']')
        print('No such depths written @downcast, take the upcast')
        print('Recording starts from '+str(depths_tmp[ind0])+' m')
        ind0 = ind_btm
        ind_btm = -1
        
        depths_tmp = depths[ind0:ind_btm]
        plt.plot(depths_tmp)
        plt.show()
        cond = np.logical_and(depths_tmp<=depth_max, 
                          depths_tmp>=depth_min)
        print(ind0, ind_btm)
        print(depths_tmp)
    
    var_array = var_array[ind0:ind_btm]
    seconds = seconds[ind0:ind_btm]
    
    plt.figure()
    plt.plot(var_array, -1*depths_tmp)
    print(st_name)
    plt.title(var_str+' @ '+st_name[0])
    plt.show()
    
    
    depths_mmin_mmax = depths_tmp[np.where(cond)]
    var_= np.asarray(var_array)[np.where(cond)]
    times = time_0 + np.asarray(seconds)[np.where(cond)]
    return var_, depths_mmin_mmax, times
    
             

def write_stations_list(fpath, st_name, data_dict):    
    tmoment = datetime.datetime.now()
    tmoment_str = datetime.datetime.strftime(tmoment, '%y%m%d_%H%M') 
    
    fobj_stations = open(fpath + 'station_list' + tmoment_str + '.txt', 'w+')
    fobj_stations.write('Station name\t CATS name\t NABOS name\t Latitude\t Longitude\t Time\t Depth\n')
    for ii, i_name in enumerate(st_name):
        fobj_stations.write(i_name+'\t'+data_dict[i_name]['st_name_CATS']+'\t'+
                            data_dict[i_name]['st_name_NABOS']+'\t'+
                            data_dict[i_name]['lat_str']+'\t'+
                            data_dict[i_name]['lon_str']+'\t'+
                            data_dict[i_name]['time_str']+'\t'+'\n')
    fobj_stations.close()

    

            
    
if __name__=="__main__":
    #fpath = '/Volumes/K-T/NABOS-2018/CTD/'
    fpath = '/Users/anastasiiatarasenko/data/CTD/'
    fpath_figs = '/Users/anastasiiatarasenko/figures/SBE911_CTD'
    
    
    clim_v = {'t1': 'deg, C', 't2':'deg, C',
              's1': 'psu', 's2': 'psu', 'turb': 'NTU', 'chla':'ml/l'}
    
    data={}; 
    for fn in sorted(os.listdir(fpath)):
#    for fn in ['AT008.cnv',]:
#        'AT009.cnv',
#               'AT013.cnv', 'AT016.cnv', 'AT017.cnv',
#               'AT020.cnv', 'AT038.cnv',]:
        if fn[-4::]=='.cnv':
            print fn

            data5_10 = {}
            
            fobj = open(fpath+fn, 'r')
            content = fobj.readlines()
            fobj.close()

            data_strings = collect_head(content)

            st_name = data_strings[0]
            data[st_name]={}
            data[st_name]['st_name_CATS'] = data_strings[1]
            data[st_name]['st_name_NABOS']= data_strings[2]
            data[st_name]['lat_str'] =      data_strings[3] 
            data[st_name]['lon_str'] =      data_strings[4]
            data[st_name]['time_str'] =     data_strings[5]
            data[st_name]['time_py'] = datetime.datetime.strptime(data_strings[5], 
                                                                  '%b %d %Y %H:%M:%S')
            

# Get all the data
            data_tmp = get_data(content)
            depths = np.asarray(data_tmp[0])
#            ind_start, v_min = filter_soak(depths, nminutes = 4)
            ind_start, v_min = filter_soak_simple(depths, 11000, st_name)

            ind_btm = up_downcasts(depths, nminutes = 3)
            print(ind_start, ind_btm)
            data[st_name]['times'] = data_tmp[6]

           
            
## Get the depths in floating numbers and find all layers between 5 and 10 m            
            for var_nm in ['t1',
                           't2', 's1', 's2', 'turb', 'chla']:

                data5_10[var_nm], data5_10['depths'], data5_10['times'] = get_10m_data(data_tmp,
                                                                                       data[st_name]['time_py'], 
                                                                                       var_nm, 5, 10, 
                                                                                       ind_start, ind_btm,
                                                                                       True,
                                                                                       st_name)
                fig_var = plt.figure(figsize = (10, 6))
                plt.plot(data5_10[var_nm], -1*data5_10['depths'], '-*')
                plt.title(var_nm +' @ '+ st_name)
                plt.xlabel(clim_v[var_nm])
                fig_var.savefig(fpath_figs+'/SBE911_'+st_name+'_'+var_nm+'.png',
                                dpi = 300)

                plt.show()
#                
            data5_10['lat'], data5_10['lon'] = latlon_str2num(data[st_name]['lat_str'],  
                                                              data[st_name]['lon_str'])
## Write down a pickle files with
            fname_time = datetime.datetime.strftime(data[st_name]['time_py'], 
                                                  '%Y%m%dT%H%M%S')
            
            with open('/Users/anastasiiatarasenko/data/CTD/pickles/'+st_name+
                      '_'+fname_time +'.pickle', 'wb') as f:
                pickle.dump(data5_10, f)
                print('/Users/anastasiiatarasenko/data/CTD/pickles/'+st_name+
                      '_'+fname_time +'.pickle is written')
#            time.sleep(10)
#            
#
## Write an ascii-file with station names and coordinates            
#    write_stations_list(fpath, sorted(data.keys()), data)


    
    

        
