#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  7 21:27:01 2018

@author: anastasiiatarasenko
"""
import numpy as np
import scipy.signal as sg
import matplotlib.pyplot as plt



def filter_soak(depths, nminutes = 5, *times):
    #12 Hz (12 measurements per minute) comes from the fact that 
    # measuring frequency of the instruments on SBE 911 is 24 Hz, and thus
    # Nyquist freq (maximum sampling rate) is 1/2 *fmax
    try:
        mm = sg.argrelmin(depths, order = 12*nminutes) #indices for DEPTHS
        depths_min = depths[mm]         #new depths
        peaks_possible = np.where(np.logical_and(depths_min>5, depths_min<40)) #indices for NEW DEPTHS
        true_peaks_min = np.asarray(mm)[0][peaks_possible] #indices back in DEPTHS
        # 10 000 is an arbitrary value due to the soaking precedure in NABOS
        # 
        true_peaks_min = true_peaks_min[true_peaks_min<10000]

        plt.plot(true_peaks_min, depths[true_peaks_min], 'o')
        plt.plot(depths)
        plt.show()
        
        if len(true_peaks_min)>0:
            ind_last_soaking = true_peaks_min[-1]
            depths_cut = depths[ind_last_soaking:10000]
            mmax = sg.argrelmax(depths_cut, order = 12*nminutes)
            mmax = ind_last_soaking + mmax[0]
            depths_min = depths[mmax]
            
            plt.plot(true_peaks_min, depths[true_peaks_min], 'o')
            plt.plot(np.asarray(mmax), depths[mmax], 'o')
            plt.plot(depths)
            plt.show()
            
            v_min = np.nanmin(depths_min)
            ind_start = ind_last_soaking+ np.where(depths_cut==v_min)
            
            print(ind_start, v_min)
    except:
        print('regard inside the file, something went wrong')
        
    return ind_start, v_min
            
            
        
    
    
def some_filter(depth_full, cutoffF):
    freqs = np.arange(len(depth_full))/(len(depth_full)*1/24.)
    fNy = 1/(2*1/24.)
    
    dfft=np.fft.fft(depth_full)
    
    m = 2.*abs(dfft/len(depth_full))
    
    b, a = sg.butter(3, cutoffF/fNy)
    y = sg.filtfilt(b, a, depth_full)
    
    
    plt.plot(dfft[0:int(len(dfft)/2)], m[0:int(len(dfft)/2)])
    plt.show()
    
 
#    Too complex staff that is sad to rejet definetely
#def filter_soak(depths, nminutes = 4, ftitle = 'depths', ncounts_max_soaking=11000):
#    #12 Hz (12 measurements per minute) comes from the fact that 
#    # measuring frequency of the instruments on SBE 911 is 24 Hz, and thus
#    # Nyquist freq (maximum sampling rate) is 1/2 *fmax
##    try:
#        if len(depths)<20000:
#            ncounts_max_soaking = int(.3*(len(depths)))
#            print('The file contains very few measurements, look only the first '+str(ncounts_max_soaking))
#        mm = sg.argrelmax(depths[0:ncounts_max_soaking], order = 12*nminutes) #indices for DEPTHS
#        
#        depths_mm = depths[mm]         #new depths
#        peaks_possible = np.where(np.logical_and(depths[mm]>5, depths[mm]<40)) #indices for NEW DEPTHS
#        print('peaks possible', peaks_possible)
#        print('depths_mm[peaks_possible]', depths_mm[peaks_possible])
#        
#        true_peaks_max = np.asarray(mm)[0][peaks_possible] #indices back in DEPTHS
#        # 10 000 is an arbitrary value due to the soaking precedure in NABOS
#        print('true_peaks_max', true_peaks_max)
#            
#        if len(true_peaks_max)>1:
#                ind_depths_max = np.where(depths[true_peaks_max]==np.nanmax(depths[true_peaks_max])) #indices in DEPTH
#                print('depths[true_peaks_max]', depths[true_peaks_max])
#                print(np.nanmax(depths[true_peaks_max]))
#                print('ind_depths_max', ind_depths_max)
#        
#                true_peaks_max_i = true_peaks_max[ind_depths_max]
#                print('true_peaks_max_i', true_peaks_max_i)
#                print(depths[true_peaks_max_i])
#                    
#                plt.plot(true_peaks_max, depths[true_peaks_max], 'o')
#                plt.plot(depths[0:ncounts_max_soaking])
#                plt.plot(true_peaks_max_i, depths[true_peaks_max_i], '*')
#                plt.show()
#    
#    #            ind_last_soaking = true_peaks_min[-1]
#                ind_last_soaking = true_peaks_max_i[0]
#                depths_cut = depths[ind_last_soaking:ncounts_max_soaking]
#                print('ind_last_soaking:ncounts_max_soaking', ind_last_soaking, ncounts_max_soaking)
#                print('depths_cut', depths_cut)
#                plt.plot(depths_cut)
#                plt.title('depths_cut')
#                plt.show()
#
#                mmin = sg.argrelmin(depths_cut, order = 12*1)
#                print('all mmin', mmin)
#
#                mmin = ind_last_soaking + mmin
#                print(mmin)
#                depths_min = depths[mmin]
#                
#                print('depths_min', depths_min)
#                
#                fig_peaks = plt.figure()
#                plt.plot(true_peaks_max, depths[true_peaks_max], 'o')
#                plt.plot(true_peaks_max_i, depths[true_peaks_max_i], '*')
#                plt.plot(np.asarray(mmin), depths[mmin], 'o')
#                plt.plot(depths[0:ncounts_max_soaking])
#                plt.title(ftitle)
#                fig_peaks.savefig(fpath_figs+'/SBE911_'+st_name+'_first_depths.png')
#                plt.show()
#                
#                v_min = np.nanmin(depths_min)
#                print('v_min', v_min)
#                
#                if v_min > 5.:
#                    print('v_min > 5.')
#                    ind_start, v_min = 0, []
#                else:
#                    ind_start = ind_last_soaking + np.where(depths_cut==v_min)
#                    ind_start = ind_start[0][0]
#                print(ind_start, v_min)
#
#        else:
#            ind_start, v_min = 0, []
#            print('Seems that nothing was written during soaking!')    
#        return ind_start, v_min
