
# coding: utf-8

# In[4]:

# Kara sea general case study:
# download all data for the formation period in 2014, eg.

import numpy as np
from numpy import ma

import matplotlib.pyplot as plt
from compiler.ast import flatten

import datetime, time
import pickle


import matplotlib.pyplot as plt

import os, sys
sys.path.append('/home/cercache/users/atarasen/codes/python/')
sys.path.append('/home3/homedir7/perso/atarasen/codes/python/')
sys.path.append('/home3/homedir7/perso/atarasen/codes/python/collocation')

import nice_functions as nf
import nice_functions2 as nf2
reload(nf2)
#import ncdfView as nv

import plot_stereo
reload(plot_stereo)

import S1_open_process_and_visu as S1_open
reload (S1_open)

import collocation as colloc
reload(colloc)

import data_analysis as dt
reload(dt)

import buoy as buoy_module
reload (buoy_module)

import data_analysis
reload(data_analysis)

sys.path.append('/home3/homedir7/perso/atarasen/codes/python/posada/')
sys.path.append('/home3/homedir7/perso/atarasen/codes/python/cmod/')


from posada.Handlers.sar import SAR
from posada.Models.swath import Swath
from posada.Readers.safe import SAFEGeoTiffFile


fpath = '/home/cercache/users/atarasen/data/Sentinel-1_analysis/Kara/'
fpathout = '/home/cercache/users/atarasen/data/Sentinel-1_analysis/Kara/mosaics/'


# In[ ]:

# part 1: download all the data for Kara Sea:

lonmin_s='52'; lonmax_s='104';
latmin_s='69'; latmax_s='81'
day_str1 ='2015-10-01'
day_str2 ='2015-11-10'
# day_str1 ='2014-10-01'
# day_str2 ='2014-11-10'

os.system('/home/cercache/users/atarasen/codes/python/peps/peps_download.py -w %s -c S1 --lonmin %s --lonmax %s --latmin %s --latmax %s -a /home/cercache/users/atarasen/codes/python/peps/peps.txt -d %s -f %s'%(fpath, lonmin_s, lonmax_s, latmin_s, latmax_s, day_str1, day_str2))


# In[ ]:

# Part 2:  check for the mosaic and collocation condition

fdir_0=sorted(os.listdir(fpath))
pol= 'hh'

s1_names_dic={}
lat_dic={}; lon_dic={}
corners_dic={}

# Step 2.1: transform all names of S1 files into datetime structures, find all potential mosaics, 
# sort the dictionary by time of sensing.

#Small additional sorting of fnames by date, not taking into account if it is SSH/SDH or GRDH/GRDM:
nf.tic()
temp_dic={}
for fname_full in fdir_0:
        if fname_full[-4:]=='.zip':
            temp_dic[fname_full[17:]]=fname_full

list_of_files = sorted(temp_dic.keys())


# Real sorting:
for counter, name_temp in enumerate(list_of_files):
        fname = temp_dic[name_temp]
        print fname,  counter, 'of', len(fdir_0)
#     try:
        corners, lats, lons = S1_open.s1_get_corners(fpath+'/'+fname, pol, excl_bord = True, )                
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
                
#     except:
#         pass
                     
nf.toc()      
dates = sorted(s1_names_dic.keys())

with open(fpathout+'/Kara.pickle', 'wb') as handle:
    pickle.dump([dates, s1_names_dic, lat_dic, lon_dic, corners_dic], handle, protocol=pickle.HIGHEST_PROTOCOL)

print 'S1 names sorted'                    


# In[ ]:

#Part 3. Real collocation
reload(colloc)
reload(plot_stereo)
reload(data_analysis)

### Open processed dictionaries
with open(fpath+'/Kara.pickle', 'rb') as handle:
    dicts = pickle.load(handle)
    
dates=dicts[0]; s1_names_dic=dicts[1]
lat_dic=dicts[2]; lon_dic=dicts[3]
corners_dic=dicts[4]  


n_of_collocations= 0
pol = 'hv'
spacing=500
deg_lat=3; deg_lon=2

coord_dic={}


# Run over the collection of dates (mosaics)
for i_date, date in enumerate(dates):
    if i_date < (len(dates)-1):
        ii=1
        while ii<(len(dates)-1) and dates[i_date+ii]-dates[i_date]<=datetime.timedelta(hours=72):
            print date, dates[i_date+ii]

# Step 3.1: Check if the images overlap:
            int_ans, coord, c = colloc.intersection_condition(corners_dic[date], corners_dic[dates[i_date+ii]],                                                            plot_cond=True)        
            if int_ans==True:
                print 'Successful overlap:' 
                print s1_names_dic[date]
                print s1_names_dic[dates[i_date+ii]]
                print 'Time difference is ', date-(dates[i_date+ii])
                
                regular_lat, regular_lon = colloc.regular_grid(c, deg_lat=deg_lat, deg_lon=deg_lon)
                
                n_of_collocations = n_of_collocations + 1


                d1 = nf2.struct2string(date, ver='YYYY-mm-ddTHHMMSS')
                d2 = nf2.struct2string(dates[i_date+ii], ver='YYYY-mm-ddTHHMMSS')
                print d1, d2
                
                key_names = [d1, d2]
                dic_overlap={}; temp_dic_overlap={}               
                total_dic_overlap={}
                total_dic_overlap[d1, d2]={}
                
#                 try:                
# Preparations for plotting the overlapped days
                fig = plt.figure(figsize=(24,24))
                fig_name = '/home3/homedir7/perso/atarasen/Figures/Kara/overlaps/'+d1+'_'+d2+'_'+pol+'.png'

# 3.2 Open files from overlapping passages one by one                
                for i_mosaic, mosaic in enumerate([s1_names_dic[date], s1_names_dic[dates[i_date+ii]]]):
                    print 'Mosaic ', i_mosaic

                    temp_dic_overlap[i_mosaic]={}
                    temp_dic_overlap[i_mosaic]['s_1'] = []
                    temp_dic_overlap[i_mosaic]['lat'] = []; temp_dic_overlap[i_mosaic]['lon'] = []

                    key_name1 =  key_names[i_mosaic] 
                    dic_overlap[key_name1]={}
                    dic_overlap[key_name1]['s1_names'] = [s1_names_dic[date], s1_names_dic[dates[i_date+ii]]]
                    dic_overlap[key_name1]['s_1_regular']= {}
                    
                    
                    nf.tic()
# Open every S1 image in "mosaic"                    
                    for fn in mosaic:
                        print fpath+fn

                        sigma0, lat, lon180, _, elev_angle = S1_open.s1(fpath+fn, pol, spacing_m=spacing, excl_bord=True)
#                         _, sigma_corr_lg = S1_open.s1_simple_processing(sigma0, elev_angle, pol=pol, corr=False)    
                        lon = plot_stereo.s1_lon360(lon180)


                        s1_fn, lat_fn, lon_fn = colloc.contains_point_withshapely_v2(sigma0, lat, lon, c)
                        if s1_fn!=[]:
                            print "sucessful overlap!"
                            temp_dic_overlap[i_mosaic]['s_1']+=s1_fn
                            temp_dic_overlap[i_mosaic]['lat']+=lat_fn
                            temp_dic_overlap[i_mosaic]['lon']+=lon_fn

#                             plt.scatter(lon_fn, lat_fn, c=s1_fn)
#                             plt.show()

# Create a common regular grid and interpolate collected sigma_0 onto this grid

                    regular_s1, regular_lat, regular_lon = colloc.regular_interp(temp_dic_overlap[i_mosaic]['s_1'], temp_dic_overlap[i_mosaic]['lat'], temp_dic_overlap[i_mosaic]['lon'],                                                                                 c, regular_lat, regular_lon,                                                                                 deg_lat=deg_lat, deg_lon=deg_lon,                                                                                 mask_land=True)
                    dic_overlap[key_name1]['s_1_regular'][pol] = regular_s1
                    dic_overlap[key_name1]['lat_regular'] = regular_lat
                    dic_overlap[key_name1]['lon_regular'] = regular_lon

# Colect all information of collocated images into one total dicrionary
                    total_dic_overlap[d1, d2][key_name1] = dic_overlap[key_name1]

                    pickle_name = 'Kara_'+d1+'_'+d2+'_lat-'+str(deg_lat)+'_lon-'+str(deg_lon)+'_'+pol+'.pickle'
                    print pickle_name

                    with open(fpathout+'/'+pickle_name, 'wb') as handle:
                        pickle.dump(total_dic_overlap, handle, protocol=pickle.HIGHEST_PROTOCOL)                    

                    print  pickle_name, 'is created'   
#Plotting                

                    ax1 = fig.add_subplot(1, 2, i_mosaic+1)
                    data_log = 10*np.log10(dic_overlap[key_name1]['s_1_regular'][pol])
                    plot_stereo.zoom(data_log, dic_overlap[key_name1]['lat_regular'], dic_overlap[key_name1]['lon_regular'],                                  fn[:-4], pol='hv_log', case='Kara',                                 fig_title="mosaic no " + str(i_mosaic),                                 buoy_case=None,                                 col_bar=True, col_map_name='bone')

                    nf.toc()
                plt.savefig(fig_name, dpi=300)
                plt.show()                                   
#                     nf.pause()
                    
                    
#                 except:
#                     print 'no HV data found'

                break
            ii=ii+1
    
print n_of_collocations
print len(dates)


# with open(fpathout+'/Kara_Pol-'+pol+'_deg-lat-'+deg_lat+'_deg-lon-'+deg_lon+'.pickle', 'wb') as handle:
#     pickle.dump(dic_overlap, handle, protocol=pickle.HIGHEST_PROTOCOL)

print 'S1 names sorted'  
        


# In[1]:

dict_test = {}
dict_test[date, dates[i_date+ii]]


# In[ ]:

reload(S1_open)
fname = 'S1A_EW_GRDH_1SDH_20141005T100711_20141005T100811_002693_003022_63A6.zip'

#     'S1A_EW_GRDH_1SDH_20141003T003636_20141003T003740_002658_002F54_ECFA.zip'
#     'S1A_EW_GRDM_1SDH_20141003T134057_20141003T134157_002666_002F83_0B57.zip'
#     
#     'S1A_EW_GRDH_1SDH_20141003T102452_20141003T102553_002664_002F76_F611.zip'
#     'S1A_EW_GRDH_1SDH_20141003T102248_20141003T102352_002664_002F76_222F.zip'
#     'S1A_EW_GRDH_1SDH_20141003T102452_20141003T102553_002664_002F76_F611.zip'
fn = '/home/cercache/users/atarasen/data/Sentinel-1_analysis/Kara/'+fname
pol = 'hv'
spacing_m=100
sigma0, lat, lon180, _, ele = S1_open.s1(fn, pol, spacing_m=spacing_m, excl_bord=True)


#     sarf = SAFEGeoTiffFile(url=fn, polarization=pol)

#     sarf.open()
#     sw = Swath()
#     sw.load(sarf)
#     sar = SAR(sw)

#     spacing = np.int0(sar.meters2pixels(spacing_m))
#     #ssr = sar.get_data('roughness', spacing=spacing)
#     sigma0 = sar.get_data('sigma0', spacing=spacing)

#     ele = sar.get_data('elevation', spacing=spacing)
#     inc = sar.get_data('incidence', spacing=spacing)
#     lat = sar.get_data('lat', spacing=spacing)
#     lon = sar.get_data('lon', spacing=spacing)



# In[ ]:

print sigma0.shape
#     sigma_mean=[]; ele_mean=[]

# #     print sigma0[:, 0]
#     for i, ix in enumerate(np.array(sigma0)):

sigma_mean=np.mean(sigma0,0) #sigma mean along the azimouth
ele_mean= np.mean(ele, 0)

plt.plot(sigma_mean)
#     plt.plot(ele_mean, sigma_mean)
#     plt.plot(ele_mean[50:-50], sigma_mean[50:-50])

plt.show()

fft_result=np.fft.fft(sigma_mean)[range(len(sigma_mean)/2)]
freq=1/len(sigma_mean)
plt.plot(fft_result[0:500])
plt.show()


print np.median(sigma_mean), np.min(sigma_mean), len(sigma_mean)


#     sigma_new = sigma0[ind_left:ind_right]
plt.pcolormesh(10*np.log10(sigma0))
plt.show()


plt.pcolormesh(10*np.log10(sigma0[:, 0:100]))
plt.show()


# In[ ]:

#backup

#Part 3. Real collocation
reload(colloc)
reload(plot_stereo)
reload(data_analysis)

### Open processed dictionaries
with open(fpathout+'/Kara.pickle', 'rb') as handle:
    dicts = pickle.load(handle)
    
dates=dicts[0]; s1_names_dic=dicts[1]
lat_dic=dicts[2]; lon_dic=dicts[3]
corners_dic=dicts[4]  



n_of_collocations= 0
pol = 'hv'
spacing=500
deg_lat=3; deg_lon=2

coord_dic={}


# Run over the collection of dates (mosaics)
for i_date, date in enumerate(dates):
    if i_date < (len(dates)-1):
        ii=1
        while ii<(len(dates)-1) and dates[i_date+ii]-dates[i_date]<=datetime.timedelta(hours=72):
            print date, dates[i_date+ii]

# Step 3.1: Check if the images overlap:
            int_ans, coord, c = colloc.intersection_condition(corners_dic[date], corners_dic[dates[i_date+ii]],                                                            plot_cond=True)        
            if int_ans==True:
                print 'Successful overlap: '
                print s1_names_dic[date] 
                print 'with', s1_names_dic[dates[i_date+ii]]
                
                coord_dic[date]={'coord':coord, 'collocated_with': dates[i_date+ii]}                             
                n_of_collocations = n_of_collocations + 1
                               
                nf.tic()
                dic_overlap={}; temp_dic_overlap={}
                try:
                
    # Preparations for plotting the overlapped days
                    fig = plt.figure(figsize=(24,24))
                    fig_name = '/home3/homedir7/perso/atarasen/Figures/Kara/overlaps/'+                               nf2.struct2string(date, ver='YYYY-mm-ddTHHMMSS')+'_'+                               nf2.struct2string(dates[i_date+ii], ver='YYYY-mm-ddTHHMMSS')+'_'+pol+'.png'

                    key_names = [date, dates[i_date+ii]]
    # 3.2 Open files from overlapping passages one by one                
                    for i_mosaic, mosaic in enumerate([s1_names_dic[date], s1_names_dic[dates[i_date+ii]]]):
                        print 'Mosaic ', i_mosaic

                        temp_dic_overlap[i_mosaic]={}
                        temp_dic_overlap[i_mosaic]['s_1']=[]
                        temp_dic_overlap[i_mosaic]['lat']=[] 
                        temp_dic_overlap[i_mosaic]['lon']=[]

                        key_name1 =  key_names[i_mosaic] 
                        dic_overlap[key_name1]={}
                        dic_overlap[key_name1]['s1_names'] = [s1_names_dic[date], s1_names_dic[dates[i_date+ii]]]


    # Open every S1 image in "mosaic"                    
                        for fn in mosaic:
                            print fpathout+fn

                            nf.tic()
                            sigma0, lat, lon180, _, elev_angle = S1_open.s1(fpathout+fn, pol, spacing_m=spacing, excl_bord=True)
    #                         _, sigma_corr_lg = S1_open.s1_simple_processing(sigma0, elev_angle, pol=pol, corr=False)    
                            lon = plot_stereo.s1_lon360(lon180)


                            s1_fn, lat_fn, lon_fn = colloc.contains_point_withshapely_v2(sigma0, lat, lon, c)
                            if s1_fn!=[]:
                                print "sucessful overlap!"
                                temp_dic_overlap[i_mosaic]['s_1']+=s1_fn
                                temp_dic_overlap[i_mosaic]['lat']+=lat_fn
                                temp_dic_overlap[i_mosaic]['lon']+=lon_fn

    #                             plt.scatter(lon_fn, lat_fn, c=s1_fn)
    #                             plt.show()

    # Create a common regular grid and interpolate collected sigma_0 onto this grid

                        regular_s1_masked, regular_lat_masked, regular_lon_masked, reg_s1 = colloc.regular_interp(temp_dic_overlap[i_mosaic]['s_1'], temp_dic_overlap[i_mosaic]['lat'], temp_dic_overlap[i_mosaic]['lon'],                                                                                                                  c, deg_lat=deg_lat, deg_lon=deg_lon,                                                                                                                  mask_land=False)


    #                     print np.array(regular_s1_masked).shape
    #                     print np.array(reg_s1).shape
    #                     plt.contourf(regular_lon_masked, regular_lat_masked, reg_s1)
    #                     plt.colorbar()                    
    #                     plt.title('regular_lon_masked/ regular_lat_masked/ reg_s1')

                        dic_overlap[key_name1]['s_1_regular'][pol] = regular_s1_masked
                        dic_overlap[key_name1]['lat_regular'] = regular_lat_masked
                        dic_overlap[key_name1]['lon_regular'] = regular_lon_masked
                        np.pause()

                        ax1 = fig.add_subplot(1, 2, i_mosaic+1)
                        data_log = 10*np.log10(dic_overlap[key_name1]['s_1_regular'][pol])
                        plot_stereo.zoom(data_log, dic_overlap[key_name1]['lat_regular'], dic_overlap[key_name1]['lon_regular'],                                      fn[:-4], pol='hv_log', case='Kara',                                     fig_title="mosaic no " + str(i_mosaic),                                     buoy_case=None,                                     col_bar=True, col_map_name='bone')

                    plt.savefig(fig_name, dpi=300)
                    plt.show()



                    def_matrix = data_analysis.corr_fft(dic_overlap[date]['s_1_regular'][pol], dic_overlap[dates[i_date+ii]]['s_1_regular'])                
                    dic_overlap[key_names[0]]['cross_corr'] = abs(def_matrix).squeeze()

                    data_corr =  abs(def_matrix).squeeze()
                    data_corr[int(np.ceil(data_corr.shape[0]/2.)), int(np.ceil(data_corr.shape[1]/2.))] = 0
                    data_corr[int(np.floor(data_corr.shape[0]/2.)), int(np.floor(data_corr.shape[1]/2.))] = 0

                    fig = plt.figure(figsize=(12,12))
                    col_map = plt.get_cmap('jet')
                    plt.pcolormesh(abs(data_corr).squeeze(), cmap = col_map)

                    plt.colorbar()
    #                 plot_stereo.zoom(abs(def_matrix).squeeze(), regular_lat_masked, regular_lon_masked, fn[:-4], pol='', \
    #                                  case='Kara', fig_title='def matrix',\
    #                                  buoy_case=None,\
    #                                  col_bar=True, col_map_name='nipy_spectral') #clim = [0, 1e2]), clim = [0, 1e-4]
    # #                                        col_bar=True, col_map_name='jet')

                    fig_name_ifft = '/home3/homedir7/perso/atarasen/Figures/Kara/overlaps/'+                                   nf2.struct2string(date, ver='YYYY-mm-ddTHHMMSS')+'_'+                                   nf2.struct2string(dates[i_date+ii], ver='YYYY-mm-ddTHHMMSS')+'def_matr.png'

                    plt.savefig(fig_name_ifft, dpi=300)
                    plt.show()
                    
                    
                    
                    np.pause()
                    

                    
                except:
                    print 'no HV data found'

                nf.toc()

        
                break
            ii=ii+1
    
print n_of_collocations
print len(dates)


# In[ ]:


                    def_matrix = data_analysis.corr_fft(dic_overlap[date]['s_1_regular'][pol], dic_overlap[dates[i_date+ii]]['s_1_regular'])                
                    dic_overlap[key_names[0]]['cross_corr'] = abs(def_matrix).squeeze()

                    data_corr =  abs(def_matrix).squeeze()
                    data_corr[int(np.ceil(data_corr.shape[0]/2.)), int(np.ceil(data_corr.shape[1]/2.))] = 0
                    data_corr[int(np.floor(data_corr.shape[0]/2.)), int(np.floor(data_corr.shape[1]/2.))] = 0

                    fig = plt.figure(figsize=(12,12))
                    col_map = plt.get_cmap('jet')
                    plt.pcolormesh(abs(data_corr).squeeze(), cmap = col_map)

                    plt.colorbar()
    #                 plot_stereo.zoom(abs(def_matrix).squeeze(), regular_lat_masked, regular_lon_masked, fn[:-4], pol='', \
    #                                  case='Kara', fig_title='def matrix',\
    #                                  buoy_case=None,\
    #                                  col_bar=True, col_map_name='nipy_spectral') #clim = [0, 1e2]), clim = [0, 1e-4]
    # #                                        col_bar=True, col_map_name='jet')

                    fig_name_ifft = '/home3/homedir7/perso/atarasen/Figures/Kara/overlaps/'+                                   nf2.struct2string(date, ver='YYYY-mm-ddTHHMMSS')+'_'+                                   nf2.struct2string(dates[i_date+ii], ver='YYYY-mm-ddTHHMMSS')+'def_matr.png'

                    plt.savefig(fig_name_ifft, dpi=300)
                    plt.show()
                    
                    
with open(fpathout+'/Kara_Pol-'+pol+'_deg-lat-'+deg_lat+'_deg-lon-'+deg_lon+'.pickle', 'wb') as handle:
    pickle.dump([dic_overlap], handle, protocol=pickle.HIGHEST_PROTOCOL)

print 'S1 names sorted'  

