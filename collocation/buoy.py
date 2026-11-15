
# coding: utf-8

# In[37]:

#==================
# This script helps to extract the necessary time-space range to find the appropriate Sentinel-1 data for ITP buoys' positions
# It uses the *.dat file with bouy positions as an input and creates several outputs .picle files with min and max latitudes/longitudes over selected time window
# Buoy positions data has one hour time resolution (in general! it can save unusual gaps fro some reasons!!),
# thus, in order to have more chances to find appropriate S1 images, one can define the desired time window dt, inside of which the max and min coordinates will be found.
# By default the dt is set up to 1 day, because "peps.py" automatic S1 downloading script recognasie the start-moment and end-moment in days


# To read the whole file with buoy positions and time moments (as a python time structure) just apply:
# year, day, lon, lat = buoy_date_region.importtext(fpath, filename)
# struct_time = buoy_date_region.yearday2struct_time(year, day)

# plot_traj() - finally you can plot the trajectory of demanded buy over the Arctic

# author Anastasia Tarasenko
# creation: oct 2016
#==================


import numpy as np
import datetime
import pickle
import os, sys

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.basemap import Basemap



# 1. Function that extracts columns from a textfile
def importtext(fpath, filename):
  f = open(fpath+filename, 'r')
  headers=f.readline().strip().split()
  content= f.readlines()[0:-1]
  year=[]; day=[]; lon=[]; lat=[]
  
  for i in xrange(len(content)):
      line = content[i].strip()
      columns =line.split()
      #print columns
      if (np.array(columns, dtype=float)<-180).any()==False:
        if filename[0:3]=='itp':
#             print 'itp case'
            ind_year=0
            ind_day=1
            ind_lon=2
            ind_lat=3       
        else:
#             print 'non itp case'
            ind_lat = headers.index('Lat')
            ind_lon = headers.index('Lon')
            ind_year= headers.index('Year')
            if any('POS_DOY' in s for s in headers):
                ind_day = headers.index('POS_DOY')
            else:
                ind_day = headers.index('DOY')        
        #print 'longitude', columns[ind_lon]
        if float(columns[ind_lon])<0 and float(columns[ind_lon])>=-180:
            lon_buoy = float(columns[ind_lon])+360
        elif float(columns[ind_lon])>=0 and float(columns[ind_lon])<=180:
            lon_buoy = float(columns[ind_lon])

#           check for the errors in data:
	if float(columns[ind_day])>=0 and day==[]:
	      year.append(int(columns[ind_year]))
	      day.append(float(columns[ind_day]))
	      lon.append(lon_buoy)
	      lat.append(float(columns[ind_lat]))
	      
	elif float(columns[ind_day])>=0 and day!=[]:
	  if float(columns[ind_day])>day[-1] or (float(columns[ind_day])-day[-1])<-364:
	      year.append(int(columns[ind_year]))
	      day.append(float(columns[ind_day]))
	      lon.append(lon_buoy)
	      lat.append(float(columns[ind_lat]))
 
  f.close()
  return year, day, lon, lat


# In[52]:

#2. Recalculate the floats 'day of the year' + 'year' into python time structure
def yearday2struct_time(year, day):
  Nsec24h=60*60*24 #total number of seconds per 1 day
  
  Nyear=int(year)
  Nday=int(day) #the N of day without decimals (number of days without hours ans minutes)
  
  a=day-Nday	#the "decimal" part of the "day"
  sec=a*Nsec24h
  Nhour=int(sec/3600) #number of hours
  Nmin=int((sec-(Nhour*3600))/60) #nb of minutes
  Nsec=int(sec-(Nhour*3600)-(Nmin*60)) #nb of seconds
       
  #print "Nyear=", Nyear, " Nday=", Nday, " a=", a, " sec=", sec, " Nhour=", Nhour, " Nmin=", Nmin, " Nsec=", Nsec      
  #struct_time=datetime.datetime(Nyear, Nmonth, Nday, Nhour, Nmin, Nsec) # Look here for explanations: http://www.tutorialspoint.com/python/time_strptime.htm
  struct_time=datetime.datetime.strptime(str(Nyear)+' '+str(Nday)+' '+str(Nhour)+':'+str(Nmin)+':'+str(Nsec),'%Y %j %H:%M:%S' ) # Look here for explanations: http://www.tutorialspoint.com/python/time_strptime.htm

  return struct_time


# In[4]:

#3. Function to represent the structure data in string format
def struct2string(struct_time_i, ver='YYYY-mm-dd'):
  if ver=='YYYY-mm-dd HH:MM:SS': 
    date_string=struct_time_i.strftime('%Y-%m-%d %H:%M:%S')
  else:
    date_string=struct_time_i.strftime('%Y-%m-%d')
  return date_string


# In[74]:

# 4. Range the latlon data (position of your buoy) into "dt" days time windows.
# INPUT:
# year, day, lon, lat -  data on drifting buoy (outputs from importtext(filename))
# dt - time "window", in days, when you describe trajectory (can be only integer or float)
# OUTPUT:
# the format of lat-lon is: latmin, latmax, lonmin, lonmax = -35, -24, 27, 42
# start_date, end_date = "2015-12-01" + "T00:00:00.000Z", str(datetime.date.today()) + "T23:59:59.999Z"
# To enlarge the area aroung bouy, the spatial range is [-0.5 + bouy loc., +0.5 degree +bouy loc.], knowing that 0.5 degree ~ 50 km 

def latlontimerange(filename, year, day, lon, lat, dt=1):   
  latmax=[]; latmin=[]; lonmax=[]; lonmin=[]
  struct_time=[]; start_date=[]; end_date=[]
  
  len_initial=len(day)
  delta_t=datetime.timedelta(days=dt)  
  delta_1h=datetime.timedelta(minutes=60)
  
  # recalculate the day of the year + year in time_structure
  for i in range(0, len_initial):
    struct_time.append(yearday2struct_time(year[i], day[i]))  
    #if len(struct_time)>1:
      #print 'struct_time:', struct_time[0], struct_time[1]  
  
  
  # Why do you need the indexlist?
  ind_list=[]; ind1_list=[]; ind2_list=[]; t1_list=[]; t2_list=[]
  
  # Start your while loop from the beginning:
  n=0
  
  # Run over all except the last value:
  while n<(len_initial-1):
      print 'n=', n, 'len(struct_time)', len(struct_time)
      #print 't1-t2=', struct_time[n+1]-struct_time[n]
      
      # if the time difference between t[n+1] and t[n] is already larger that 1.5 delta_t,
      # move on to the next n, the time difference is too large to consider these positions as neibough (the buoy might hqve been switched off)
      if (struct_time[n+1]-struct_time[n])>(datetime.timedelta(dt*1.5)):
	   n=n+1
	   
      # if you're so lucky, and the t[n+1] - t[n] == delta_t
      # append and move to the next n
      elif (struct_time[n+1]-struct_time[n])==(delta_t):
		ind1_list.append(n)
		ind2_list.append(n+1)
		t1_list.append(struct_time[n])
		t2_list.append(struct_time[n+1])
		n=n+1
		
      # in any other case (when  t[n+1] - t[n] <1.5 delta_t
      # add ii parameter, and run the script only up to the moment when ii is less than (L-n-1)
      else:
	ii=2
	while ii<(len(struct_time)-n-1):
	  print 'ii=', ii, 'n+ii=', n+ii, len(struct_time)
	  
	  # if the difference (t[n+ii] - t[n]) is too large (larger than dt+1 hour)
	  # move to the next n and restart looking for neighbours - between the t(n+1) and t(n+2) the buoy stopped writing or crashed, or whatever...
	  
	  if (struct_time[n+ii]-struct_time[n])>(delta_t+delta_1h):
	    n=n+ii
	    break
	  
	  elif (struct_time[n+ii]-struct_time[n])<=(delta_t+delta_1h) and (struct_time[n+ii]-struct_time[n])>=(delta_t-delta_1h): #the difference btwn two neigbour "days" should be dt days
	    ind1_list.append(n)
	    ind2_list.append(n+ii)	    	    
	    t1_list.append(struct_time[n])
	    t2_list.append(struct_time[n+ii])
            print ' ' 
            print 't[n+ii]-t[n]', struct_time[n+ii]-struct_time[n]
            print struct_time[n+ii], ',', n+ii, 'line is appended'          
            print ' '   
     
            n=n+ii
            break
	  else:
	    ii=ii+1
        n=n+1	
  #print 'last n'
  ind1_list.append(len(struct_time)-1)
  ind2_list.append(len(struct_time))
  t1_list.append(struct_time[-2])
  t2_list.append(struct_time[-1])
  
  for ind in range(len(ind1_list)):
    ind_1=ind1_list[ind]
    ind_2=ind2_list[ind]
    #print ind_1, ind_2
    
    #print start_date[-1], end_date[-1]    
    
    latmin.append(min(lat[ind_1:ind_2])) #don't add+/-0.5 degree, or ~50 km from the buoy location
    latmax.append(max(lat[ind_1:ind_2]))
    
    lonmin.append(min(lon[ind_1:ind_2]))
    lonmax.append(max(lon[ind_1:ind_2]))
    
  return latmin, latmax, lonmin, lonmax, t1_list, t2_list 




# In[85]:

# 5. Average daily trajectories 
# recalculate the day of the year + year in time_structure

def buoy_daily(year, day, lon, lat):    
    struct_time=[]; d_n =[]; lon_new=[]
    
    for i in range(0, len(day)):
        struct_time.append(yearday2struct_time(year[i], day[i]))
        d_n.append(struct_time[i].replace(hour = 0, minute = 0, second = 0))
            
# Deal with 0-360 meridian passage: create a negative value for the 359:260 degree interval
# **The interval is so small, cause you compare two neighbour buoy positions, ans the bouy can not drift faster
# than 1 degree per day (~111 km)

        if lon[i]>330 and lon[i]<360:
            lon_new.append(lon[i] - 360)
        else:
            lon_new.append(lon[i])
            
# Write a new set of daily means of lats and lons      
    lat_n=[]; lon_n=[];
    lats=[]; lons=[]; days=[]
    for n in xrange(0, len(day)-1):
#         print 'n=', n, 'from', len(struct_time) 

# For all, but the last value:
        if n<len(day)-1:
            lat_n.append(lat[n])
            lon_n.append(lon_new[n])
            
            if (d_n[n+1]-d_n[n])>=datetime.timedelta(days=1):
                days.append(d_n[n]+datetime.timedelta(days=0.5))
                lats.append(np.mean(lat_n))
                lon_n_mean=np.mean(lon_n)
#  Go back to 0:360 degrees system for longitudes:      
                if lon_n_mean<0:
                    lons.append(lon_n_mean+360)
#                     print lons[-1]
                else:
                    lons.append(lon_n_mean)
                lat_n=[]; lon_n=[]
                
#  Append last values: mean of collected lat and lon
        else: 
            days.append(d_n[n]+datetime.timedelta(days=0.5)) 
            lats.append(np.mean(lat_n)); lat_n=[]            
            lon_n_mean=np.mean(lon_n); lon_n=[]
#  Go back to 0:360 degrees system for longitudes:
            if lon_n_mean<0:
                lons.append(lon_n_mean+360)
                print  'lon mean', lons[-1]
            else:
                lons.append(lon_n_mean)
                
# #    Last check: neigbour lons should not be different more than 1 degree:
#     days_new=[days[0],]; lons_new=[lons[0],]; lats_new=[lats[0],]
#     for i_new in range(1, len(days)):
#         if abs(lons[i_new]-lons[i_new-1])<1:
#             days_new.append(days[i_new])
#             lons_new.append(lons[i_new])
#             lats_new.append(lats[i_new])
            
#     return lats_new, lons_new, days_new
    return lats, lons, days


def plot_traj(case, buoy_case, fig_path_out='', daily=False):
    fpath='/home/cercache/users/atarasen/data/Buoys/Arctic-buoy-data/'+case+'/'
    if fig_path_out=='':
        fig_path_out='/home3/homedir7/perso/atarasen/Figures/'+case+'/'
        
    years, days, lons, lats = importtext(fpath, buoy_case+'.dat')
    
    if daily==True:
        lats, lons, buoy_days = buoy_daily(years, days, lons, lats)
    else:
        buoy_days=[]
        for ii, day in enumerate(days):
            buoy_days.append(yearday2struct_time(years[ii], day))
                
    buoy_lats=np.array(lats)
    buoy_lons=np.array(lons)

    
    fig_name=fig_path_out+buoy_case+'.png'
    fig = plt.figure(figsize=(12,12))

    proj='npstere'
    lat_min = -80.;  lat_max=81.;  lat_period=20.
    lon_min = -180.; lon_max=181.; lon_period=20.
    
    m = Basemap(projection=proj,boundinglat=60,lon_0=0,resolution='l')
    m.drawcoastlines()
    m.fillcontinents(color='white',lake_color='white')
    m.drawparallels(np.arange(0, 85, 10), labels=[True, False, False, False]) # labels = [left,right,top,bottom]
    m.drawmeridians(np.arange(-180, 180, 20),labels=[False, False, False, True])
    m.drawmapboundary(fill_color='white')
    
    x_b, y_b = m(buoy_lons, buoy_lats)
    m.plot(x_b, y_b)    

    fig_title = buoy_case+'\n'+struct2string(buoy_days[0])+' : '+struct2string(buoy_days[-1])
    plt.title(fig_title)
    
    plt.savefig(fig_name, dpi=300,  bbox_inches="tight")
    print fig_name, 'is created'
    plt.show()

    
    
    
    
    

# In[86]:

if __name__ == "__main__":
#     listnames=['itp37rawlocs.dat','itp76rawlocs.dat','itp83rawlocs.dat', 'itp91rawlocs.dat']
#     for areas=['Greenland', 'Svalbard', 'Kara']: 
 
    fpath='/home/cercache/users/atarasen/data/Buoys/Arctic-buoy-data/Greenland/'
  #listnames=['itp83rawlocs.dat']
  
#     fpath='/home/cercache/users/atarasen/data/Buoys/Arctic-buoy-data/EastSiberian_Chukchi/'
#     listnames=['300234060436000.dat']

 
  #fpath='/home/cercache/users/atarasen/data/Buoys/Arctic-buoy-data/Svalbard/'
  #listnames=['300234062428050.dat','itp37rawlocs.dat']
  #listnames=['itp37rawlocs.dat']

    dt=2
    fdir = os.listdir(fpath)
    for filename in fdir:
        if filename[-3:]=='dat':
            print filename
            year, day, lon, lat  = importtext(fpath, filename)
            lats, lons, days = buoy_daily(filename, year, day, lon, lat)
#             for i in range(0, len(days)):
#                 print 'Date:', days[i], 'lats', lats[i], 'lons', lons[i] 
    
#           start_date=[]; start_date_hours=[]; end_date=[]; end_date_hours=[]
#           for ind in range(len(t1_list)):
#             start_date.append(struct2string(t1_list[ind], 'YYYY-mm-dd'))
#             start_date_hours.append(struct2string(t1_list[ind], 'YYYY-mm-dd %H:%M:%S'))
#             end_date.append(struct2string(t2_list[ind], 'YYYY-mm-dd'))
#             end_date_hours.append(struct2string(t2_list[ind], 'YYYY-mm-dd %H:%M:%S'))

    ## Saving the objects:		    
            filenameout=fpath+filename[0:-4]+"-lat-lon-days.pickle"
            with open(filenameout, 'w') as f:
                pickle.dump([lats, lons, days], f)

#                 pickle.dump([latmin, latmax, lonmin, lonmax, t1_list, t2_list, start_date, start_date_hours, end_date, end_date_hours], f)
  
  
# Getting back the objects:
#import pickle
#with open('itp76-lat-lon-time.pickle') as f:
    #latmin, latmax, lonmin, lonmax, start_date, end_date = pickle.load(f)
    






