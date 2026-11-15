# s1 open basics and plot
# setenv PYTHONPATH /opt/python_libs/sar/gmf/:/opt/python_libs/sar/:/opt/python_libs/sar/sar/data 
import numpy as np
import nice_functions as nf

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import pdb

from cerbere.mapper.safegeotifffile import SAFEGeoTiffFile
from sar.data.sarimage import SARImage

import zipfile
import os
import datetime


os.sys.path.append('/opt/python_libs/sar/gmf')
os.sys.path.append('/opt/python_libs/sar')
os.sys.path.append('/opt/python_libs/sar/sar/data')

#Go to the specific library, check the data, extract if it is in zipfile
#path1='/media/atarasen/Transcend1/PhD_fr-ru/data/Sentinel-1_analysis/Kara-sea/Collection-of-features_Kara/Karskiye_vorota/automn/'

def checkzip(path1):
    data_dir=[]
    fdir=os.listdir(path1)
    for il, i in enumerate(fdir):
      if i[-4:]=='.zip':
        if os.path.isfile(i[-4:]+'.SAFE')==False:
#          print i
          zipfile.ZipFile(path1+i, 'r').extractall(path1)
#          data_dir.append(path1+i[0:-4]+'.SAFE/measurement/')
        else:
          print 'all files are unzipped' 
      elif i[-5:]=='.SAFE':
#          print i
          data_dir.append(path1+i+'/measurement/')
      else:
          pass       
    return data_dir           
          
path1='/export/home/atarasen/data/S1/Karskiye_vorota/case1/'
data_dir=checkzip(path1)

#Tiff open
for img in data_dir:  
    for ii in os.listdir(img): 
      print ii
#    for ii in enumerate(data_dir[1]): #case for check of one  hv imaige
      imagefile = SAFEGeoTiffFile(url=img+ii)
      surcouche = SARImage(imagefile)
      surcouche.get_datanames()
      
      pixelsspacing = surcouche.meters2pixels(1000) #1km grid spacing
      lons = surcouche.get_data('lon',spacing=pixelsspacing)
      lats = surcouche.get_data('lat',spacing=pixelsspacing)
      inc_angle = surcouche.get_data('incidence',spacing=pixelsspacing)
      sigma0 = surcouche.get_data('sigma0',spacing=pixelsspacing)
      
#pdb.set_trace()  
#%%      
#correct it    
#1. Find the mean along the azimuth to define the positions of the peaks
      plt.imshow(sigma0)
      plt.show()
      
      plt.imshow(inc_angle)
      plt.show()
      
      for ix in sigma0[:, 0].shape:
          sigma_mean=np.mean(sigma0,0)
      plt.plot(np.mean(inc_angle,0),sigma_mean)
      
      plt.title('Mean sigma_0, '+ii)
      plt.xlabel('incidence angle')
      plt.ylabel('sigma_0')
      plt.show()

#  
      #zone of interest
      cond1=np.logical_and(lons>55.4, lons<60.4)
      cond2=np.logical_and(lats>69.3, lats<70.5)
      [indx, indy]=np.where(np.logical_and(cond1, cond2))
      
      Sigma0_zone=nf.nans([sigma0.shape[0], sigma0.shape[1]])
      Sigma0_zone[indx, indy] = sigma0[indx, indy]
#      print np.shape(Sigma0_zone)
         
      
#      #plot the zone of interest
      plt.imshow(Sigma0_zone)
      plt.show()
      
#      plt.imshow(inc_angle)
#      plt.show()
#      
      #raw_input("Press the <ENTER> key to continue...")