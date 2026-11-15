# s1 open basics and plot
# setenv PYTHONPATH /opt/python_libs/sar/gmf/:/opt/python_libs/sar/:/opt/python_libs/sar/sar/data://home/cercache/users/atarasen/codes/python/
# setenv PYTHONPATH /home/cercache/users/atarasen/codes/python/ 

import os
import sys
sys.path.append('/opt/python_libs/sar/gmf/:/opt/python_libs/sar/:/opt/python_libs/sar/sar/data')
sys.path.append('/home/cercache/users/atarasen/codes/python/general_staff')
sys.path.append('/home/cercache/users/atarasen/codes/python/')

from cerbere.mapper.safegeotifffile import SAFEGeoTiffFile
from sar.data.sarimage import SARImage

import zipfile

import datetime

from nice_functions import tic
from nice_functions import toc

import matplotlib.pyplot as plt
import matplotlib.image as mpimg

#Go to the specific library, check the data, extract if it is in zipfile

def S1_raw(path1, fname_input, p=500):
    print fname_input
    if fname_input[-4:]=='.zip':
        tic()
        if os.path.isfile(fname[-4:]+'.SAFE')==False:
            zipfile.ZipFile(path1+fname, 'r').extractall(path1)
        fname=fname_input[:-4]+'.SAFE'
    else:
        print 'all files are unzipped'
        fname=fname_input
    toc()
    
  #Find all available measurements
  #elif fname[-4:]=='SAFE':
    path2=path1+fname+'/measurement/'
    data_dir=os.listdir(path2)
    #Tiff image open  
    for ii in data_dir:
      tic()
      tiffimg=path2+ii
      print tiffimg
      imagefile = SAFEGeoTiffFile(url=path2+ii)
      surcouche = SARImage(imagefile)
      surcouche.get_datanames()
      
      pixelsspacing = surcouche.meters2pixels(p) #500 m grid spacing
      lons = surcouche.get_data('lon',spacing=pixelsspacing)
      lats = surcouche.get_data('lat',spacing=pixelsspacing)
      inc_angle = surcouche.get_data('incidence',spacing=pixelsspacing)
      sigma0 = surcouche.get_data('sigma0')
      print sigma0
      
      fig = plt.figure()
      ax = fig.add_axes([0.1,0.1,0.8,0.8])
      
      plt.imshow(sigma0)
      
      nameimg=path1+ii[:-4]+'png'
      print "image name=",  nameimg
      plt.savefig(nameimg, bbox_inches='tight')
      toc()    
      return sigma0, lons, lats

def S1_draw(path1, data, lats, lons):
  fig = plt.figure()
  ax = fig.add_axes([0.1,0.1,0.8,0.8])
   
  plt.imshow(sigma0)
   
  nameimg=path1+ii[:-4]+'png'
  print "image name=",  nameimg
  plt.savefig(nameimg, bbox_inches='tight')
 
   #raw_input("Press the <ENTER> key to continue...")

if __name__ == "__main__": 
  path1='/home/cercache/users/atarasen/data/Sentinel-1_analysis/S1/Kara/Karskiye_vorota/case1/'
 #path1='/media/atarasen/Transcend2/PhD_fr-ru/data/Sentinel-1_analysis/Kara-sea/Collection-of-features_Kara/Karskiye_vorota/automn/case1/'
 #path1='/export/home/atarasen/data/S1/Karskiye_vorota/case1/'
 
  fdir=os.listdir(path1)
  for fname in fdir:
    sigma0, lons, lats = S1_raw(path1, fname)
    #S1_draw(path1,sigma0, lons, lats)
  