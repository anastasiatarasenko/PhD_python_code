#!/usr/bin/env python
import os, sys
import os, sys
sys.path.append('/home3/homedir7/perso/atarasen/codes/python/posada/')
sys.path.append('/home3/homedir7/perso/atarasen/codes/python/cmod/')


import numpy as np

from posada.Handlers.sar import SAR
from posada.Models.swath import Swath
from posada.Readers.safe import SAFEGeoTiffFile
import matplotlib.pyplot as plt


dirname = '/home/cercache/users/atarasen/data/Sentinel-1_analysis/S1/Greenland/itp76rawlocs/'
filename ='S1A_EW_GRDM_1SDH_20141110T065604_20141110T065704_003216_003B66_E5F5.zip'
fn = dirname + filename
print(fn)

sarf = SAFEGeoTiffFile(url=fn, polarization='hh')
sarf.open()
sw = Swath()
sw.load(sarf)
sar = SAR(sw)

spacing_m = 400
spacing = np.int0(sar.meters2pixels(spacing_m))
#ssr = sar.get_data('roughness', spacing=spacing)
sigma0 = sar.get_data('sigma0', spacing=spacing)
ele = sar.get_data('elevation', spacing=spacing)
inc = sar.get_data('incidence', spacing=spacing)
lat = sar.get_data('lat', spacing=spacing)
lon = sar.get_data('lon', spacing=spacing)


