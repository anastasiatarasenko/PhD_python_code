#Retrieval of wind data from NCEP GFS database (NOAA forecast system)
# These data are provided 4 timed per day: at 00, 06, 12, 18 UTC
import numpy as np
from nansat import *

date='20151029'
T='0000'
n = Nansat('ncep_wind_online:'+date+T)
print n

minlon=65
maxlon=90

minlat=75
maxlat=80

imres_x=0.1
imres_y=0.1

d = Domain('+proj=longlat', '-te ' +str(minlon)+' '+str(minlat)+' '+str(maxlon)+' '+str(maxlat)+' -tr '+str(imres_x)+' '+str(imres_y))
n.reproject(d)

u = n['U']
v = n['V']
w = np.hypot(u,v)

# create canvas to draw a map
nmap = Nansatmap(n, resolution='l')

# make plot of wind speed (color) and add colorbar
nmap.pcolormesh(w)
nmap.add_colorbar(fontsize=10)

# add vectors showing wind directions
nmap.quiver(u, v, step=10)

# add continents on top
nmap.draw_continents()

# draw grid of parallels and meridians
# at given locations
parallels=np.arange(minlat, maxlat+1, 2.5)
meridians=np.arange(minlon, maxlon+1, 5)

nmap.drawparallels(parallels, labels=[1,0,0,0])
nmap.drawmeridians(meridians, labels=[0,0,1,0])

# add title
plt.suptitle('Wind speed and direction, '+str(date)+'T'+str(T))

# set size of the figure (inches)
nmap.fig.set_figheight(8)
nmap.fig.set_figwidth(20)

# save figure to a PNG file
figname='wind_'+str(date)+'T'+str(T)
nmap.save(figname+'.png')
