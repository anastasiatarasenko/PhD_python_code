# coding: utf-8

import datetime
import numpy as np

__author__   = 'Alexander Myasoedov'
__email__    = 'mag@rshu.ru'
__created__  = datetime.datetime(2016, 5, 11)
__modified__ = datetime.datetime(2016, 5, 11)
__version__  = "1.0"
__status__   = "Development"


# Перечислить все входные параметры:
# Разработать простую систему управления по типу Odyssey.
# Управление загрузками и обработкой данных.
# Возможность изменять:
#  фильтр по имени файла
#  входной/выходной пути
#  временной диапазон
#  разрешение
#  район по широтам/долготам
#  максимальный размер снимка

iPath = '/nfs1/store/satellite/asar/2012/036'
fileName = 'ASA_WSM_1PNPDK20120205_093002_000000923111_00223_51961_2779.N1'
pxlRes   = 800
image_size_skip = 140*1e6
area_bbox = 35
units     = 'dB'
oPath     = '/media/SOLabNFS/hyrax-cluster/data/public/allData/tzh/'
# When False append to NC wind and roughness and skip appending sigma0?
overwrite_nc = True
asar_conf = '/home/mag/Documents/repos/solab/PySOL/notebooks/pySAR/ASAR.json'
proj      = 'EPSG:3413'
numProcs  = 8

self, iPath, fileName, pxlRes=800.0, image_size_skip=40*1e6, area_bbox=45, units=None, oPath=None, overwrite_nc=False, asar_conf=None, proj='EPSG:3413', numProcs=None



class BaseHandler:
    def __init__(self):
        pass

    def findTimeIndices(self, time_variable):
        shape = time_variable.shape
        time_min = 0
        time_max = shape[0] - 1

        return (time_min, time_max)

    def buildClippingDapQuery(self, link, var_names, var_indices, variables_list):
        result = link + '.nc?'

        result += "%(lat_name)s[%(min)i:%(max)i]," % {'lat_name': var_names['lat'], 'min': var_indices['lat_min'], 'max': var_indices['lat_max']}
        result += "%(lon_name)s[%(min)i:%(max)i]," % {'lon_name': var_names['lon'], 'min': var_indices['lon_min'], 'max': var_indices['lon_max']}

        if var_names['time']:
            result += "%(time_name)s[%(min)i:%(max)i]," % { 'time_name': var_names['time'], 'min': var_indices['time_min'], 'max': var_indices['time_max'] }

        for idx, var in enumerate(variables_list):
            result += "%(var_name)s" % { 'var_name': var }

            if var_names['time']:
                result += "[%(time_min)i:%(time_max)i]" % { 'time_min': var_indices['time_min'], 'time_max': var_indices['time_max'] }

            result += "[%(lat_min)i:%(lat_max)i]" % { 'lat_min': var_indices['lat_min'], 'lat_max': var_indices['lat_max'] }
            result += "[%(lon_min)i:%(lon_max)i]" % { 'lon_min': var_indices['lon_min'], 'lon_max': var_indices['lon_max'] }

            if idx < len(variables_list) - 1:
                result += ","

        return result

    def get_lat_lon(self, dataset, variables_list=[]):
        if 'lat' in dataset.keys():
            lat_var = 'lat'
        elif 'latitude' in dataset.keys():
            lat_var = 'latitude'
        else:
            lat_var = ''

        if 'lon' in dataset.keys():
            lon_var = 'lon'
        elif 'longitude' in dataset.keys():
            lon_var = 'longitude'
        else:
            lon_var = ''

        if 'time' in dataset.keys():
            time_var = 'time'
        else:
            time_var = ''

        # remove variables missing in dataset
        for var in variables_list:
            if var not in dataset.keys():
                variables_list.remove(var)
                print "Missing {} variable. Removed.".format(var)

        if not variables_list:
            variables_list = dataset.keys()
            variables_list.remove(lat_var)
            variables_list.remove(lon_var)
            if variables_list.count(time_var):
                variables_list.remove(time_var)

        return lat_var, lon_var, time_var, variables_list

    def findSubsetIndices(self, min_lat, max_lat, min_lon, max_lon, lats, lons):
        """Array to store the results returned from the function"""
        res=np.zeros((4),dtype=np.float64)
        minLon=min_lon; maxLon=max_lon

        distances1 = []; distances2 = []
        indices=[]; index=1

        for point in lats:
            s1 = max_lat-point # (vector subtract)
            s2 = min_lat-point # (vector subtract)
            distances1.append((np.dot(s1, s1), point, index))
            distances2.append((np.dot(s2, s2), point, index-1))
            index=index+1

        distances1.sort()
        distances2.sort()
        indices.append(distances1[0])
        indices.append(distances2[0])

        distances1 = []; distances2 = []; index=1

        for point in lons:
            s1 = maxLon-point # (vector subtract)
            s2 = minLon-point # (vector subtract)
            distances1.append((np.dot(s1, s1), point, index))
            distances2.append((np.dot(s2, s2), point, index-1))
            index=index+1

        distances1.sort()
        distances2.sort()
        indices.append(distances1[0])
        indices.append(distances2[0])

        # Save final product: max_lat_indices,min_lat_indices,max_lon_indices,min_lon_indices
        minJ = indices[1][2]
        maxJ = indices[0][2]
        minI = indices[3][2]
        maxI = indices[2][2]

        res[0] = min(minI, maxI)
        res[1] = max(minI, maxI)
        res[2] = min(minJ, maxJ)
        res[3] = max(minJ, maxJ)

        return res

    def create_thumb(self, filepath):
        basewidth = 263
        img = Image.open(filepath)
        wpercent = (basewidth / float(img.size[0]))
        hsize = int(float(img.size[1]) * float(wpercent))
        img = img.resize((basewidth, hsize), PIL.Image.ANTIALIAS)
        savepath = filepath[:-4] + '_thumb' + filepath[-4:]
        img.save(savepath)

    def mkdir_p(self, path_dir):
        try:
            os.makedirs(path_dir)
        except OSError as exc: # Python >2.5
            if exc.errno == errno.EEXIST and os.path.isdir(path_dir):
                pass
            else: raise


    def setColorbar(self, gtmName):

        cm = gmtColormap(fileName=gtmName, GMTPath='api/handlers/gmtcolorbar')
        return cm
        # cmap=LevelColormap(levels, cmap=cm)

        # cs = m.contourf(x,y,bathy,levels,
        #                    cmap=LevelColormap(levels,cmap=cm),
        #                    alpha=1.0,
        #                    extend='both')

    def getColorbar(self, product, pp_name):
        colorbar_path = os.path.join('api/handlers/palette', product+'_'+pp_name.lower()+'.pal.png')

        print "colorbar path: ", colorbar_path
        if not os.path.isfile(colorbar_path):
            return None
        img = Image.open(colorbar_path)

        img = img.convert("RGBA")
        datas = img.getdata()

        newData = []
        for item in datas:
            if item[0] == 255 and item[1] == 255 and item[2] == 255:
                newData.append((255, 255, 255, 180))
            else:
                newData.append(item)

        img.putdata(newData)
        return img

    def save_composed(self, full_img, cb_img, sol_img, bg_img, comp_file, small=False):
        if not cb_img:
            return

        if small:
            full_img = full_img.resize((1440, 720))

        w, h = full_img.size
        bg_img = bg_img.resize((w, h))

        # instancemethod im.paste(self, im, box=None, mask=None)
        bg_img.paste(full_img, (0, 0), full_img)
        bg_img.paste(sol_img, (10, full_img.size[1] - sol_img.size[1] - 10), sol_img)
        bg_img.paste(cb_img, ( full_img.size[0] - cb_img.size[0] - 30, full_img.size[1]/2 - cb_img.size[1]/2), cb_img)

        bg_img.save(comp_file)

    def get_bg_image(self, month):
        cb_img = Image.open(os.path.join('api/handlers/resources', "%d.jpg" % int(month)))
        return cb_img

    def fix_longutude(self, lon):
        if lon > 180:
            return lon - 360
        return lon
