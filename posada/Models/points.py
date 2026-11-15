# coding: utf-8

import datetime

__author__ = 'Alexander Myasoedov'
__email__ = 'mag@rshu.ru'
__created__ = datetime.datetime(2016, 6, 29)


def swap(array):
    return concatenate([array[:, array.shape[1] / 2:],
                        array[:, :array.shape[1] / 2]], axis=1)


def wind(wind_source='NCEP',
         wind_iPath='/media/SOLabNFS2/store/model/ncep/gfs/',
         wind_FileName=None,
         dateTime=None, lats_2=None, lons_2=None, bbox=None
         logger=None):
    """
    HANDLER for any model/scatterometer wind for given time, lat/lon crop
    """
    if logger:
        logger.info("Reading %s wind.......", % wind_source)

    # Read
    if 'ERA_INTERIM' in wind_source:
        wind = windERAinterim(dateTime, lats_2, lons_2)
    elif 'NCEP_GFS' in wind_source:
        wind = windNCEPgfs(dateTime, lats_2, lons_2)

    # Crop/Process
    # find subset
    if lats_2 is not None and lons_2 is not None:
        bbox = [lats_2.min(), lats_2.max(), lons_2.min(), lons_2.max()]

    if bbox is not None:
        res = findSubsetIndices(bbox[0], bbox[1], bbox[2], bbox[3],
                                lats_wind[:, 0], lons_wind[0, :])
        # expand subset for better further pyresample
        # if close to the North Pole, crop all the data covering Pole
        if (res[0] == 0 and res[1] == lats_wind.shape[1] or
                lons_2.min() < -178 and lons_2.max() > 178):
            res[2] = 0
            res[3] = res[3] + 2
        else:
            res[0] = res[0] - 2
            res[1] = res[1] + 2
            res[2] = res[2] - 2
            res[3] = res[3] + 2

        # crop the data
        u = u[int(res[2]):int(res[3]), int(res[0]):int(res[1])]
        v = v[int(res[2]):int(res[3]), int(res[0]):int(res[1])]
        lats_wind = lats_wind[
            int(res[2]):int(res[3]), int(res[0]):int(res[1])]
        lons_wind = lons_wind[
            int(res[2]):int(res[3]), int(res[0]):int(res[1])]
    # Resample to swath/area_def

    # Export/Tile

    wind_def = pr.geometry.GridDefinition(lons=wind['lons_wind'],
                                          lats=wind['lats_wind'])
    swath_def = pr.geometry.SwathDefinition(lons=lons_2,
                                            lats=lats_2)
    pxlResWind = asarray(distancelib.getPixelResolution(wind['lats_wind'],
                                                        wind['lons_wind'],
                                                        wind['lats_wind'].shape, 'km'))
    windSwath = {}

    logger.info("Resampling wind to swath lats/lons.......")

    if 'gauss' in method:
        windSwath['wind_speed'] = pr.kd_tree.resample_gauss(wind_def, wind['wind_speed'].ravel(), swath_def,
                                                            radius_of_influence=2 * pxlResWind.max() * 1e3, neighbours=12,
                                                            sigmas=pxlResWind.max() * 1e3, fill_value=None, nprocs=numProcs)
        windSwath['wind_dir'] = pr.kd_tree.resample_gauss(wind_def, wind['wind_dir'].ravel(), swath_def,
                                                          radius_of_influence=2 * pxlResWind.max() * 1e3, neighbours=12,
                                                          sigmas=pxlResWind.max() * 1e3, fill_value=None, nprocs=numProcs)
        windSwath['u'] = pr.kd_tree.resample_gauss(wind_def, wind['u'].ravel(), swath_def,
                                                   radius_of_influence=2 * pxlResWind.max() * 1e3, neighbours=12,
                                                   sigmas=pxlResWind.max() * 1e3, fill_value=None, nprocs=numProcs)
        windSwath['v'] = pr.kd_tree.resample_gauss(wind_def, wind['v'].ravel(), swath_def,
                                                   radius_of_influence=2 * pxlResWind.max() * 1e3, neighbours=12,
                                                   sigmas=pxlResWind.max() * 1e3, fill_value=None, nprocs=numProcs)
    elif 'nearest' in method:
        windSwath['wind_speed'] = pr.kd_tree.resample_nearest(wind_def, wind['wind_speed'].ravel(), swath_def,
                                                              radius_of_influence=2 * pxlResWind.max() * 1e3,
                                                              fill_value=None, nprocs=numProcs)
        windSwath['wind_dir'] = pr.kd_tree.resample_nearest(wind_def, wind['wind_dir'].ravel(), swath_def,
                                                            radius_of_influence=2 * pxlResWind.max() * 1e3,
                                                            fill_value=None, nprocs=numProcs)
        windSwath['u'] = pr.kd_tree.resample_nearest(wind_def, wind['u'].ravel(), swath_def,
                                                     radius_of_influence=2 * pxlResWind.max() * 1e3,
                                                     fill_value=None, nprocs=numProcs)
        windSwath['v'] = pr.kd_tree.resample_nearest(wind_def, wind['v'].ravel(), swath_def,
                                                     radius_of_influence=2 * pxlResWind.max() * 1e3,
                                                     fill_value=None, nprocs=numProcs)

    windNCEPgfs = {}  # empty dict for windNCEPgfs

    if wind_FileName is None and dateTime:
        # find the ncep gfs filename to open from dateTime
        baseHour = floor((dateTime.hour + 3 / 2) / 6) * 6
        baseHour = min(18, baseHour)
        if dateTime.hour - baseHour > 1.5:
            forecastHour = 3
        else:
            forecastHour = 0

        if dateTime <= datetime.datetime(2014, 8, 19):
            wind_FileName = 'gfs' + dateTime.strftime("%Y%m%d") + '/gfs.t' + '%.2d' % (
                baseHour) + 'z.master.grbf' + '%.2d' % (forecastHour)

            grbs = pygrib.open(wind_iPath + wind_FileName)

            u_wind = None
            v_wind = None

            # wind contains u=u_wind.values[:], Lats=u_wind.latlons()[0],
            # Lons=u_wind.latlons()[1]
            for idx, msg_info in enumerate(grbs.select()):
                if msg_info['short_name'] == '10u':
                    u_wind = grbs.message(idx + 1)
                elif msg_info['short_name'] == '10v':
                    v_wind = grbs.message(idx + 1)

            u = u_wind.values[:]
            v = v_wind.values[:]
            lats_wind = u_wind.latlons()[0]
            lons_wind = u_wind.latlons()[1]
        else:
            try:
                wind_FileName = 'gfs.' + dateTime.strftime("%Y%m%d") + '%.2d' % (baseHour) + '/gfs.t' + '%.2d' % (
                    baseHour) + 'z.master.grbf' + '%.2d' % (forecastHour) + '.10m.uv.grib2'

                grbs = pygrib.open(wind_iPath + wind_FileName)

                u_wind = grbs.message(1)
                v_wind = grbs.message(2)
                u = u_wind['values']
                v = v_wind['values']
                lats_wind = u_wind['latitudes']
                lons_wind = u_wind['longitudes']
                lons_wind = reshape(lons_wind, (lons_wind.shape[0] / 720, 720))
                lats_wind = reshape(lats_wind, (lats_wind.shape[0] / 720, 720))
            except Exception:
                wind_FileName = 'gfs.' + dateTime.strftime("%Y%m%d") + '%.2d' % (
                    baseHour) + '/gfs.t' + '%.2d' % (baseHour) + 'z.pgrb2.0p25.f' + '%.3d' % (forecastHour)

                grbs = pygrib.open(wind_iPath + wind_FileName)

                u_wind = grbs.message(1)
                v_wind = grbs.message(2)
                u = u_wind['values']
                v = v_wind['values']
                lats_wind = u_wind['latitudes']
                lons_wind = u_wind['longitudes']
                lons_wind = reshape(
                    lons_wind, (lons_wind.shape[0] / 1440, 1440))
                lats_wind = reshape(
                    lats_wind, (lats_wind.shape[0] / 1440, 1440))

    # Make sure the longitude is between -180.00 .. 179.9
    lons_wind = map(
        lambda x:
        (lons_wind.ravel()[x] + 180) -
        int((lons_wind.ravel()[x] + 180) / 360) * 360 - 180,
        range(0, lons_wind.size)
    )
    lons_wind = reshape(lons_wind, lats_wind.shape)

    lats_wind = swap(lats_wind)
    lons_wind = swap(lons_wind)
    u = swap(u)
    v = swap(v)

    windNCEPgfs['lats_wind'] = lats_wind
    windNCEPgfs['lons_wind'] = lons_wind
    windNCEPgfs['wind_dir'] = windDirection(u, v)
    windNCEPgfs['wind_speed'] = sqrt(u**2 + v**2)
    windNCEPgfs['u'] = u
    windNCEPgfs['v'] = v
    windNCEPgfs['baseHour'] = baseHour
    windNCEPgfs['forecastHour'] = forecastHour

    return windNCEPgfs
