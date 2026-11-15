# coding: utf-8

import datetime

__author__ = 'Alexander Myasoedov'
__email__ = 'mag@rshu.ru'
__created__ = datetime.datetime(2016, 1, 1)


def findSubsetIndices(min_lat, max_lat, min_lon, max_lon, lats, lons):
    """Array to store the results returned from the function"""
    res = zeros((4), dtype=float64)
    minLon = min_lon
    maxLon = max_lon

    distances1 = []
    distances2 = []
    indices = []
    index = 1

    for point in lats:
        s1 = max_lat - point  # (vector subtract)
        s2 = min_lat - point  # (vector subtract)
        distances1.append((dot(s1, s1), point, index))
        distances2.append((dot(s2, s2), point, index - 1))
        index = index + 1

    distances1.sort()
    distances2.sort()
    indices.append(distances1[0])
    indices.append(distances2[0])

    distances1 = []
    distances2 = []
    index = 1

    for point in lons:
        s1 = maxLon - point  # (vector subtract)
        s2 = minLon - point  # (vector subtract)
        distances1.append((dot(s1, s1), point, index))
        distances2.append((dot(s2, s2), point, index - 1))
        index = index + 1

    distances1.sort()
    distances2.sort()
    indices.append(distances1[0])
    indices.append(distances2[0])

    # Save final product:
    # max_lat_indices,min_lat_indices,max_lon_indices,min_lon_indices
    minJ = indices[1][2]
    maxJ = indices[0][2]
    minI = indices[3][2]
    maxI = indices[2][2]

    res[0] = min(minI, maxI)
    res[1] = max(minI, maxI)
    res[2] = min(minJ, maxJ)
    res[3] = max(minJ, maxJ)

    return res


def crop():
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
