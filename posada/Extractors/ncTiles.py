import os
import sys

from ..Common.imports import \
    mean, ones, pr, MaskedArray, masked_equal, \
    sign, where, uint8, uint16, uint32, uint64, concatenate, array, asarray, \
    Polygon

import math
from netCDF4 import Dataset as ncDataset

try:
    # For Python 3.0 and later
    from urllib.request import urlopen
except ImportError:
    # Fall back to Python 2's urllib2
    from urllib2 import urlopen

import json

import warnings
warnings.simplefilter("ignore")


sys.path.append(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
)
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(
    os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__)))
    )
)
sys.path.append(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))))
    )
)

__author__ = 'Alexander Myasoedov'


def area_def_crop(area_def, new_shape):
    """
    Trimming the area_def correspondingly to the normalized array
    (x_ll, y_ll, x_ur, y_ur)
    """
    area_id = area_def.area_id
    name = area_def.name
    proj_id = area_def.proj_id
    proj4_args = area_def.proj_dict
    xsize = new_shape[1]
    ysize = new_shape[0]

    # area_def.get_lonlats()[0].mean()
    lon = mean([area_def.area_extent_ll[i] for i in [0, 2]])
    if (lon >= 90 and lon < 180) or (lon >= -270 and lon < -180):
        area_extent = (
            abs((new_shape[1]) * area_def.pixel_size_y) +
            area_def.area_extent[2],
            abs((new_shape[0]) * area_def.pixel_size_x) +
            area_def.area_extent[3],
            area_def.area_extent[2],
            area_def.area_extent[3])
    else:
        area_extent = (
            abs((new_shape[1]) * area_def.pixel_size_y) +
            area_def.area_extent[0],
            abs((new_shape[0]) * area_def.pixel_size_x) +
            area_def.area_extent[1],
            area_def.area_extent[0],
            area_def.area_extent[1])

    area_def_crop = pr.utils.get_area_def(area_id, name, proj_id, proj4_args,
                                          xsize, ysize, area_extent)

    return area_def_crop


def ncTilesMetainfo(model):

    pxlRes = asarray(
        (abs(model.area_def.pixel_size_x), abs(model.area_def.pixel_size_y))
    )

    cardinal_names = ['lon_min', 'lat_min', 'lon_max', 'lat_max']
    bbox = dict(zip(cardinal_names, model.get_bbox()))
    geo_extent = (
        (bbox['lon_min'], bbox['lat_min']),
        (bbox['lon_max'], bbox['lat_min']),
        (bbox['lon_max'], bbox['lat_max']),
        (bbox['lon_min'], bbox['lat_max']),
    )
    geospatial_extent_ll = Polygon(geo_extent)

    # get list of available resolutions
    avail_resolution_list = [int(256 * pow(2, i)) for i in range(15)]

    # get first guess size of area
    # size = int(256*math.pow(2, 14)/res)
    _size = abs(int(-5000000 - 5000000) / model.area_def.pixel_size_x)

    # get closest size of area from available resolutions list
    size = min(filter(lambda x: _size <= x,
                      avail_resolution_list))

    # Extent is large enough, not to get errors about points outside area
    _area_def = pr.geometry.AreaDefinition(
        'epsg_3413_crude', 'NSIDC Polar Stereographic North EPSG:3413',
        'epsg_3413_crude',
        {'proj': 'stere', 'lat_0': '90',
            'lon_0': '-45', 'lat_ts': '70', 'ellps': 'WGS84',
            'datum': 'WGS84', 'units': 'm'},
        size, size, [-7000000, -7000000, 7000000, 7000000])

    try:
        geospatial_extent = (
            _area_def.get_proj_coords(data_slice=(
                _area_def.get_xy_from_lonlat(geo_extent[0][0], geo_extent[0][1])[::-1])),
            _area_def.get_proj_coords(data_slice=(
                _area_def.get_xy_from_lonlat(geo_extent[1][0], geo_extent[1][1])[::-1])),
            _area_def.get_proj_coords(data_slice=(
                _area_def.get_xy_from_lonlat(geo_extent[2][0], geo_extent[2][1])[::-1])),
            _area_def.get_proj_coords(data_slice=(
                _area_def.get_xy_from_lonlat(geo_extent[3][0], geo_extent[3][1])[::-1])),
            _area_def.get_proj_coords(data_slice=(
                _area_def.get_xy_from_lonlat(geo_extent[0][0], geo_extent[0][1])[::-1]))
        )
    except:
        # if most southern lat is less 45 degrees
        # we take the most southern point available from the area_def
        geospatial_extent = (
            _area_def.get_proj_coords(data_slice=(
                _area_def.get_xy_from_lonlat(geo_extent[0][0], geo_extent[0][1])[::-1])),
            _area_def.get_proj_coords(data_slice=(
                _area_def.get_xy_from_lonlat(geo_extent[1][0], geo_extent[1][1])[::-1])),
            _area_def.get_proj_coords(data_slice=(
                (_area_def.shape[1] / 2, _area_def.shape[0]))),
            _area_def.get_proj_coords(data_slice=(
                (_area_def.shape[1] / 2, _area_def.shape[0]))),
            _area_def.get_proj_coords(data_slice=(
                _area_def.get_xy_from_lonlat(geo_extent[0][0], geo_extent[0][1])[::-1]))
        )

    return pxlRes, geospatial_extent, geospatial_extent_ll


def write_attrib_to_nc(nc_path, tile_extent, start_time, max_zoom_level,
                       resolution,
                       geospatial_extent=None,
                       geospatial_extent_ll=None,
                       polarizations=None,
                       variables=None,
                       dynamic_range=None):

    if os.path.isfile(nc_path):
        dataset = ncDataset(nc_path, 'a', format='NETCDF4')

    # minx, miny, maxx, maxy
    # Data extent
    if geospatial_extent:
        dataset.setncattr('geospatial_extent', geospatial_extent)
    if geospatial_extent_ll:
        dataset.setncattr('geospatial_extent_ll', str(geospatial_extent_ll))

    # Check extent to be float
    tile_extent = map(float, tile_extent)
    # Cropped tile extent
    dataset.setncattr('bbox_lat_min',
                      min(tile_extent[0], tile_extent[2]))
    dataset.setncattr('bbox_lon_min',
                      min(tile_extent[1], tile_extent[3]))
    dataset.setncattr('bbox_lat_max',
                      max(tile_extent[0], tile_extent[2]))
    dataset.setncattr('bbox_lon_max',
                      max(tile_extent[1], tile_extent[3]))

    dataset.setncattr('start_date', start_time.isoformat())

    if polarizations is not None:
        dataset.setncattr('polarizations', ','.join(polarizations))
    if variables is not None:
        dataset.setncattr('variables', ','.join(variables))
    dataset.setncattr('max_zoom_level', str(max_zoom_level))
    dataset.setncattr('resolution', str(resolution))
    if dynamic_range is not None:
        dataset.setncattr('dynamic_range', str(dynamic_range))

    dataset.close()


def array_size_normalize(input_array, is_land_masked=False):
    avail_resolution_list = [int(256 * math.pow(2, i)) for i in range(15)]
    # avail_resolution_list = [int(256*i) for i in range(15)]
    # logger.debug(avail_resolution_list)
    shape = input_array.shape
    # logger.debug('shape: ', shape)
    current_pixels = max(shape)
    # logger.debug('current pixels: ', current_pixels)
    necessary_pixels = min(filter(lambda x: current_pixels <= x,
                                  avail_resolution_list))
    # logger.debug('necessary_pixels: ', necessary_pixels)

    if current_pixels == necessary_pixels:
        return input_array

    if type(input_array) == MaskedArray:
        m = input_array.fill_value
    else:
        m = input_array.max() + 999999

    top_array = m * \
        ones((necessary_pixels - shape[0],
              shape[1]), dtype=input_array.dtype)
    right_array = m * \
        ones((necessary_pixels, necessary_pixels -
              shape[1]), dtype=input_array.dtype)
    input_array = concatenate((top_array, input_array), axis=0)
    input_array = concatenate((input_array, right_array), axis=1)

    input_array = masked_equal(input_array, m)

    # if is_land_masked:
    #     top_array   = m*ones((necessary_pixels-shape[0], shape[1]))
    #     right_array = m*ones((necessary_pixels, necessary_pixels-shape[1]))
    #     reprojected_array = concatenate((top_array, reprojected_array), axis=0)
    #     reprojected_array = concatenate((reprojected_array, right_array), axis=1)
    #     reprojected_array = masked_equal(reprojected_array, m)
    # if not is_land_masked:
    #     top_array   = m*ones((necessary_pixels-shape[0], shape[1]))
    #     right_array = m*ones((necessary_pixels, necessary_pixels-shape[1]))
    #     reprojected_array = concatenate((top_array, reprojected_array), axis=0)
    #     reprojected_array = concatenate((reprojected_array, right_array), axis=1)
    #     reprojected_array = masked_equal(reprojected_array, m)

    del top_array, right_array
    return input_array


def decrease_zoom(reprojected_array):
    """
    Scale image by 1/2 resolution for next zoom level
    # https://en.wikipedia.org/wiki/Pyramid_(image_processing)#Gaussian_pyramid
    """
    rows, cols = reprojected_array.shape
    rows_2 = rows / 2
    cols_2 = cols / 2
    sh = rows_2, rows // rows_2, cols_2, cols // cols_2
    return reprojected_array.reshape(sh).mean(-1).mean(1)


def generate_tiles_for_nc(input_array, zoom_level, _min=0, _max=4, nc_data_type='u1', logger=None):

    out_array = []

    for zoom in range(zoom_level, -1, -1):
        logger.debug('Tiling zoom level: %d' % zoom)
        # logger.debug(reprojected_array.shape)

        m = input_array
        m = where(m <= _min, _min + 0.0001, m)
        m = where(m >= _max, _max - sign(_max) * 0.0001, m)

        if nc_data_type == 'u1':
            m = (m - _min) / float(_max - _min) * (2**8 - 2)
            m = where(m.mask, (2**8 - 1), m)
            m = uint8(m)
        elif nc_data_type == 'u2':
            m = (m - _min) / float(_max - _min) * (2**16 - 2)
            m = where(m.mask, (2**16 - 1), m)
            m = uint16(m)
        elif nc_data_type == 'u4':
            m = (m - _min) / float(_max - _min) * (2**32 - 2)
            m = where(m.mask, (2**32 - 1), m)
            m = uint32(m)
        elif nc_data_type == 'u8':
            m = (m - _min) / float(_max - _min) * (2**64 - 2)
            m = where(m.mask, (2**64 - 1), m)
            m = uint64(m)
        out_array.append(m)

        # del m

        if zoom:
            input_array = decrease_zoom(input_array)

    return out_array


def create_dataset(output_path, max_zoom, nc_data_type='u1',
                   variables_list=None,
                   polarizations_list=None):
    if os.path.isfile(output_path):
        return ncDataset(output_path, 'a', format='NETCDF4')
    else:
        dataset = ncDataset(output_path, 'w', format='NETCDF4')

    if variables_list is not None:
        dataset.createDimension('vars', len(variables_list))

    if polarizations_list is not None:
        dataset.createDimension('polarizations', len(polarizations_list))

    dataset.createDimension('zoom', max_zoom + 1)
    dataset.createDimension('row', int(256 * math.pow(2, max_zoom)))
    dataset.createDimension('col', int(256 * math.pow(2, max_zoom)))

    if polarizations_list is not None and variables_list is not None:
        dims = ('vars', 'polarizations', 'zoom', 'row', 'col')
    elif polarizations_list is None and variables_list is not None:
        dims = ('vars', 'zoom', 'row', 'col')
    else:
        dims = ('zoom', 'row', 'col')

    if variables_list is not None:
        vars_var = dataset.createVariable('Variables', 'S1', ('vars',),
                                          zlib=True, complevel=6)
        vars_var[:] = array(variables_list, dtype='string')[:]

    if polarizations_list is not None:
        polar_var = dataset.createVariable('Polarizations', 'S1',
                                           ('polarizations',),
                                           zlib=True, complevel=6)
        polar_var[:] = array(polarizations_list, dtype='string')[:]

    # u1 = NC_UBYTE 0-255
    dataset.createVariable('Data', nc_data_type, dims, zlib=True, complevel=6)
    return dataset


def create_nc_tiles(reprojected_array, nc_path, nc_data_type=None, variable=None, polarization=None, configpath=None, logger=None):
    ''' Create NetCDF tile from preliminary reprojected array of data in given EPSG projection (e.g. 3413, 4326)

        Arguments:
        - reprojected_array: array to create tiles from.
        - variable:      variable name, e.g. sigma0, roughness, wind_speed, etc.
        - nc_path:       path to the output NetCDF
        - nc_data_type: data type, default 'u1', could be one of the following
                        'u1' (NC_UBYTE)  2**8 =0-255
                        'u2' (NC_USHORT  2**16=0-65535
                        'u4' (NC_UINT)   2**32=0-4294967295
                        'u8' (NC_UINT64) 2**64=0-18446744073709551615
        - polarization: for SAR data one can create tiles with several polarizations in one NetCDF file
        - configpath:   path to a JSON-type config file of a product

    '''

    # if no logger is specifyed, use the default one fromposada
    if not logger:
        from posada_logger import get_logger
        logger = get_logger()

    # Add possibility to read all the parameters from the JSON file!!!
    if not configpath:
        logger.warning(
            "No configuration path specified.\n"
            "Please specify it as: configpath = '/path/to/ASAR.json'")
        # default min_max_value
        min_max_value = (0, 255)
        polarization_num = None
        variable_num = None

    if configpath:
        try:
            config_dict = json.load(urlopen(configpath))
            variables_list = config_dict.get('variables_list', {})
            min_max_values = config_dict.get('min_max_values', {})
            # is_land_masked = config_dict.get('is_land_masked', {})
        except:
            logger.error('Could not load config')
            return

        if nc_data_type is None:
            nc_data_type = config_dict.get('nc_data_type', {})
            if len(nc_data_type) == 1:
                nc_data_type = nc_data_type[0]
            else:
                logger.warning(
                    "One have specified more than 1 datatype: %s" % nc_data_type)
        if 'polarizations_list' in config_dict:
            polarizations_list = config_dict.get('polarizations_list', {})

        if variable in variables_list:
            variable_num = variables_list.index(variable)
        else:
            logger.error('Variable not found in variables_list')
            return False

        if polarization is not None:
            polarization_num = polarizations_list.index(polarization)
        else:
            pass

        min_max_value = min_max_values[variable]

    zoom_level = int(math.log(reprojected_array.shape[0] / 256, 2))

    logger.debug('Number of zoom levels: %s', str(zoom_level))

    out_array = generate_tiles_for_nc(reprojected_array, zoom_level,
                                      _min=min_max_value[0],
                                      _max=min_max_value[1],
                                      nc_data_type=nc_data_type, logger=logger)

    logger.info('Writing netCDF-file')

    if configpath and 'polarizations_list' in config_dict and 'variables_list' in config_dict:
        dataset = create_dataset(
            nc_path, zoom_level, nc_data_type, variables_list, polarizations_list)
    elif configpath and 'variables_list' in config_dict and 'polarizations_list' not in config_dict:
        dataset = create_dataset(
            nc_path, zoom_level, nc_data_type, variables_list)
    else:
        dataset = create_dataset(nc_path, zoom_level, nc_data_type)

    datasetVar = dataset.variables['Data']

    if polarization is not None and variable is not None:
        for zoom in range(zoom_level, -1, -1):
            r = out_array[zoom_level - zoom].shape[0]
            c = out_array[zoom_level - zoom].shape[1]
            datasetVar[variable_num, polarization_num, zoom,
                       0:r, 0:c] = out_array[zoom_level - zoom][:]
    elif variable is not None and polarization is None:
        for zoom in range(zoom_level, -1, -1):
            r = out_array[zoom_level - zoom].shape[0]
            c = out_array[zoom_level - zoom].shape[1]
            datasetVar[variable_num, zoom, 0:r,
                       0:c] = out_array[zoom_level - zoom][:]
    else:
        for zoom in range(zoom_level, -1, -1):
            r = out_array[zoom_level - zoom].shape[0]
            c = out_array[zoom_level - zoom].shape[1]
            datasetVar[zoom, 0:r, 0:c] = out_array[zoom_level - zoom][:]

    dataset.close()

    return zoom_level
