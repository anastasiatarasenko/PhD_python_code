# coding: utf-8

"""
.. module::posada.Extractors.BaseExtractor

This module contains the abstract class for resampling/reprojecting and
extracting data to different file formats, as NetCDF, MAT, GeoTIFF etc...

:copyright: Copyright 2011-2016 SOLab.

.. sectionauthor:: Alexander Myasoedov <mag@rshu.ru>
.. codeauthor:: Alexander Myasoedov <mag@rshu.ru>
.. author:: SOLab Dev Team
"""

import os

from .ncTiles import \
    array_size_normalize, area_def_crop, create_nc_tiles, write_attrib_to_nc, \
    ncTilesMetainfo

from ..Models.Field import Field
from ..Extractors import netcdf

from ..Common.imports import \
    sign, where, uint8, uint16, uint32, uint64

import json

try:
    # For Python 3.0 and later
    from urllib.request import urlopen
except ImportError:
    # Fall back to Python 2's urllib2
    from urllib2 import urlopen

from copy import copy

from ..Common.logger import get_logger
logger = get_logger()

__all__ = ['extract']


def extract(model, fieldnames, url,
            configpath, nc_data_type='u1',
            append=True, mode='nctiles'):
    """
    Generic class for extracting data to different file formats, as
    NetCDF, MAT, GeoTIFF etc...

    Args:
        field (:class:'~posada.Models.Field.Field'): field to export.

        append (bool): if True then 'append', else - 'overwrite'.
            When 'append' mode is selected field(s) will be appended
            to the NetCDF file.
    Returns:
        None if exporting to file was a sucsess.
    """

    # check if nc-file exists
    if os.path.isfile(url) and append:
        logger.info("Appending to the NC-file")
        pass
    elif os.path.isfile(url) and not append:
        try:
            os.remove(url)
            logger.info("Overwriting the NC-file")
        except Exception:
            logger.error("Could not remove NC-file")
            return

    if 'nctiles' in mode:
        if not hasattr(model, 'area_def'):
            logger.error('Please specify the area_def in model.'
                         'Check that using resampleModel2proj() correctly')
        else:
            export2nctiles(model, fieldnames, url,
                           configpath, nc_data_type, append)
    else:
        return NotImplementedError


def export2nctiles(model,
                   fieldnames,
                   url,
                   configpath,
                   nc_data_type='u1', append=True):
    """
    Saves the data from field(s) class to SOLab ncTiler NetCDF file format.
    Data could be exported as uint8, uint16, uint32, uint64 and floats with
    specified scale and offset.

    Args:
        model (:class:'~posada.Models'): model to export.

        fieldnames (str): list of field names.

        url (str): full path to the NetCDF file
            where model fields will be saved.

        configpath (str): full path to the product json config file.

        nc_data_type (str): data type, default 'u1',
            could be one of the following:
            'u1' (NC_UBYTE)  2**8 =0-255
            'u2' (NC_USHORT  2**16=0-65535
            'u4' (NC_UINT)   2**32=0-4294967295
            'u8' (NC_UINT64) 2**64=0-18446744073709551615

        append (bool): if True then 'append', else - 'overwrite'.
            When 'append' mode is selected field(s) will be appended to
            the existing NetCDF file.
    Returns:
        None if exporting to file was a sucsess.
    """

    if not type(fieldnames) is list:
        fieldnames = [fieldnames]
    for fieldname in fieldnames:
        data = model.get_field(fieldname).get_values()

        # Normalize array size to be a multyply of the tile size
        data = array_size_normalize(data, False)
        new_shape = data.shape

        # Trimming the area_def correspondingly to the normalized array
        tile_area_def = area_def_crop(model.area_def, new_shape)

        max_zoom_level = create_nc_tiles(data,
                                         url,
                                         nc_data_type=nc_data_type,
                                         variable=fieldname,
                                         configpath=configpath,
                                         logger=logger)

    pxlRes, geospatial_extent, geospatial_extent_ll = ncTilesMetainfo(model)

    write_attrib_to_nc(url,
                       tile_area_def.area_extent,
                       model.get_datetimes()[0],
                       max_zoom_level,
                       pxlRes,
                       geospatial_extent, geospatial_extent_ll,
                       polarizations=None,
                       variables=fieldnames,
                       dynamic_range=None)

    return NotImplementedError


def compute_scale_and_offset(_min, _max, nc_data_type='u1'):
    """

    """
    if nc_data_type == 'u1':
        n = 8
    elif nc_data_type == 'u2':
        n = 16
    elif nc_data_type == 'u4':
        n = 32
    elif nc_data_type == 'u8':
        n = 64

    # stretch/compress data to the available packed range
    scale_factor = float(_max - _min) / (2**n - 2)
    # translate the range to be symmetric about zero
    add_offset = float(_min)
    # _FillValue
    _FillValue = float(2**n - 1)
    return (scale_factor, add_offset, _FillValue)


def export2SOLabnc(model,
                   fieldnames,
                   url,
                   configpath=None,
                   nc_data_type='u2',
                   attrs={'acknowledgement': 'MAGiCAN'},
                   append=True):
    """
    Saves the data from model class to specific NetCDF-file format
    used at SOLab.
    Data could be exported as uint8, uint16, uint32, uint64 and floats with
    specified scale and offset.

    Args:
        model (:class:'~posada.Models'): model to export.

        fieldnames (str): list of field names.

        url (str): full path to the NetCDF file
            where model fields will be saved.

        configpath (str): full path to the product json config file.

        nc_data_type (str): data type, default 'u1',
            could be one of the following:
            'u1' (NC_UBYTE)  2**8 =0-255
            'u2' (NC_USHORT  2**16=0-65535
            'u4' (NC_UINT)   2**32=0-4294967295
            'u8' (NC_UINT64) 2**64=0-18446744073709551615

        append (bool): if True then 'append', else - 'overwrite'.
            When 'append' mode is selected field(s) will be appended to
            the existing NetCDF file.
    Returns:
        None if exporting to file was a sucsess.

    """

    if not configpath:
        logger.warning(
            "No configuration path specified. Using model min_max_values.\n "
            "Other, specify configpath as: configpath='/path/to/ASAR.json'")
        return

    if not type(fieldnames) is list:
        fieldnames = [fieldnames]
    for fieldname in fieldnames:

    #

        # Getting valid min, max values for further data packing
        if configpath:
            try:
                # Fix. Not loaded in configpath. Python3
                with open(configpath) as json_file:
                    config_dict = json.load(json_file)
                #config_dict = json.load(urlopen(configpath))
                min_max_values = config_dict.get('min_max_values', {})
            except:
                logger.error('Could not load config')
                return
            try:
                min_max_value = min_max_values[fieldname]
            except:
                logger.error('Could not get "min_max_value" from config'
                             'for fieldname: %s' % fieldname)
                continue
            _min = float(min_max_value[0])
            _max = float(min_max_value[1])
        else:
            _min = float(model.get_field(fieldname).valid_min)
            _max = float(model.get_field(fieldname).valid_max)

        # Calculating scale_factor, offset and FillValue
        scale_factor, add_offset, _FillValue = compute_scale_and_offset(
            _min, _max, nc_data_type)

        # Get the data from field to be packed
        m = model.get_values(fieldname)

        # Packing values
        m = where(m <= _min, _min + 0.0001, m)
        m = where(m >= _max, _max - sign(_max) * 0.0001, m)

        # DO NOT apply scale and offset here,
        # as it is done during NetCDF write!
        # m = (m - add_offset) / scale_factor
        # m = where(m.mask, _FillValue, m)

        if nc_data_type == 'u1':
            m = uint8(m)
        elif nc_data_type == 'u2':
            m = uint16(m)
        elif nc_data_type == 'u4':
            m = uint32(m)
        elif nc_data_type == 'u8':
            m = uint64(m)

        # print '\n', _min, _max, '\n'
        # print m
        # print '\n', scale_factor, add_offset, _FillValue, '\n'
        # print m.astype(float) * scale_factor + add_offset

        sourcefield = model.get_field(fieldname)

        field = Field(copy(sourcefield.variable),
                      copy(sourcefield.dimensions),
                      values=m,
                      units=copy(sourcefield.units),
                      attributes=copy(sourcefield.attributes),
                      fillvalue=_FillValue,
                      valid_min=(sourcefield.valid_min -
                                 add_offset) / scale_factor,
                      valid_max=(sourcefield.valid_max -
                                 add_offset) / scale_factor,
                      )
        # outparameter = model.extract_subset(fields=[])
        # outparameter.add_field(field)

        field.scale_factor = scale_factor
        field.add_offset = add_offset

        export2file(url, 'NETCDF4', field, attrs, append)
        

def export2file(url,
                file_format,
                field,
                attrs={'acknowledgement': 'MAGiCAN'},
                append=False):
    """
    """

    if 'NETCDF4' in file_format:
        if not url.endswith('.nc'):
            url = url + '.nc'
        nc_handler = netcdf.open_nc(url, file_format, append)
        netcdf.create_nc_var(nc_handler, field)
        netcdf.write_nc_var(nc_handler, field)
        netcdf.write_global_attributes(nc_handler, attrs, append)
        nc_handler.close()

    return NotImplementedError
