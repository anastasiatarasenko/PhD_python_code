# -*- coding: utf-8 -*-
"""
.. module::posada.Extractors.netcdf

Class to write netCDF files using CF conventions.

:copyright: Copyright 2011-2016 SOLab.

:copyright: Copyright 2013 Ifremer / Cersat / SOLab.
:license: Released under GPL v3 license, see :ref:'license'.

.. sectionauthor:: Jeff Piolle <jfpiolle@ifremer.fr>, Alexander Myasoedov <mag@rshu.ru>
.. codeauthor:: Jeff Piolle <jfpiolle@ifremer.fr>, Alexander Myasoedov <mag@rshu.ru>
.. author:: SOLab Dev Team
"""

import os
import numpy
import datetime

from netCDF4 import Dataset as ncDataset
from netCDF4 import default_fillvals, getlibversion

from ..Readers.BaseReader import CorruptFileException
from ..Models.Field import Field
from ..Models.Variable import Variable

from ..Common.imports import OrderedDict, asarray
from copy import copy

from ..Common.logger import get_logger
logger = get_logger()

CONVENTIONS = 'CF-1.6, Unidata Observation Dataset v1.0'

STANDARD_ATTRIBUTES_VALUES = OrderedDict([
    ('id', ''),
    ('naming_authority', 'solab.rshu.ru'),
    ('Metadata_Conventions', 'Unidata Dataset Discovery v1.0'),
    ('standard_name_vocabulary',
     'NetCDF Climate and Forecast (CF) Metadata Convention'),
    ('institution', 'Russian State '
        'Hydrometeorological University/Satellite Oceanography Laboratory'),
    ('institution_abbreviation', 'rshu/solab'),
    ('title', ''),
    ('summary', ''),
    ('cdm_data_type', ''),
    ('keywords', ''),
    ('keywords_vocabulary', 'NASA Global Change Master Directory (GCMD)'
        ' Science Keywords'),
    ('standard_name_vocabulary', 'NetCDF Climate and Forecast (CF) '
        'Metadata Convention'),
    ('scientific_project', ''),
    ('acknowledgement', ''),
    ('license', ''),
    ('format_version', ''),
    ('processing_software', 'SOLab/POSADA'),
    ('product_version', ''),
    ('uuid', ''),
    ('processing_level', ''),
    ('history', ''),
    ('publisher_name', 'rshu/solab'),
    ('publisher_url', 'http://solab.rshu.ru'),
    ('publisher_email', 'sol@rshu.ru'),
    ('creator_name', ''),
    ('creator_url', ''),
    ('creator_email', ''),
    ('references', ''),
    ('metadata_link', ''),
    ('source', ''),
    ('source_version', ''),
    ('platform', ''),
    ('platform_type', ''),
    ('sensor', ''),
    ('sensor_type', ''),
    ('band', ''),
    ('spatial_resolution', ''),
    ('geospatial_lat_min', ''),
    ('geospatial_lat_max', ''),
    ('geospatial_lat_units', 'degrees'),
    ('geospatial_lat_resolution', ''),
    ('geospatial_lon_min', ''),
    ('geospatial_lon_max', ''),
    ('geospatial_lon_units', 'degrees'),
    ('geospatial_lon_resolution', ''),
    ('geospatial_vertical_min', ''),
    ('geospatial_vertical_max', ''),
    ('geospatial_vertical_units', 'meters above mean sea level'),
    ('geospatial_vertical_positive', 'up'),
    ('time_coverage_start', ''),
    ('time_coverage_end', ''),
    ('time_coverage_resolution', '')
])

TIME_FORMAT = '%Y%m%dT%H%M%SZ'


def open_nc(url, nc_format='NETCDF4', append=True):
    """
    Open the netCDF file for writing

    Args:
        url (str): full path to the NetCDF file
            where model fields will be saved.

        nc_format (str): NetCDF file format, could be one of
            'NETCDF3_CLASSIC' or 'NETCDF4'

        append (bool): if True then 'append', else - 'overwrite'.
            When 'append' mode is selected field(s) will be appended
            to the existing NetCDF file.
    Returns:
        a handler of the opened file
    """
    logger.debug("opening NetCDF file for writing")
    logger.debug("url : %s", url)

    # netCDF4 fails if filename is longer than 256 char
    chdir = False
    if len(url) >= 256:
        chdir = True
    try:
        # move to file directory if filename is too long
        if chdir:
            localpath = os.path.curdir
            os.chdir(os.path.dirname(url))
            url = os.path.basename(url)

        if os.path.isfile(url) and append:
            nc_handler = ncDataset(url, 'a', format=nc_format)
        elif os.path.isfile(url) and not append:
            try:
                os.remove(url)
            except:
                pass
            nc_handler = ncDataset(url, 'w', format=nc_format)
        else:
            nc_handler = ncDataset(url, 'w', format=nc_format)

        # move back to original folder
        if chdir:
            os.chdir(localpath)

        return nc_handler
    except:
        logger.error("Could not open file %s", url)
        raise CorruptFileException("Could not open file %s", url)


def write_global_attributes(nc_handler, attrs, append=False):
    """Write the global attributes of the file.

    Includes a set of predefined standard attributes for CF compliant
    files. Default attributes will be filled automatically with default
    values unless you override them in 'attrs' argument. You can also
    add additional attributes in 'attrs' dictionary.

    Args:
        attrs (dict<string, string or number or datetime>): a dictionary
            containing the attributes names and values to be written.
    """

    attributes = copy(attrs)
    if not append:
        nc_handler.Conventions = CONVENTIONS
        nc_handler.netcdf_version_id = getlibversion()
        nc_handler.date_created\
            = datetime.datetime.now().strftime(TIME_FORMAT)
        nc_handler.date_modified\
            = datetime.datetime.now().strftime(TIME_FORMAT)
        # convert datetime objects to string
        for att in attributes:
            if isinstance(attributes[att], datetime.datetime):
                attributes[att]\
                    = attributes[att].strftime(TIME_FORMAT)
        # we try to keep some order in the attributes
        for att in STANDARD_ATTRIBUTES_VALUES:
            if att in attributes:
                nc_handler.setncattr(att, attributes[att])
            else:
                if STANDARD_ATTRIBUTES_VALUES[att] != '':  # if dict not empty
                    nc_handler.setncattr(att,
                                         STANDARD_ATTRIBUTES_VALUES[att])
        # additional attributes
        #keys = attributes.keys()
        keys = list(attributes.keys())
        keys.sort()
        for att in keys:
            if att not in STANDARD_ATTRIBUTES_VALUES:
                nc_handler.setncattr(att, attributes[att])
    else:
        nc_handler.date_modified\
            = datetime.datetime.now().strftime('%Y%m%dT%H%M%SZ')
        for att in attributes:
            nc_handler.setncattr(att, attributes[att])
    return


def get_netcdf_fillvalue(datatype):
    if datatype == str:
        return str('')
    elif datatype == numpy.dtype('S1'):
        return ' '
    elif datatype == numpy.dtype(numpy.complex64):
        # take type of real and imag components
        comptype = numpy.dtype('float32')
        strtype = '%s%s' % (comptype.kind, comptype.itemsize)
    elif datatype == numpy.dtype(numpy.complex128):
        # take type of real and imag components
        comptype = numpy.dtype('float64')
        strtype = '%s%s' % (comptype.kind, comptype.itemsize)
    else:
        strtype = '%s%s' % (datatype.kind, datatype.itemsize)
    return default_fillvals[strtype]


def _create_nc_var(nc_handler,
                   name,
                   datatype,
                   dims,
                   comp,
                   field,
                   fillvalue,
                   scale_factor,
                   add_offset,
                   ):
    dropAttributes = ['dtype']

    if nc_handler.file_format != 'NETCDF3_CLASSIC':
        try:
            if (nc_handler.file_format != "NETCDF4" and
                    datatype == numpy.dtype(numpy.int64)):
                datatype = numpy.dtype(numpy.int32)
                if fillvalue is None:
                    fillvalue = get_netcdf_fillvalue(datatype)
                else:
                    fillvalue = numpy.int32(fillvalue)

            if nc_handler.file_format == 'NETCDF4':
                # cast from int32 to int64 or the file appeared corrupted
                # later (bug in netcdf library?).
                if datatype == 'int32':
                    datatype = numpy.dtype('int64')
                    fillvalue = get_netcdf_fillvalue(datatype)
                    fillvalue = numpy.int64(fillvalue)
                else:
                    fillvalue = get_netcdf_fillvalue(datatype)

            if datatype == str:
                fillvalue = None

            # assume for numpy object array that it contains strings
            # netCDF4 does not accept objet type => cast to string
            dtype = datatype
            if dtype == object:
                dtype = str

            for dim in dims.items():
                if dim[0] in nc_handler.dimensions:
                    pass
                else:
                    nc_handler.createDimension(dim[0], dim[1])

            result = nc_handler.createVariable(
                name,
                datatype=dtype,
                dimensions=list(dims),
                zlib=True,
                complevel=6,
                shuffle=True,
                fletcher32=False,
                fill_value=fillvalue)

        except:
            logger.error("Could not create field %s", name)
            raise
    else:
        result = nc_handler.createVariable(
            name,
            datatype=datatype,
            dimensions=dims.keys(),
            fill_value=fillvalue
        )

    # add variable attributes
    if comp.variable.description is not None:
        result.long_name = comp.variable.description
    if comp.variable.standardname is not None:
        result.standard_name = comp.variable.standardname
    if comp.units is not None:
        result.units = comp.units
    if comp.valid_min is not None:
        result.valid_min = asarray(comp.valid_min, dtype=dtype)
    if comp.valid_max is not None:
        result.valid_max = asarray(comp.valid_max, dtype=dtype)

    result.scale_factor = scale_factor
    result.add_offset = add_offset

    for att in comp.attributes:
        # Stop transferring some attributes which are redundant or
        # incorrect.
        if att.lower() not in dropAttributes:
            result.setncattr(att, comp.attributes[att])
    # create quality fields if defined
    if field.qc_levels is not None:
        result = nc_handler.createVariable(
            comp.variable.shortname + '_qc_level',
            'i1',
            dimensions=tuple(field.qc_levels.dimensions.keys()),
            fill_value=-127
        )
        result.flag_values = field.qc_levels.levels
        result.flag_meanings = field.qc_levels.meanings
    if field.qc_details is not None:
        result = nc_handler.createVariable(
            comp.variable.shortname + '_qc_details',
            'i2',
            dimensions=tuple(field.qc_details.dimensions.keys()),
            fill_value=-32767
        )
        result.flag_masks = field.qc_details.mask
        result.flag_meanings = field.qc_details.meanings


def _create_complex(nc_handler,
                    varName,
                    dType,
                    dims,
                    comp,
                    field,
                    fillValue, scale_factor, add_offset):
    """Internal wrapper for the creation of a complex field."""
    if fillValue is None:
        fillValue = get_netcdf_fillvalue(dType)
    names = (varName + '_real', varName + '_imag')
    for name in names:
        _create_nc_var(nc_handler, name, dType, dims, fillValue, comp, field)


def create_nc_var(nc_handler, field, dim_translation=None, model_type=None):
    """
    Create netCDF variable for provided field.

    Args:
        field (:class:'~posada.Models.Field.Field'): field to create.

        dim_translation (dict, optional): dictionary, alternative dimension
             names to be used.

        model_type (string, optional): class name of the
            :class:'~posada.Models'. Used to implement specific writing
            conventions for some :class:'~posada.Models'
            (ex: additional dummy time dimension for grids).
    .. Warning::
        Actual writing of data to the file is done by calling the function
        :func:'create_nc_var'.

    """
    try:
        varname = field.variable.shortname

        for comp in field.get_components():
            datatype = comp.datatype
            dims = comp.dimensions
            if dim_translation:
                # renaming of default dimensions
                dims = list(comp.dimensions)
                for idim, dim in enumerate(dims):
                    if dim in dim_translation:
                        dims[idim] = dim_translation[dim]

            # time variable
            if varname == 'time':
                # Need to convert datetime object to int32
                datatype = comp.datatype
                if datatype == datetime.datetime:
                    datatype = numpy.dtype(numpy.int32)
                    logger.warning("Datetime converted to int32")
                fillvalue = get_netcdf_fillvalue(datatype)
                scale_factor = float(1)
                add_offset = float(0.)

                # if time is two-dimensional (ex: L3 with regridded passes)
                # then the time coordinate variable can not have multiple
                # dimensions. We need two variables.
                if model_type == 'Grid' and len(field.dimensions) > 1:
                    # single dimension time variable. We use start sensing
                    # time as reference, by convention.
                    _create_nc_var(nc_handler, varname, datatype, ('time',),
                                   comp, field,
                                   fillvalue,
                                   scale_factor,
                                   add_offset)
                    # and an offset 'dtime' variable, in seconds
                    var = Variable(
                        'dtime',
                        description='time difference from reference time',
                    )
                    comment = 'time plus sst_dtime gives pixel time'
                    dims = field.dimensions
                    if 'time' not in dims:
                        dims = collections.OrderedDict([('time', 1)])
                        for dim in field.dimensions:
                            dims[dim] = field.dimensions[dim]
                    units = field.units.split('since')[0].strip()
                    dtimefield = Field(
                        var,
                        dimensions=dims,
                        units=units,
                        attributes={'comment': comment}
                    )
                    _create_nc_var(
                        nc_handler,
                        'dtime', datatype,
                        ('time', 'lat', 'lon',),
                        dtimefield, dtimefield,
                        get_netcdf_fillvalue(datatype),
                        scale_factor,
                        add_offset)
                    return

            # add a dummy time dimension for grids (CF convention) to all
            # fields except lat/lon
            atemporal_fields = ['lat', 'lon', 'lat_bnds', 'lon_bnds']
            if (model_type == 'Grid' and
                    varname not in atemporal_fields and
                    'time' not in dims):
                dims.insert(0, 'time')

            fillvalue = None
            scale_factor = float(1.)
            add_offset = float(0.)

            if comp.fillvalue is not None:
                fillvalue = comp.fillvalue
            if hasattr(comp, 'scale_factor'):
                if comp.scale_factor is not None:
                    scale_factor = comp.scale_factor
            if hasattr(comp, 'add_offset'):
                if comp.add_offset is not None:
                    add_offset = comp.add_offset

            if comp.datatype == numpy.dtype(numpy.complex128):
                # netcdf can not store complex and the field has to be
                # split into two netcdf variables
                _create_complex(nc_handler, comp.variable.shortname,
                                numpy.dtype(numpy.float64),
                                dims, comp, field,
                                fillvalue,
                                scale_factor,
                                add_offset)
            elif comp.datatype == numpy.dtype(numpy.complex64):
                # netcdf can not store complex and the field has to be
                # split into two netcdf variables
                _create_complex(nc_handler, comp.variable.shortname,
                                numpy.dtype(numpy.float32),
                                dims, comp, field,
                                fillvalue,
                                scale_factor,
                                add_offset)
            elif comp.datatype == numpy.str_:
                _create_nc_var(nc_handler, varname, 'S1', dims,
                               comp, field,
                               fillvalue,
                               scale_factor,
                               add_offset)
            else:
                # print(nc_handler)
                # print(varname)
                # print(datatype)
                # print(dims)
                # print(fillvalue)
                # print(comp)
                # print(field)
                _create_nc_var(nc_handler, varname, datatype, dims,
                               comp, field,
                               fillvalue,
                               scale_factor,
                               add_offset)
        return True
    except:
        logger.error("Error when creating field %s",
                     field.variable.shortname)
        raise


def write_nc_var(nc_handler, field):
    """
    Write the data of a field to the netCDF file.

    Args:
        field (:class:'~posada.Models.Field.Field'): the field which
            contents has to be written to the file.

    .. Warning::
       The netcdf variable with all its metadata (attributes) must have been
       explicitely created by the :func:'create_nc_var' before,
       by calling the :func:'create_nc_var' function.
    """
    try:
        values = field.get_values()
        # print values

        if values is None:
            raise Exception(
                "No values provided for field ", field.name)

        # forbids int64
        if nc_handler.file_format != "NETCDF4":
            if values.dtype == numpy.dtype(numpy.int64):
                logger.warning(
                    "Forced type of data from int64 to int32")
                values = values.astype(numpy.int32)

        # Need to convert incoming datetime object to a floating point
        # time as netcdf cannot handle time properly.
        if field.name == 'time':
            if isinstance(values.reshape(-1)[0], datetime.datetime):
                units = DEFAULT_TIME_UNITS
                if field.units is not None:
                    units = field.units
                values = netCDF4.date2num(values, units)
            # Else already in proper units (we hope).

        logger.debug("Writing field %s", field.name)

        # netcdf can not store complex and the field is
        # splitted in two netcdf variables
        if field.datatype == numpy.dtype(numpy.complex64)\
                or field.datatype == numpy.dtype(numpy.complex128):
            ncvar_real = nc_handler.variables[
                field.name + '_real']
            ncvar_imag = nc_handler.variables[
                field.name + '_imag']
            if ncvar_real.shape == (0,)\
                    or values.shape == ()\
                    or ncvar_real.shape == values.shape:
                ncvar_real[:] = values.real
                ncvar_imag[:] = values.imag
            else:
                raise IndexError(
                    str(
                        "Your input data shape {} is not compatible with "
                        "your field shape {}."
                    ).format(
                        values.shape, ncvar_real.shape
                    )
                )
        # for Grid objet with 2-dimensional time fields, two variables have
        # to be written : time and dtime
        elif (field.name == 'time' and
              'dtime' in nc_handler.variables):
            ncvar = nc_handler.variables['time']
            # start time is set as reference by convention
            ncvar[0] = values.min()
            ncvar = nc_handler.variables['dtime']
            ncvar[0, ...] = values - values.min()
        else:
            ncvar = nc_handler.variables[field.name]
            if ncvar.shape == (0,)\
                    or ncvar.shape == ()\
                    or values.shape == ()\
                    or ncvar.shape == values.shape\
                    or (ncvar.shape == (1,) and values.size == 1):
                ncvar[:] = values[:]

            elif (len(ncvar.shape) - 1 == len(values.shape) and
                  (ncvar.shape[0] == 1 or ncvar.shape[0] == 0)):
                # case of grids
                ncvar[0, ...] = values

            else:
                raise IndexError(
                    str(
                        "Your input data shape {} is not compatible with"
                        " your field shape {}."
                    ).format(
                        values.shape, ncvar.shape
                    )
                )

        if field.qc_levels is not None:
            nc_handler.variables[field.name + '_qc_level'][:]\
                = field.qc_levels.values[:]

        if field.qc_details is not None:
            fieldname = field.name + '_qc_details'
            nc_handler.variables[
                fieldname].set_auto_maskandscale(True)
            nc_handler.variables[fieldname][:]\
                = field.qc_details.values[:]

    except:
        logger.error("Could not write field {}".format(field.name))
        raise
