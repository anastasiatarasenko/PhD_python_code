# coding: utf-8

"""
.. module::posada.Extractors.resampler

A module for resampling/reprojecting POSADA Products/Parameters.

:copyright: Copyright 2013 Ifremer / Cersat / SOLab.
:license: Released under GPL v3 license, see :ref:'license'.

.. sectionauthor:: Jeff Piolle <jfpiolle@ifremer.fr>, Alexander Myasoedov <mag@rshu.ru>
.. codeauthor:: Jeff Piolle <jfpiolle@ifremer.fr>, Alexander Myasoedov <mag@rshu.ru>
.. author:: SOLab Dev Team
"""

from copy import copy
import numpy
from collections import OrderedDict
from pyresample import kd_tree, geometry, utils
from scipy import interpolate, ndimage
from sys import version_info
#from scipy.stats import nanmean

from ..Common.imports import array, size, pow,\
    pr, SpatialReference, Proj, cpu_count

from ..Models.Field import Field
from ..Models.Variable import Variable
from ..Models.swath import Swath
from ..Models.grid import Grid
from ..Models.satpass import Pass
# from ..Models.trajectory import Trajectory

from ..Common.logger import get_logger
logger = get_logger()

__all__ = ['Resampler']

try:
    from .fastregrid import fastregrid
except ImportError:
    logger.warning(
        "This POSADA instance doesn't include a compiled Cython "
        "module 'pyx'. This is probably because you are running "
        "cerbere from a repository, and not as an installed package. "
        "Building 'pyx' implace...")

if version_info[0] == 3:
    pass
else:
    import pyximport
    pyximport.install(setup_args={'include_dirs': numpy.get_include()})
    from .fastregrid import fastregrid
    logger.warning("'pyx' is available but built inplace.")


def swath_area_def(cs_str,
                   lonlim=(-180, 180), latlim=(-90, 90),
                   res=111.2e3,
                   area_extent=None):
    """
    Convert given swath coordinates to pyresample area definition.

    Args:
        cs_str (string): coordinate system name in form 'EPSG:3575' on which
                to resample the source parameter.

    """

    if size(res) == 1:
        res = array((res, res))

    up = max(latlim)
    down = min(latlim)
    left = min(lonlim)
    right = max(lonlim)

    srs = SpatialReference()
    proj4_args = srs.SetFromUserInput(cs_str)
    proj4_args = srs.ExportToProj4()
    proj4_args = proj4_args.replace(
        ' +no_defs', '')  # no need for '+no_defs'

    area_name = srs.GetAttrValue('projcs')
    if not area_name:
        area_name = cs_str

    area_id = area_name.replace(
        " ", "_").replace(
        "_/_", "_").replace(
        ":", "_").lower()
    proj_id = area_id

    p = Proj(proj4_args)

    left_ex1, up_ex1 = p(left, up)
    right_ex1, up_ex2 = p(right, up)
    left_ex2, down_ex1 = p(left, down)
    right_ex2, down_ex2 = p(right, down)

    if not area_extent:
        if 'stere' in proj4_args:
            lon = (left + right) / 2.0
            if (lon >= 0 and lon < 90) or (lon >= -360 and lon < -270):
                logger.debug("11111111111")
                area_extent = array(
                    (
                        min(left_ex1, left_ex2, right_ex1, right_ex2),
                        min(down_ex1, down_ex2, up_ex1, up_ex2),
                        max(left_ex1, left_ex2, right_ex1, right_ex2),
                        max(down_ex1, down_ex2, up_ex1, up_ex2)
                    )
                )
            elif (lon >= 90 and lon < 180) or (lon >= -270 and lon < -180):
                logger.debug("2222222222222")
                area_extent = array(
                    (
                        max(left_ex1, left_ex2, right_ex1, right_ex2),
                        max(down_ex1, down_ex2, up_ex1, up_ex2),
                        min(left_ex1, left_ex2, right_ex1, right_ex2),
                        min(down_ex1, down_ex2, up_ex1, up_ex2)
                    )
                )
            elif (lon >= 180 and lon < 270) or (lon >= -180 and lon < -90):
                logger.debug("333333333333")
                area_extent = array(
                    (
                        min(left_ex1, left_ex2, right_ex1, right_ex2),
                        min(down_ex1, down_ex2, up_ex1, up_ex2),
                        max(left_ex1, left_ex2, right_ex1, right_ex2),
                        max(down_ex1, down_ex2, up_ex1, up_ex2)
                    )
                )
            else:
                logger.debug("44444444444444444")
                area_extent = array(
                    (
                        min(left_ex1, left_ex2, right_ex1, right_ex2),
                        min(down_ex1, down_ex2, up_ex1, up_ex2),
                        max(left_ex1, left_ex2, right_ex1, right_ex2),
                        max(down_ex1, down_ex2, up_ex1, up_ex2)
                    )
                )
        else:
            area_extent = array(
                (
                    min(left_ex1, left_ex2, right_ex1, right_ex2),
                    min(up_ex1, up_ex2, down_ex1, down_ex2),
                    max(left_ex1, left_ex2, right_ex1, right_ex2),
                    max(up_ex1, up_ex2, down_ex1, down_ex2)
                )
            )

        modulox = (area_extent[2] - area_extent[0]) % res[0]
        moduloy = (area_extent[3] - area_extent[1]) % res[1]
        area_extent[0] = area_extent[0] + modulox
        area_extent[1] = area_extent[1] + moduloy
        area_extent[2] = area_extent[2]
        area_extent[3] = area_extent[3]

        # Using abs() to avoid negative numbers of
        # coloumns/rows as for epsg3413 for example
        xsize = abs(int(round((area_extent[2] - area_extent[0]) / res[0])))
        ysize = abs(int(round((area_extent[3] - area_extent[1]) / res[1])))
    elif area_extent:
        # get list of available resolutions
        avail_resolution_list = [int(256 * pow(2, i)) for i in range(15)]

        # get first guess size of area
        # size = int(256*math.pow(2, 14)/res)
        _size = abs(int(area_extent[0] - area_extent[2]) / min(res))

        # get closest size of area from available resolutions list
        xsize = min(filter(lambda x: _size <= x,
                           avail_resolution_list))
        ysize = xsize

    area_def = pr.utils.get_area_def(area_id, area_name,
                                     proj_id, proj4_args,
                                     xsize, ysize, area_extent)

    return area_def


class Resampler(object):
    """
    Toolbox for resampling/reprojecting of POSADA Products/Parameters
    on top of each others or on any other grid/projection.
    Providing several methods for resampling/interpolating.
    """

    def __init__(self):
        pass

    @classmethod
    def __monotonic(cls, longitudes):
        values = longitudes
        """make longitude values monotonically increasing or decreasing"""
        if len(values.shape) == 1:
            diff = values[1:] - values[0:-1]
        else:
            # take some values along x axis
            # not bullet proof... to be improved
            diff = values[0, 1:] - values[0, 0:-1]
        if diff.min() * diff.max() < 0:
            logger.warning("Trying to make longitudes monotonic")
            values = longitudes.copy()
            values[values < 0] += 360.
        return values

    @classmethod
    def extrapolate(cls, data, extrapolation=4, copy=False):
        """Fill in missing values with a mean of the neighbours.

        Args:
            data (:class:'numpy.ma.MaskedArray'): array to fill in.
            extrapolation (optional, int): number of iterations of the
                extrapolation process (how far to extrapolate).
            copy (optional, boolean): if True, copy the result in a new array.
        """
        # Get array filled with NaN and mask
        # NOTE : in a global case, it is possible here to duplicate the fist/last
        # columns at the end/beginning then remove it before returning
        array = data.filled(numpy.NaN)
        if copy:
            mask = numpy.ma.getmaskarray(data).copy()
        else:
            mask = numpy.ma.getmaskarray(data)

        for dummy in range(extrapolation):

            # Find points to inpaint
            eroded_mask = ndimage.morphology.binary_erosion(
                mask, border_value=1)
            indexinp = numpy.where(mask & ~eroded_mask)

            # Find neighbours
            xneighb = numpy.repeat(indexinp[1], 8).reshape((-1, 8)) + \
                numpy.array([-1, 0, 1, 1, 1, 0, -1, -1])
            yneighb = numpy.repeat(indexinp[0], 8).reshape((-1, 8)) + \
                numpy.array([-1, -1, -1, 0, 1, 1, 1, 0])

            # Avoid neighbours to be out of bounds
            # V1 : clip (some neighbours will be duplicated)
            # numpy.clip(xneighb, 0, array.shape[1]-1, out=xneighb)
            # numpy.clip(yneighb, 0, array.shape[0]-1, out=yneighb)
            # V2 : replace with the point to inpaint which has to be NaN
            indexout = numpy.where((xneighb < 0) | (xneighb > array.shape[1] - 1) |
                                   (yneighb < 0) | (yneighb > array.shape[0] - 1))
            xneighb[indexout] = indexinp[1][indexout[0]]
            yneighb[indexout] = indexinp[0][indexout[0]]

            # Replace with neighbours mean and update mask
            try:
                newpts = numpy.nanmean(array[[yneighb, xneighb]], axis=1)
            except RuntimeWarning:
                raise Exception('Unexpectedly, some points have all neighbours '
                                'set to NaN')
            array[indexinp] = newpts
            mask[indexinp] = False

        data = numpy.ma.MaskedArray(array, mask, copy=False)
        return data

    @classmethod
    def interp2highResGrid(cls,
                           source,
                           target,
                           fields=None,
                           extrapolation=0,
                           fillvalue=None,
                           new=False,
                           add_reference=False,
                           add_geolocation=True,
                           mask=None,
                           prefix='interpolated_'
                           ):
        """
        Use interp2d to interpolate source grid over a target high resolution
        grid. On contrary to fastReGrid.

        Usually, the source grid should be lower or equal resolution than the
        target grid.

        Note that the interpolation will not perform well if the longitudes are
        not monotonically ordered or in the same range of values between the
        source and the target. For instance, wrong interpolation will happen
        if:
          * ''source'' is a global grid starting at 0 with longitudes expressed
            between -180/180 (not monotone because of meridian 180). It is
            recommended to shift the source by 180 (cerbere convention is to
            start global grids at -180.)
          * ''source'' is expressed between -180/180 and target between 0/360
            (cerbere convention is to express longitudes between -180/180).

        Args:
            source (:class:'..Models.grid.Grid'): source grid.

            fieldname (string or list of string): name of the field(s) to
                interpolate.

            target (:class:'..Models.grid.Grid'): target grid.

            extrapolation (int): if different from 0, extrapolate values of
                source over edges (gaps, land, sea-ice,...). Default is 0. The
                provided value is the extent, in pixels, of the extrapolation
                from edges.

            fillvalue (float): value to replace the missing values. The
                interpolation considers the missing values as ordinary values.
                These can be replaced with a more realistic value using this
                argument (for instance water freezing temperature when
                interpolating SST along sea ice edges).

            mask (ndarray):  a boolean mask to apply on the interpolated
                values after they have been calculated. The mask must match the
                dimensions of the target grid. This is meant for instance to
                mask land data over which interpolation may have been performed
                in the resampling.
        """
        if fields is None:
            logger.warning("No field names were provided")
            return None
        if type(fields) is not list:
            fields = [fields]
        if new:
            # create a clone of target without any field
            outparameter = target.extract_subset(fields=[])
        else:
            outparameter = target
        target_lons = target.get_lon()
        target_lats = target.get_lat()
        source_lons = source.get_lon()
        source_lats = source.get_lat()
        # interpolate each field independently
        for fieldname in fields:
            # extrapolating low resolution analysis over edges
            # ------------------------------------------------
            data = source.get_values(fieldname)
            if fillvalue is None:
                missing = source.get_field(fieldname).fillvalue
            else:
                missing = fillvalue
            logger.debug('missing value: %s', missing)
            if mask is None:
                srcmask = numpy.array(data.mask, copy=True)
            if extrapolation != 0:
                logger.debug('Read source and extrapolate over missing data')
                # call extrapolation procedure (replace NaN with interpolated
                # values)
                logger.debug('data %s', data.shape)
                data = cls.extrapolate(data, extrapolation)
            # masked data to missing value
            data = data.filled(missing)
            # re-interpolate source over target
            # ---------------------------------
            source_lons_new = Resampler.__monotonic(source_lons)
            test_val = source_lons_new != numpy.unique(source_lons_new)
            # means that there is not unique values in source_lons_new
            if isinstance(test_val, bool):
                logger.warning(
                    'longitudes of the sources contains not unique values ( interpolation wont be bijective')
            finterp = interpolate.interp2d(
                source_lons_new,
                source_lats,
                data,
                kind='linear',
                bounds_error=False,
                fill_value=None)
            interpolated_data = finterp(
                Resampler.__monotonic(target_lons).ravel(),
                target_lats.ravel())
            # mask the result with provided mask, if any
            if mask is not None:
                interpolated_data = numpy.ma.masked_where(
                    mask,
                    interpolated_data,
                    copy=False
                )
            else:
                # extrapolate the source mask
                extrapolate_mask = Resampler.__closest_neighbour(source_lats,
                                                                 source_lons,
                                                                 srcmask.astype(
                                                                     'int8'),
                                                                 source.get_model_name(),
                                                                 target_lats,
                                                                 target_lons,
                                                                 target.get_model_name(),
                                                                 )
                interpolated_data = numpy.ma.masked_where(
                    extrapolate_mask == 1,
                    interpolated_data,
                    copy=False
                )
            #==================================================================
            # interpolated_data = numpy.ma.masked_values(
            #     interpolated_data,
            #     missing,
            #     copy=False)
            #==================================================================
            interpolated_data = numpy.ma.masked_where(
                numpy.isnan(interpolated_data.data),
                interpolated_data,
                copy=False
            )
            target_shape = target.get_geolocation_dimsizes().values()
            if isinstance(target, Grid):
                interpolated_data = interpolated_data.reshape(target_shape)
            if outparameter.has_field(fieldname):
                raise Exception('Field already existing in target grid')
            else:
                var = copy(source.get_field(fieldname).variable)
                sourcefield = source.get_field(fieldname)
                var.shortname = prefix + var.shortname
                new_dims = target.get_geolocation_dimsizes()
                field = Field(var,
                              new_dims,
                              values=interpolated_data,
                              units=sourcefield.units,
                              attributes=copy(sourcefield.attributes),
                              fillvalue=sourcefield.fillvalue,
                              valid_min=sourcefield.valid_min,
                              valid_max=sourcefield.valid_max,
                              )
                outparameter.add_field(field)

        # Add geolocation
        if add_geolocation:
            # add time and lat/lon
            for fieldname in ['time']:
                field = source.get_geolocation_field(fieldname).clone()
                if len(field.dimensions) == 1:
                    field.variable.shortname = 'source_time'
                    field.set_values(
                        source.get_geolocation_field(
                            fieldname).get_values().copy()
                    )
                    outparameter.add_field(field)
        return outparameter

    @classmethod
    def __get_pyresample_def(cls,
                             srclat,
                             srclon,
                             srcmodel,
                             tgtlat,
                             tgtlon,
                             tgtmodel,
                             radius=None,
                             ):
        """
        Returns the required parameters for pyresample closest neighbout
        resampling.

        Refer to:
        'pyresample <https://pyresample.readthedocs.org/en/latest/index.html>'_
        package documentation.
        """
        if (srclon > 180).any():
            srclon = utils.wrap_longitudes(srclon)
        if srcmodel == 'Swath' or srcmodel == 'Pass' or \
                srcmodel == 'Trajectory' or srcmodel == 'PointCollection':
            source_def = geometry.SwathDefinition(lons=srclon,
                                                  lats=srclat
                                                  )
        elif srcmodel == 'Grid':
            if len(srclon.shape) == 1:
                srclon, srclat = numpy.meshgrid(srclon, srclat)

            source_def = geometry.GridDefinition(lons=srclon,
                                                 lats=srclat
                                                 )
        if tgtmodel == 'Swath' or tgtmodel == 'Pass' or \
                tgtmodel == 'Trajectory':
            target_def = geometry.SwathDefinition(lons=tgtlon,
                                                  lats=tgtlat
                                                  )
        elif tgtmodel == 'Grid':
            if len(tgtlon.shape) == 1:
                tgtlon, tgtlat = numpy.meshgrid(tgtlon, tgtlat)
            target_def = geometry.GridDefinition(lons=tgtlon,
                                                 lats=tgtlat
                                                 )
        if radius is None:
            if srcmodel == 'PointCollection':
                #                 resolution = abs(srclat[0] - srclat[1])
                # resolution = numpy.amin(abs(numpy.diff(srclat)))
                resolution = 1
            else:
                y, x = srclat.shape
                resolution = abs(srclat[y / 2, x / 2] -
                                 srclat[y / 2 + 1, x / 2])
            radius = resolution * 100000. / 2. * 1.6
        logger.debug("Radius for resampling : %s", radius)
        return source_def, target_def, radius

    @classmethod
    def __get_closest_neighbour_info(cls,
                                     source_def,
                                     target_def,
                                     radius=None,
                                     ):
        """Returns the indexes for a closest neighbour resampling.

        This is a wrapper for :func:'get_neighbour_info' in the
        'pyresample <https://pyresample.readthedocs.org/en/latest/index.html>'_
        package.
        """
        valid_input_index, valid_output_index, index_array, distance_array\
            = kd_tree.get_neighbour_info(source_def,
                                         target_def,
                                         radius_of_influence=float(radius),
                                         reduce_data=False,
                                         neighbours=1
                                         )
        return (valid_input_index, valid_output_index,
                index_array, distance_array)

    @classmethod
    def __closest_neighbour(cls,
                            srclat,
                            srclon,
                            srcdata,
                            srcmodel,
                            tgtlat,
                            tgtlon,
                            tgtmodel,
                            radius=None,
                            ):
        """Resample an array using closest neighbour

        This is a wrapper for
        'pyresample <https://pyresample.readthedocs.org/en/latest/index.html>'_
        package.
        """
        source_def, target_def, radius = Resampler.__get_pyresample_def(
            srclat,
            srclon,
            srcmodel,
            tgtlat,
            tgtlon,
            tgtmodel,
            radius)
        valid_input_index, valid_output_index, index_array, distance_array\
            = Resampler.__get_closest_neighbour_info(source_def,
                                                     target_def,
                                                     radius=radius)
        if not valid_input_index.any():
            # return None if no source data can be resampled on the target
            return None
        logger.debug('target %s srcdtaa %s', target_def.shape, srcdata)
        res = kd_tree.get_sample_from_neighbour_info('nn',
                                                     target_def.shape,
                                                     srcdata,
                                                     valid_input_index,
                                                     valid_output_index,
                                                     index_array,
                                                     fill_value=None
                                                     )
        return res

    @classmethod
    def resampleModel2Model(cls,
                            source,
                            target,
                            fields,
                            radius=None,
                            new=False,
                            add_reference=False,
                            add_geolocation=True,
                            prefix='resampled_'
                            ):
        """
        Return new parameter corresponding to the resampled 'source'
        parameter on the 'target' parameter.

        Args:
            source (parameter class from 'Models' package): parameter to
                resample.

            target (parameter class from 'Models' package): parameter on which
                to resample the source parameter.

            fields (list): names of the fields from source to add to the
                returned parameter.

            new (bool): if new is True, creates a new parameter matching the
                'target' characteristics and add the resampled fields in this
                new parameter. if False, the resampled fields are added to
                'target' parameter.

            add_reference (bool): if True, add the original indice of the
                resampled data (only when nearest neighbour is used).

            add_geolocation (bool): if True, add the geolocation information
                from the resampled pixels.
        """
        if not type(fields) is list:
            fields = [fields]
        if fields == []:
            logger.warning("No fields requested? Skipping resampling.")
            return

        # Get source and target grid definitions
        source_lons = source.get_lon().astype('float32')
        source_lats = source.get_lat().astype('float32')
        source_times = source.get_times()
        if isinstance(source, Grid):
            source_lons, source_lats = numpy.meshgrid(source_lons, source_lats)
            source_times = numpy.resize(source_times, source_lons.shape)
        target_lons = target.get_lon().astype('float32')
        target_lats = target.get_lat().astype('float32')
        if isinstance(target, Grid) and len(target_lons.shape) == 1:
            target_lons, target_lats = numpy.meshgrid(target_lons, target_lats)
        source_def, target_def, radius = Resampler.__get_pyresample_def(
            source_lats,
            source_lons,
            source.__class__.__name__,
            target_lats,
            target_lons,
            target.__class__.__name__,
            radius)
        if new:
            # create a clone of target without any field
            outparameter = target.extract_subset(fields=[])
        else:
            outparameter = target
        first = True

        resampled_fields = []
        if add_geolocation:
            resampled_fields.extend(['lat', 'lon', 'time'])
        resampled_fields.extend(fields)
        for fieldname in resampled_fields:
            if first:
                valid_input_index, valid_output_index,\
                    index_array, distance_array\
                    = Resampler.__get_closest_neighbour_info(source_def,
                                                             target_def,
                                                             radius)
                first = False
                if not valid_input_index.any():
                    # return None if no source data can be resampled on the
                    # target
                    return None
            if fieldname == 'lat':
                values = source_lats
            elif fieldname == 'lon':
                values = source_lons
            elif fieldname == 'time':
                values = source_times
            else:
                values = source.get_values(fieldname)
            res = kd_tree.get_sample_from_neighbour_info(
                'nn',
                target_def.shape,
                values,
                valid_input_index,
                valid_output_index,
                index_array,
                fill_value=None
            )
            target_shape = target.get_geolocation_dimsizes().values()
            if isinstance(target, Grid):
                res = res.reshape(target_shape)
            if fieldname in ['time', 'lat', 'lon']:
                if not res.any():
                    # return None if no valid source data can be resampled on
                    #  the target
                    if not (values.min() == 0. and values.max() == 0):
                        return None

            newfieldname = prefix + fieldname
            if outparameter.has_field(newfieldname):
                raise Exception(
                    'Field {} already existing in target grid'.format(newfieldname))
            else:
                if fieldname in ['lat', 'lon', 'time']:
                    var = copy(source._geolocation_fields[
                        fieldname].variable)
                    sourcefield = source._geolocation_fields[fieldname]
                else:
                    var = copy(source.get_field(fieldname).variable)
                    sourcefield = source.get_field(fieldname)
                var.shortname = newfieldname
                new_dims = target.get_geolocation_dimsizes()
                field = Field(var,
                              new_dims,
                              values=res,
                              units=sourcefield.units,
                              attributes=copy(sourcefield.attributes),
                              fillvalue=sourcefield.fillvalue,
                              valid_min=sourcefield.valid_min,
                              valid_max=sourcefield.valid_max,
                              )
                outparameter.add_field(field)

        if add_reference:
            # Add the indices (in the source parameter data) of the
            # selected pixels
            target_dims = target.get_geolocation_dimsizes()
            if isinstance(source, Swath) or isinstance(source, Pass) or isinstance(source, Trajectory):
                xdimname = 'cell'
                ydimname = 'row'
            elif isinstance(source, Grid):
                xdimname = 'x'
                ydimname = 'y'
            if isinstance(source, Trajectory) == False:
                ny, nx = source_lons.shape
            else:
                ny = len(source_lons)
                nx = 1
            # the array_index indices correspond to a reduced version
            # of the data where invalid data have been reduced
            # (reducing over valid_input_index). We need to recover the
            # original indices in the full source data.
            true_index = numpy.arange(
                ny * nx, dtype='uint32')[valid_input_index]
            # mask where no data resampled (pyresample use the length of the
            # reduced data for this)
            index_array = numpy.ma.masked_equal(index_array.ravel(),
                                                true_index.size)
            index_array[valid_output_index] = numpy.ma.array(
                true_index[index_array[valid_output_index].filled(0)],
                mask=index_array[valid_output_index].mask
            )
            index_array = index_array.reshape(target_shape)
            reference_x = numpy.ma.array(index_array % nx)
            reference_y = numpy.ma.array(index_array / nx)
            varx = Variable(
                shortname='resampled_%s' % xdimname,
                description='%s index of the selected pixel' % xdimname
            )
            fieldx = Field(varx,
                           target_dims,
                           values=reference_x,
                           datatype=numpy.dtype(numpy.int32),
                           fillvalue=-2363645994848887)
            vary = Variable(
                shortname='resampled_%s' % ydimname,
                description='%s index of the selected pixel' % ydimname
            )
            fieldy = Field(vary,
                           target_dims,
                           values=reference_y,
                           datatype=numpy.dtype(numpy.int32),
                           fillvalue=-2363645994848887)
            outparameter.add_field(fieldx)
            outparameter.add_field(fieldy)
            # Add the distance
            vardist = Variable(
                shortname='distance_to_pixel',
                description='distance of the selected pixel to the pixel center'
            )
            fielddist = Field(
                vardist,
                target_dims,
                values=numpy.ma.array(distance_array.reshape(target_shape),
                                      mask=res.mask
                                      ),
                units='m',
                datatype=numpy.dtype(numpy.float32)
            )
            outparameter.add_field(fielddist)
        return outparameter

    @classmethod
    def fastReGrid(cls, src, target, fields, data=None,
                   neighbours=1,
                   use_src_time=False,
                   add_empty_fields=False):
        """
        Resample the source parameter to a regular grid, binning with
        arithmetic mean.
        Use fastReGrid to interpolate source grid over a target low resolution
        grid. On contrary to interp2highResGrid.

        .. warning:
           works only for now for fields with spatial and temporal dimensions
           only. This will be improved in a later version.

        .. warning:
           arithmetic mean is not possible on discrete values. Integer values
           will be converted to double then back to integer taking the floor
           value after regridding. Be extra cautious when providing integer
           values to this function.

        Args:
            src (Parameter): the parameter to be resampled over ''target''
                parameter.
            neighbours (int): the minimum number of values in a bin required to
                calculate a value (the bin is not filled otherwise).
            use_src_time (bool): use the source parameter as timestamp for the
                target parameter.
            add_empty_fields (bool): if True, always create a field even when
                no valid data were remapped. If false, the remapped field is
                not addedd to ''target'' when it is empty.

        Returns:
            Parameter: the parameter on which the remapped fields have
                been added. None if there was no valid data in the remapped
                parameter and 'add_empty_field' was set to False when calling the
                function.
        """
        logger.debug("Resample to regular grid")
        srclons = src.get_lon()
        srclats = src.get_lat()
        logger.debug('srclons: %s', srclons)
        logger.debug('srclats: %s', srclats)
        if not target.projection.is_cylindrical():
            raise Exception(
                'Projection is not valid for this type of regridding')
        tarlon = target.get_lon()
        tarlat = target.get_lat()
        lonbin = tarlon[1] - tarlon[0]
        latbin = tarlat[1] - tarlat[0]
        lonmin = tarlon[0] - lonbin / 2.
        lonmax = tarlon[-1] + lonbin / 2.
        if lonmin == lonmax:
            lonmax += 360.
        latmin = tarlat[0] - latbin / 2.
        latmax = tarlat[-1] + latbin / 2.
        if isinstance(src, Grid) and len(srclons.shape) == 1:
            srclons, srclats = numpy.meshgrid(srclons, srclats)
        # dimensions of a field in target grid
        geodims = OrderedDict([('lat', len(tarlat)),
                               ('lon', len(tarlon))])
        empty = True
        for i, fieldname in enumerate(fields):
            if data is None:
                val = src.get_values(fieldname)
            else:
                val = data[i]
            logger.debug('val %s %s', val.shape, val.mean())
            lons = numpy.ma.masked_where(val.mask, srclons).compressed()
            lats = numpy.ma.masked_where(val.mask, srclats).compressed()
            # create empty field if add_empty_fields is set when there are
            # no valid values to remap.
            if len(lats) == 0 and add_empty_fields:
                dims = (len(tarlat), len(tarlon))
                result = numpy.ma.masked_all(shape=dims,
                                             dtype=val.dtype)
            # there are valid values.Regrid with arithmetic mean.
            elif len(lats) != 0:
                empty = False
                # if values are integer, they are converted to float
                # arithmetic mean is not possible on discrete values, be extra
                # cautious when providing integer values.
                valdtype = val.dtype
                if val.dtype not in ['float32', 'float64']:
                    logger.warning("Converting to double this field: %s" %
                                   fieldname)
                    val = val.astype('float64')
                logger.debug('latmin %s latmax %s lonmin %s lonmax %s',
                             latmin, latmax, lonmin, lonmax)
                val_comp = val.compressed()
                logger.debug('shape val_compressed: %s', val_comp.shape)
                logger.debug('lats compressed :%s', lats.shape)
                result, nbsampl = fastregrid(
                    lats,
                    lons,
                    val_comp,
                    latmin, latmax, lonmin, lonmax,
                    lonbin, latbin
                )
                nbsampl = numpy.ma.masked_where(nbsampl < neighbours,
                                                nbsampl,
                                                copy=False)
                result = numpy.ma.masked_where(nbsampl.mask,
                                               result,
                                               copy=False)
                # only arithmetic mean is supported for now
                result = result / nbsampl.filled(1)
                # convert back to initial type
                if result.dtype != valdtype:
                    if valdtype not in ['float32', 'float64']:
                        result = numpy.ma.floor(result).astype(valdtype)
                    else:
                        result = result.astype(valdtype)
            # they are valid data in remapped field (or we force returning
            # an empty field)
            if len(lats) != 0 or add_empty_fields:
                if target.has_field(fieldname):
                    raise Exception('Variable %s already existing in target'
                                    ' grid' % fieldname)
                else:
                    logger.debug("create target grid field for %s"
                                 % fieldname)
                    if src.has_field(fieldname):
                        field = src.get_field(fieldname).clone()
                        field.dimensions = geodims
                        field.set_values(result)
                    else:
                        var = Variable(shortname=fieldname)
                        field = Field(
                            var,
                            dimensions=geodims,
                            values=result
                        )
                    field.attributes['binning_method'] = 'mean'
                    target.add_field(field)
        # return None if there is valid remapped data
        if empty:
            logger.warning("There were no valid data in any remapped fields."
                           " No field was returned.")
            return None
        if use_src_time:
            times = src.get_times()
            if times.shape == (1,):
                target.set_times(times, units=src.get_time_units())
            else:
                valdtype = times.dtype
                if valdtype not in ['float32', 'float64']:
                    logger.warning("Converting to double this field: times")
                    times = times.astype('float64')
                lons = numpy.ma.masked_where(times.mask, srclons).compressed()
                lats = numpy.ma.masked_where(times.mask, srclats).compressed()
                result, nbsampl = fastregrid(
                    lats,
                    lons,
                    times.compressed(),
                    latmin, latmax, lonmin, lonmax,
                    lonbin, latbin
                )
                nbsampl = numpy.ma.masked_where(nbsampl < neighbours,
                                                nbsampl,
                                                copy=False)
                result = numpy.ma.masked_where(nbsampl.mask,
                                               result,
                                               copy=False)
                times = result / nbsampl.filled(1)
                # convert back to initial (integer) type
                if result.dtype != valdtype:
                    times = numpy.ma.floor(times).astype(valdtype)
                timefield = target.get_geolocation_field('time')
                if timefield is None:
                    timefield = Field(
                        variable=copy(
                            src.get_geolocation_field('time').variable),
                        dimensions=geodims,
                        units=src.get_time_units(),
                        values=times
                    )
                    target._geolocation_fields['time'] = timefield
        return target

    @classmethod
    def resampleModel2proj(cls,
                           source,
                           target_area_def,
                           fields,
                           radius=None,
                           add_reference=False,
                           add_geolocation=True,
                           prefix='resampled_'
                           ):
        """
        Return new parameter corresponding to the resampled 'source'
        parameter on the 'target' parameter.

        Args:
            source (parameter class from 'Models' package): parameter to
                resample.

            target_area_def (object of type
                'pyresample.geometry.AreaDefinition'): The cartographic
                definition of grid or swath areas used by Pyresample.

            fields (list): names of the fields from source to add to the
                returned parameter.

            radius (number): the radius around each grid pixel in meters
                to search for neighbours in the swath.

            add_geolocation (bool): if True, add the geolocation information
                from the resampled pixels.

            prefix (string): prefix to the resampled field name/
        """

        # Create target model from target_area_def
        target_lons, target_lats = target_area_def.get_lonlats(
            nprocs=cpu_count())
        target = Swath(longitudes=target_lons, latitudes=target_lats)

        del target_lons, target_lats  # clean up

        outparameter = cls.resampleModel2Model(source,
                                               target,
                                               fields,
                                               radius=radius,
                                               new=True,
                                               add_reference=add_reference,
                                               add_geolocation=add_geolocation,
                                               prefix=prefix
                                               )
        outparameter.area_def = target_area_def
        outparameter._geolocation_fields[
            'time'] = source.get_geolocation_field('time')

        return outparameter
