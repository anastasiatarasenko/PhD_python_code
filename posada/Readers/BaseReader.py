# coding: utf-8

"""
.. module::posada.Readers.BaseReader

This module contains the abstract class for any reader. A reader is used to
read the content of a storage (file). It is designed to allow the mapping of
this content to a Model (also refered as a 'Parameter'), i.e. a class from
:mod:'posada.Models' package.

It can also be used independently, providing the user with a single API
to access any data format.

:copyright: Copyright 2013 Ifremer / Cersat / SOLab.
:license: Released under GPL v3 license, see :ref:'license'.

.. sectionauthor:: Jeff Piolle <jfpiolle@ifremer.fr>, Alexander Myasoedov <mag@rshu.ru>
.. codeauthor:: Jeff Piolle <jfpiolle@ifremer.fr>, Alexander Myasoedov <mag@rshu.ru>
.. author:: SOLab Dev Team
"""

import datetime
import os
from numpy import diff, mean, abs

from abc import ABCMeta, abstractmethod
from collections import OrderedDict

from ..Models import Field
from .. import READ_ONLY, WRITE_NEW, READ_WRITE, SAVED, NOTSAVED

from ..Common.logger import get_logger
logger = get_logger()

__author__ = 'Alexander Myasoedov'
__email__ = 'mag@rshu.ru'
__created__ = datetime.datetime(2016, 5, 11)


class CorruptFileException(Exception):
    """Exception class for corrupted or unreadable files."""
    pass


class FieldHandler(object):
    """
    Storage pointer to the data of a particular field of one single feature
    (files can store several features).
    """

    def __init__(self, reader, fieldname, index=None, status=NOTSAVED):
        self.reader = reader
        self.fieldname = fieldname
        self.index = index
        self.status = status
        return

    def set_unsaved(self):
        self.status = NOTSAVED

    def set_saved(self):
        self.status = SAVED

    def is_saved(self):
        return (self.status == SAVED)


class BaseReader(object):
    """
    An abstract class for any reader.

    This does not actually open the corresponding file.
    Explicit opening must be done by calling the :func:'open()' method.

    Args:
        url (str): full path to the file.

        mode (enum): access mode (READ_ONLY, WRITE_NEW, READ_WRITE)
    """
    __metaclass__ = ABCMeta

    #@profile
    def __init__(self,
                 url=None,
                 urlseries=None,
                 mode=READ_ONLY,
                 **kwargs):
        """
        """
        if not url and not urlseries:
            raise Exception('url or urlseries must be provided')
        object.__init__(self)
        # memorize opening arguments
        self.args = kwargs
        self._feature_type = None
        self._feature_count = None
        self._handler = None
        self.view = None
        self._mode = mode
        self._url = url
        self._urlseries = urlseries
        self.datamodel_geolocation_dims = None
        if url and urlseries:
            raise Exception('BaseReader: url and fileSeries are exclusive')
        elif url and 'http://' not in url:
            existant = self.exists(url)
            if mode == READ_ONLY and not existant:
                logger.error(
                    "File {} not existing and can not be read".format(url))
                raise IOError(
                    "File {} not existing and can not be read".format(url))
            elif mode == WRITE_NEW and existant:
                logger.error("File '%s' already existing and can not be\
                             created", url)
                raise IOError("File '%s' already existing and can not be\
                    created", url)
            elif mode == READ_WRITE and not existant:
                logger.warning(
                    "File {} is not existing. Will be created".format(url))
        elif urlseries:
            self._urlseries.storageClass = self.__class__
        return

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()

    def sync(self):
        """force physical writing of the content on disk."""
        raise NotImplementedError

    @classmethod
    def exists(cls, url):
        """tests if 'url' is an existing resource"""
        if os.path.exists(url):
            return True
        else:
            return False

    def check_geolocation_dimensions(self, dimensions):
        """Returns True if the requested dimensions exist.

        Used for internal purpose and should not be called directly.
        """
        return all([self.get_matching_dimname(d) for d in dimensions])

    def check_geolocation_fields(self, fields):
        """Returns True if the requested geolocation fields exist.

        Used for internal purpose and should not be called directly.
        """
        return all([self.has_geolocation_field(f) for f in fields])

    def has_geolocation_field(self, fieldname):
        """Returns True if the field 'fieldname' exists."""
        return (self.get_geolocation_field(fieldname) is not None)

    @abstractmethod
    def get_geolocation_field(self, fieldname):
        """Returns the equivalent field name in the file format for a standard
        geolocation field (lat, lon, time, z).

        Used for internal purpose and should not be called directly.

        Args:
            fieldname (str): name of the standard geolocation field (lat, lon
                or time)

        Return:
            str: name of the corresponding field in the native file format.
                Returns None if no matching is found
        """
        raise NotImplementedError

    @abstractmethod
    def get_matching_dimname(self, dimname):
        """Returns the equivalent name in the native format for a standard
        dimension.

        This is a translation of the standard names to native ones. It is used
        for internal purpose only and should not be called directly.

        The standard dimension names are:

        * x, y, time for :class:'~posada.datamodel.grid.Grid'
        * row, cell, time for :class:'~posada.datamodel.swath.Swath' or
          :class:'~posada.datamodel.image.Image'

        To be derived when creating an inherited data reader class. This is
        mandatory for geolocation dimensions which must be standard.

        Args:
            dimname (str): standard dimension name.

        Returns:
            str: return the native name for the dimension. Return 'dimname' if
                the input dimension has no standard name.

        See Also:
            see :func:'get_standard_dimname' for the reverse operation
        """
        raise NotImplementedError

    @abstractmethod
    def get_standard_dimname(self, dimname):
        """
        Returns the equivalent standard dimension name for a
        dimension in the native format.

        This is a translation of the native names to standard ones. It is used
        for internal purpose and should not be called directly.

        To be derived when creating an inherited data reader class. This is
        mandatory for geolocation dimensions which must be standard.

        Args:
            dimname (string): native dimension name

        Return:
            str: the (translated) standard name for the dimension. Return
            'dimname' if the input dimension has no standard name.

        See Also:
            see :func:'get_matching_dimname' for the reverse operation
        """
        raise NotImplementedError

    def get_url(self):
        return self._url

    def get_basename(self):
        return os.path.basename(self._url)

    #@profile
    def get_handler(self):
        """Returns the handler to the physical file or resource from which data
        are to be read.

        In case the reader was initialised with a URLSeries, the handler to the
        reference time URL of the series is returned.
        """
        if self._handler is None:
            if self._url is None and self._urlseries is not None:
                ref_url = self._urlseries.get_example_url()
                self._handler = self.open(ref_url)
            elif self._url is not None and self._urlseries is None:
                if self._mode in [READ_WRITE, READ_ONLY]\
                        and not self.exists(self._url):
                    raise Exception("File does not exist: ", self._url)
                self._handler = self.open()
        return self._handler

    def is_opened(self):
        """
        """
        if self._handler is None:
            return False
        if self._handler == False:
            return False
        return True

    def is_writable(self):
        """Returns True if the storage is opened in write mode"""
        return self._mode == WRITE_NEW

    def is_readonly(self):
        """Returns True if the storage is opened in write mode"""
        return self._mode == READ_ONLY

    def get_size(self):
        """return the product file size, in octets"""
        return os.path.getsize(self._url)

    def get_creation_date(self):
        """return the date the product was generated (NOT the file date!)"""
        return datetime.datetime.fromtimestamp(os.path.getctime(self._url))

    @abstractmethod
    def open(self, view=None, datamodel=None, datamodel_geolocation_dims=None):
        """Open the file (or any other type of storage)

        Args:
            view (dict, optional): a dictionary where keys are dimension names
                and values are slices. A view can be set on a file, meaning
                that only the subset defined by this view will be accessible.
                This view is expressed as any subset (see :func:'get_values').
                For example:

                view = {'time':slice(0,0), 'lat':slice(200,300),
                'lon':slice(200,300)}

            datamodel (str): type of feature read or written. Internal argument
                only used by the classes from :mod:'~posada.Models'
                package. Can be 'Grid', 'Swath', 'Pass' etc...

            datamodel_geolocation_dims (list, optional): list of the name of
                the geolocation dimensions defining the data model to be read 
                in the file. Optional argument, only used by the datamodel
                classes, in case the reader class can store different types of
                data models.

        Returns:
            a handler on the opened file
        """
        self.view = view
        if datamodel is not None:
            self._feature_type = datamodel

    @abstractmethod
    def close(self):
        """Close handler on storage"""
        raise NotImplementedError

    @abstractmethod
    def read_values(self, fieldname, slices=None, **kwargs):
        """Read the data of a field.

        Args:
            fieldname (str): name of the field which to read the data from

            slices (list of slice, optional): list of slices for the field if
                subsetting is requested. A slice must then be provided for each
                field dimension. The slices are relative to the opened view (see
                :func:open) if a view was set when opening the file.

        Return:
            MaskedArray: array of data read. Array type is the same as the
                storage type.
        """
        raise NotImplementedError

    @abstractmethod
    def read_field(self, fieldname):
        """
        Returns the :class:'posada.Field.Field' object corresponding to
        the requested fieldname.

        The :class:'posada.Field.Field' class contains all the metadata
        describing a field (equivalent to a variable in netCDF).

        Args:
            fieldname (str): name of the field

        Returns:
            :class:'posada.Field.Field': the corresponding field object
        """
        raise NotImplementedError

    @abstractmethod
    def write_field(self, fieldname):
        """Writes the field data on disk.

        Args:
            fieldname (str): name of the field to write.
        """
        raise NotImplementedError

    @abstractmethod
    def read_fillvalue(self, fieldname):
        """Read the fill value of a field.

        Args:
            fieldname (str): name of the field.

        Returns:
            number or char or str: fill value of the field. The type is the
                as the type of the data in the field.
        """
        raise NotImplementedError

    @abstractmethod
    def create_field(self, field, dim_translation=None):
        """Creates a new field in the reader.

        Creates the field structure but don't write yet its values array.

        Args:
            field (Field): the field to be created.

        See also:
            :func:'write_field' for writing the values array.
        """
        raise NotImplementedError

    @abstractmethod
    def create_dim(self, dimname, size=None):
        """Add a new dimension.

        Args:
            dimname (str): name of the dimension.
            size (int): size of the dimension (unlimited if None)
        """
        raise NotImplementedError

    @abstractmethod
    def get_fieldnames(self):
        """Returns the list of geophysical fields stored for the feature.

        The geolocation field names are excluded from this list.

        Returns:
            list<string>: list of field names
        """
        raise NotImplementedError

    @abstractmethod
    def get_dimensions(self, fieldname=None):
        """Returns the dimension's standard names of a file or a field in the
        file.

        Args:
            fieldname (str): the name of the field from which to get the
                dimensions. For a geolocation field, use the posada.standard
                name (time, lat, lon), though native field name will work too.

        Returns:
            tuple<str>: the standard dimensions of the field or file.
        """
        raise NotImplementedError

    def get_dimsize(self, dimname):
        """Returns the size of a dimension.

        Args:
            dimname (str): name of the dimension.

        Returns:
            int: size of the dimension.
        """
        raise NotImplementedError

    def get_full_dimensions(self, fieldname):
        """Returns the dimension names and sizes of a field.

        Args:
            fieldname (str): name of the field

        Returns:
            OrderedDict: ordered dictionary where keys are the dimension names
                and values are their respective size
        """
        dims = OrderedDict()
        for d in self.get_dimensions(fieldname):
            dims[d] = self.get_dimsize(
                self.get_matching_dimname(d)
            )
        return dims

    @abstractmethod
    def read_field_attributes(self, fieldname):
        """Returns the specific attributes of a field.

        Args:
            fieldname (str): name of the field.

        Returns:
            dict<string, string or number or datetime>: a dictionary where keys
                are the attribute names.
        """
        raise NotImplementedError

    @abstractmethod
    def read_global_attributes(self):
        """Returns the names of the global attributes.

        Returns:
            list<str>: the list of the attribute names.
        """
        raise NotImplementedError

    @abstractmethod
    def write_global_attributes(self, attrs):
        """Write the global attributes of the file.

        Args:
            attrs (dict<string, string or number or datetime>): a dictionary
                containing the attributes names and values to be written.
        """
        raise NotImplementedError

    @abstractmethod
    def read_global_attribute(self, name):
        """Returns the value of a global attribute.

        Args:
            name (str): name of the global attribute.

        Returns:
            str, number or datetime: value of the corresponding attribute.
        """
        raise NotImplementedError

    @abstractmethod
    def get_start_time(self):
        """Returns the minimum date of the file temporal coverage.

        Returns:
            datetime: start time of the data in file.
        """
        raise NotImplementedError

    @abstractmethod
    def get_end_time(self):
        """Returns the maximum date of the file temporal coverage.

        Returns:
            datetime: end time of the data in file.
        """
        raise NotImplementedError

    def get_spatial_resolution_in_deg(self):
        """Returns the average spatial resolution in degrees.

        Returns:
            float: resolution, in degrees.
        """
        return abs(mean(diff(self.get_geolocation_field('lat'))))

    def get_product_version(self):
        """return the product version"""
        raise NotImplementedError

    @abstractmethod
    def get_bbox(self):
        """Returns the bounding box of the feature, as a tuple.

        Returns:
            tuple: bbox expressed as (lonmin, latmin, lonmax, latmax)
        """
        raise NotImplementedError

    def get_collection_id(self):
        """return the identifier of the product collection"""
        raise NotImplementedError

    def get_orbit_number(self):
        """In the case of a satellite orbit file, returns the orbit number.

        Returns:
            int: the orbit number
        """
        raise NotImplementedError

    def get_cycle_number(self):
        """In the case of a satellite orbit file, returns the cycle number.

        Returns:
            int: the cycle number
        """
        raise NotImplementedError

    def get_naming_authority(self):
        globalattrs = self.read_global_attributes()
        auth = None
        if 'Conventions' in globalattrs:
            auth = self.read_global_attribute('Conventions')
        else:
            auth = None
        del globalattrs
        return auth

    def get_field_handler(self, fieldname):
        if fieldname in self.get_fieldnames()\
                or self.has_geolocation_field(fieldname):
            status = SAVED
        else:
            status = NOTSAVED
        descr = FieldHandler(self, fieldname, status=status)
        return descr

    @classmethod
    def _fill_slices(cls, slices, dimsizes):
        """Fill (and check) slices."""
        # If you change the policy here, please make sure it is consistent
        # with slices management in read_values() and get_dimsize()
        if len(slices) != len(dimsizes):
            raise Exception('Unexpected slices list length : {}'
                            .format(slices))
        filled_slices = []
        for sli, dimsize in zip(slices, dimsizes):
            if sli.start is not None and abs(sli.start) > dimsize:
                raise Exception('Unexpected slice start : {}'.format(sli))
            if sli.stop is not None and abs(sli.stop) > dimsize:
                raise Exception('Unexpected slice stop : {}'.format(sli))
            start, stop, step = sli.indices(dimsize)
            if (step > 0 and start > stop) or (step < 0 and start < stop):
                raise Exception('Unexpected slice start Vs stop : {}'
                                .format(sli))
            # Special case : negative step until the first element
            if stop == -1:
                stop = None
            filled_slices.append(slice(start, stop, step))
        return filled_slices

    @classmethod
    def _combine_slices(cls, view, slices, dimnames, dimsizes):
        """Combines the slices from the view and within the view"""
        # view slicing
        if view is not None:
            # make sure viewslice info is complete
            viewslices = Field.Field.format_slices(view, dimnames)
            viewslices = BaseReader._fill_slices(viewslices, dimsizes)
        # standard slicing
        if slices is not None:
            if view is None:
                finalslices = BaseReader._fill_slices(slices, dimsizes)
            else:
                # slices are relative to the opened view
                viewdimsizes = dimsizes
                finalslices = BaseReader._fill_slices(slices, viewdimsizes)
                for index, (sli, vwsli) in enumerate(zip(finalslices,
                                                         viewslices)):
                    vw_first = vwsli.start
                    first = sli.start
                    if sli.stop is None:
                        last = 0
                    else:
                        last = sli.stop - sli.step / abs(sli.step)
                    first = vw_first + first * vwsli.step
                    last = vw_first + last * vwsli.step
                    start = first
                    step = vwsli.step * sli.step
                    stop = last + step / abs(step)
                    if stop == -1:
                        stop = None
                    finalslices[index] = slice(start, stop, step)
        elif view is not None:
            finalslices = viewslices
        else:
            finalslices = [slice(None) for i in range(len(dimnames))]
            finalslices = BaseReader._fill_slices(finalslices, dimsizes)
        return finalslices

    @classmethod
    def _adjust_dimsize(cls, view, dimname, dimsize):
        """correct the size of a dimension with respect to a view"""
        if view is not None and dimname in view:
            viewslice = BaseReader._fill_slices([view[dimname]],
                                                [dimsize])[0]
            start, stop, step = viewslice.start, viewslice.stop, viewslice.step
            if stop is None:
                stop = -1
            dimsize = 1 + (abs(stop - start) - 1) / abs(step)
        return dimsize
