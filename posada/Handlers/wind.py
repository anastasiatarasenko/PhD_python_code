# coding: utf-8

"""
.. module::posada.Handlers.wind

This module contains the abstract class for WIND Products/Parameters.

:copyright: Copyright 2011-2016 SOLab.

.. sectionauthor:: Alexander Myasoedov <mag@rshu.ru>
.. codeauthor:: Alexander Myasoedov <mag@rshu.ru>
.. author:: SOLab Dev Team
"""

from collections import OrderedDict

from copy import copy

from ..Models.grid import Grid
from ..Models.swath import Swath

from ..Models.Field import Field
from ..Models.Variable import Variable

from ..Common.wind import windSpeed, windDirection
from ..Common.logger import get_logger
logger = get_logger()

__all__ = ['Wind']


class Wind(object):
    """
    WIND Products generic storage class
    """

    def __init__(self, model):
        if isinstance(model, Grid):
            self.model = model
        elif isinstance(model, Swath):
            self.model = model
        else:
            raise NotImplementedError
        return

    def add_wind_speed_field(self):
        """
        """

        # renaming wind speed components to fit 'u', 'v'
        stringsU = ('10u', 'u10')
        stringsV = ('10v', 'v10')
        for field in self.model._fields:
            if field not in self.model._geolocation_fields:
                for s in stringsU:
                    if s in field:
                        self.model._fields[field].set_name('u')
                        self.model._fields = OrderedDict([
                            ('u' if k == field else k, v)
                            for k, v in self.model._fields.items()])
                for s in stringsV:
                    if s in field:
                        self.model._fields[field].set_name('v')
                        self.model._fields = OrderedDict([
                            ('v' if k == field else k, v)
                            for k, v in self.model._fields.items()])

        u = self.model.get_field('u').get_values()
        v = self.model.get_field('v').get_values()

        # create a field for wind_speed and wind_direction
        globalattrs = self.model.get_metadata()
        auth = None
        if 'Conventions' in globalattrs:
            auth = self.read_global_attribute('Conventions')
        else:
            auth = 'solab.rshu.ru'

        wind_speed = windSpeed(u, v)

        uField = self.model.get_field('u')

        # Adding Wind Speed
        variable = Variable(
            shortname='wind_speed',
            description='10 metre wind speed',
            authority=auth,
            standardname='wind speed'
        )
        field = Field(
            variable,
            dimensions=copy(uField.dimensions),
            datatype=uField.datatype,
            fillvalue=uField.fillvalue,
            values=wind_speed,
            units=uField.units,
            valid_min=0, valid_max=35
        )

        self.model.add_field(field)

    def add_wind_dir_field(self):
        """
        """

        # renaming wind speed components to fit 'u', 'v'
        stringsU = ('10u', 'u10')
        stringsV = ('10v', 'v10')
        for field in self.model._fields:
            if field not in self.model._geolocation_fields:
                for s in stringsU:
                    if s in field:
                        self.model._fields[field].set_name('u')
                        self.model._fields = OrderedDict([
                            ('u' if k == field else k, v)
                            for k, v in self.model._fields.items()])
                for s in stringsV:
                    if s in field:
                        self.model._fields[field].set_name('v')
                        self.model._fields = OrderedDict([
                            ('v' if k == field else k, v)
                            for k, v in self.model._fields.items()])

        u = self.model.get_field('u').get_values()
        v = self.model.get_field('v').get_values()

        # create a field for wind_speed and wind_direction
        globalattrs = self.model.get_metadata()
        auth = None
        if 'Conventions' in globalattrs:
            auth = self.read_global_attribute('Conventions')
        else:
            auth = 'solab.rshu.ru'

        wind_dir = windDirection(u, v)

        uField = self.model.get_field('u')

        # Adding Wind Direction
        variable = Variable(
            shortname='wind_from_direction',
            description='Direction from which wind is blowing',
            authority=auth,
            standardname='wind direction'
        )
        field = Field(
            variable,
            dimensions=copy(uField.dimensions),
            datatype=uField.datatype,
            fillvalue=uField.fillvalue,
            values=wind_dir,
            units='degree',
            valid_min=0, valid_max=360
        )

        self.model.add_field(field)

    def get_wind_speed_values(self):
        """
        """

        if 'wind_speed' in self.model._fields:
            logger.warning("Field wind_speed already exists in Parameter. " +
                           "Will use existing!")
        else:
            self.add_wind_speed_field()

        return self.model.get_field('wind_speed').get_values()

    def get_wind_dir_values(self):
        """
        """

        if 'wind_from_direction' in self.model._fields:
            logger.warning("Field wind_dir already exists in Parameter. " +
                           "Will use existing!")
        else:
            self.add_wind_dir_field()

        return self.model.get_field('wind_from_direction').get_values()
