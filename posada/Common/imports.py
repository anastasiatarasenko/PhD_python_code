# coding: utf-8

"""
.. module::posada.Common.imports

This module contains the main functions imports used in POSADA.

:copyright: Copyright 2011-2016 SOLab.

.. sectionauthor:: Alexander Myasoedov <mag@rshu.ru>
.. codeauthor:: Alexander Myasoedov <mag@rshu.ru>
.. author:: SOLab Dev Team
"""

from collections import OrderedDict

from numpy import floor, size, linspace, arange, row_stack, delete, round, \
    double, asarray, array, fliplr, flipud, mean, pi, exp, \
    sin, cos, radians, log10, ma, maximum, minimum, repeat, zeros, reshape, \
    float64, dot, sqrt, diff, isnan, roll, sign, ones, \
    uint8, uint16, uint32, uint64

from numpy.ma import masked_where, where, MaskedArray, masked_equal, concatenate

from math import atan, pow

from multiprocessing import cpu_count

from scipy.interpolate import RectBivariateSpline, RectSphereBivariateSpline

import pyresample as pr
from pyproj import Proj
from osgeo.osr import SpatialReference
from shapely.geometry import Polygon
