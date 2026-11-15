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

from ..Models.grid import Grid
from ..Models.swath import Swath

from .. import READ_ONLY, WRITE_NEW, READ_WRITE, SAVED, NOTSAVED

from ..Common.logger import get_logger
logger = get_logger()

__all__ = ['BaseExtractor']


def BaseResampler(projIn, projOut):
    """
    Generic class for resampling/reprojecting data to different projections
    """

    if isinstance(projIn, Grid):
        self.projIn = projIn
    elif isinstance(projIn, Swath):
        self.projIn = projIn
    else:
        raise NotImplementedError

    if isinstance(projOut, Grid):
        self.projOut = projOut
    elif isinstance(projOut, Swath):
        self.projOut = projOut
    else:
        raise NotImplementedError
    return


def resample(self):
    """
    """
    return NotImplementedError
