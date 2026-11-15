__import__('pkg_resources').declare_namespace(__name__)

# encoding: utf-8
__all__ = ['BaseReader',
           'netcdf', 'grib', 'safe']

from pkg_resources import WorkingSet, DistributionNotFound
from ..Common.logger import get_logger
logger = get_logger()

# Import what readers are avaialble.
packages = WorkingSet()

try:
    packages.require('netCDF4')
except DistributionNotFound:
    logger.warning('Python netCDF4 package is required for netCDF. '
                   'No netCDF compatible modules are available on this posada instance')
except ImportError:
    logger.exception('A netCDF reader failed to load, and is unavailable:')

try:
    packages.require('pygrib')
    # from . import gribfile
except DistributionNotFound:
    logger.warning('Python pygrib package is required for GRIB. '
                   'No GRIB compatible modules are available on this posada instance')
except ImportError:
    logger.exception('A GRIB reader failed to load, and is unavailable:')
