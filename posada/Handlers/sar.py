# coding: utf-8

"""
.. module::posada.Handlers.sar

This module contains the abstract class for SAR Products/Parameters.

:copyright: Copyright 2011-2016 SOLab.


:copyright: Copyright 2013 Ifremer / Cersat / SOLab.
:license: Released under GPL v3 license, see :ref:'license'.

.. sectionauthor:: Alexander Myasoedov <mag@rshu.ru>
.. codeauthor:: Alexander Myasoedov <mag@rshu.ru>
.. author:: SOLab Dev Team
"""

import numpy as np

from ..Models.grid import Grid
from ..Models.swath import Swath

from ..Common.wind import windSpeed, windDirection
from ..Common.logger import get_logger
logger = get_logger()

try:
    from cmod.cmod_vect import cmod5n_forward as cmod5
except ImportError:
    logger.warning(
        "This POSADA instance doesn't include an installed CMOD. \n"
        "Setting up using pip inplace...")
    import pip
    pip.main(['install', '--user', '--upgrade',
              'https://github.com/SOLab/cmod/archive/master.zip'])

    from cmod.cmod_vect import cmod5n_forward as cmod5


__all__ = ['SAR']

DEFAULTSTYLE = {'lon': 'point',
                'lat': 'point',
                'incidence': 'point',
                'elevation': 'point',
                'digital_number': 'mean',
                'gamma0_lut': 'point',
                'beta0_lut': 'point',
                'sigma0_lut': 'point',
                'gamma0': 'mean',
                'beta0': 'mean',
                'sigma0': 'mean',
                'complex': 'mean',
                'azimuth_time': 'point',
                'slant_range_time': 'point',
                'doppler_centroid': 'point',
                'roughness': 'mean'}

# class SAR(Grid, Swath):


class SAR(object):
    """
    SAR Products generic storage class
    """

    def __init__(self, model):
        if isinstance(model, Grid):
            self.model = model
        elif isinstance(model, Swath):
            self.model = model
        return

    def get_infos(self):
        """Get all SAR informations as a ``dict``."""
        return self.model._datastorage.read_global_attributes()

    def print_infos(self):
        """Print formatted SAR informations."""
        for key, value in self.get_infos().items():
            if isinstance(value, dict):
                value2print = 'not printed dict'
            else:
                value2print = value
            print "%30s : %30s %s" % (key, value2print, type(value))

    def get_info(self, name):
        """Get one SAR information given a key name."""
        return self.get_infos()[name]

    def get_datanames(self):
        """
        """
        datanames = self.model._datastorage.get_fieldnames()
        if 'time' in datanames:
            datanames.remove('time')
        return datanames

    def get_data(self, name, extent=None, spacing=[1, 1], midazimuth=False,
                 midrange=False, style=None, blocksize=None):
        """Get numerical values of a SAR swath variable.

        Parameters
        ----------
        name : str
            Variable name
        extent : array_like, optional
            [azimuth_start, range_start, azimuth_stop, range_stop] pixel extent.
            Default to [0, 0, azimuth_size-1, range_size-1]. Note that the
            limits defined by extent are all included in output unlike python's
            slice notation.
        spacing : array_like, optional
            [azimuth, range] pixel spacing. Default to [1, 1]. If only one value
            is provided, spacing is assumed equal for azimuth and range.
        midazimuth : bool, optionnal
            Default to False.
        midrange : bool, optionnal
            Default to False.
        style : str, optionnal
            Setting style='point' means data are read/computed every spacing
            pixel. Setting style='mean' means data are read/computed at full
            resolution then averaged every spacing pixel.
            Note that when style='point', the values are approximately at the
            center of the bins used for averaging when style='mean'.
            By default, style is automatically guessed from the given variable
            name (see global DEFAULTSTYLE).
        blocksize : int, optional. Get data by blocks for memory issues

        Returns
        -------
        values : ndarray
            [azimuth, range] values.
        """
        # Manage inputs and make slices for mapper
        if extent is None:
            extent = self.extent_max()
        (ext, spa) = self._format_extent_spacing(extent, spacing,
                                                 midazimuth=midazimuth,
                                                 midrange=midrange)
        if style is None:
            style = DEFAULTSTYLE[name]
        if style == 'point':
            start = ext[0:2] + spa // 2
            stop = ext[2:4] - (spa - 1) + spa // 2 + 1
            step = spa
        elif style == 'mean':
            start = ext[0:2]
            stop = ext[2:4] + 1
            step = (1, 1)
        else:
            raise Exception('Unknown style :' + str(style))
        slices = [slice(int(start[0]), int(stop[0]), int(step[0])),
                  slice(int(start[1]), int(stop[1]), int(step[1]))]
        # Get data by blocks for memory issues
        # for style='point', blocksize is handled by mapper
        if blocksize is not None and style == 'mean' and spa.max() > 1:
            azsize = slices[0].stop - slices[0].start
            rasize = slices[1].stop - slices[1].start
            if azsize * rasize > blocksize:
                blrasize = rasize
                blazsize = np.maximum(blocksize // blrasize, spa[0])
                blazsize -= blazsize % spa[0]
                nblocks = np.ceil(azsize / float(blazsize))
                blext = np.array(ext, copy=True)
                for ibl in np.arange(nblocks):
                    blext[0] = ext[0] + ibl * blazsize
                    blext[2] = np.minimum(blext[0] + blazsize - 1, ext[2])
                    blvalues = self.get_data(name, extent=blext, spacing=spa,
                                             midazimuth=False, midrange=False,
                                             style=style, blocksize=None)
                    if ibl == 0:
                        dtype = blvalues.dtype
                        values = np.empty((azsize, rasize) / spa, dtype=dtype)
                    az0 = int(ibl * blazsize / spa[0])
                    values[az0:az0 + blvalues.shape[0], :] = blvalues
                return values
        # Read values
        if name == 'roughness':
            values = self.model._datastorage.read_values('sigma0',
                                                         slices=slices,
                                                         blocksize=blocksize)
        else:
            values = self.model._datastorage.read_values(name,
                                                         slices=slices,
                                                         blocksize=blocksize)
        # Rebin
        if style == 'mean' and spa.max() > 1:
            sha = (values.shape[0] / spa[0], spa[0],
                   values.shape[1] / spa[1], spa[1])
            values = values.reshape(sha).mean(-1).mean(1)
        # Roughness special case
        if name == 'roughness':
            inc = self.get_data('incidence', spacing=spacing, extent=extent,
                                midazimuth=midazimuth, midrange=midrange,
                                blocksize=blocksize)
            pol = self.get_info('polarisation')
            values = compute_roughness(values, inc, pol)
        return values

    def _format_extent_spacing(self, extent, spacing, midazimuth=False,
                               midrange=False):
        """Format (and check) extent and spacing."""
        # Check extent
        ext = np.round(extent).flatten()
        if ext.size != 4:
            raise Exception('extent must contain 4 elements')
        extmax = self.extent_max()
        if (ext[0:2] < extmax[0:2]).any() or (ext[2:4] > extmax[2:4]).any():
            exttmp = np.array(ext)
            ext[0:2] = np.maximum(ext[0:2], extmax[0:2])
            ext[2:4] = np.minimum(ext[2:4], extmax[2:4])
            print 'Warning : extent is outside SAR image, ' + str(exttmp) +\
                ' becomes ' + str(ext)
        if (ext[0:2] > ext[2:4]).any():
            raise Exception('extent[0:2] must be less or equal than ' +
                            'extent[2:4]')
        # Check spacing
        spa = np.round(spacing).flatten()
        if spa.size == 1:
            spa = np.repeat(spa[0], 2)
        elif spa.size == 2:
            pass
        else:
            raise Exception('spacing must contain 1 or 2 elements')
        if (spa < [1, 1]).any():
            spatmp = np.array(spa)
            spa = np.maximum(spa, [1, 1])
            print 'Warning : spacing too small, ' + str(spatmp) + ' becomes ' +\
                str(spa)
        if (spa > ext[2:4] - ext[0:2] + 1).any():
            spatmp = np.array(spa)
            spa = np.minimum(spa, ext[2:4] - ext[0:2] + 1)
            print 'Warning : spacing too large, ' + str(spatmp) + ' becomes ' +\
                str(spa)
        # Make extent to be spacing modulo
        ext[2:4] = ext[2:4] - (ext[2:4] - ext[0:2] + 1) % spa
        # 1D extent
        if midazimuth == True:
            dim = (ext[2] - ext[0] + 1) / spa[0]
            ext[0:3:2] = ext[0] + (dim - 1) // 2 * spa[0] + [0, spa[0] - 1]
        if midrange == True:
            dim = (ext[3] - ext[1] + 1) / spa[1]
            ext[1:4:2] = ext[1] + (dim - 1) // 2 * spa[1] + [0, spa[1] - 1]
        return (ext, spa)

    def extent_max(self):
        """Get extent for the whole SAR image."""
        return np.array((0, 0, self.get_info('number_of_lines') - 1,
                         self.get_info('number_of_samples') - 1))

    def extent_burst(self, burst, valid=True):
        """Get extent for a SAR image burst."""
        nbursts = self.get_info('number_of_bursts')
        if nbursts == 0:
            raise Exception('No bursts in SAR image')
        if burst < 0 or burst >= nbursts:
            raise Exception('Invalid burst index number')
        if valid == True:
            burst_list = self.get_info('burst_list')
            extent = np.copy(burst_list['valid_location'][burst, :])
        else:
            extent = self.extent_max()
            nlines = self.get_info('lines_per_burst')
            extent[0:3:2] = [nlines * burst, nlines * (burst + 1) - 1]
        return extent

    def extent_from_lonlat(self, lonlat, dist):
        """Get extent from a geographic point and a distance."""
        raise Exception('Not yet implemented !')

    def ground_spacing(self, range_index=None, extent=None):
        """Get SAR image ground spacing.

        Parameters
        ----------
        range_index : int, optionnal
            Range index (pixel) at which the range ground spacing is computed
            for SLC products case. By default, corresponds to middle swath.
        extent : array_like, optionnal
            If provided, range_index is computed from the middle range extent.

        Returns
        -------
        ground_spacing : ndarray
            [azimuth, range] ground spacing in meters.

        Notes
        -----
        For GRD products, range_index and extent are ignored.
        """
        ground_spacing = np.array((self.get_info('azimuth_pixel_spacing'),
                                   self.get_info('range_pixel_spacing')))
        if self.get_info('product') == 'SLC':
            ext = self.extent_max()
            if extent is not None:
                ext[1:4:2] = extent[1:4:2]
            else:
                if range_index is not None:
                    ext[1:4:2] = range_index
            inc = self.get_data('incidence', extent=ext, midazimuth=True,
                                midrange=True)
            ground_spacing[1] /= np.sin(inc * np.pi / 180)
        return ground_spacing

    def pixels2meters(self, pixels, **kwargs):
        """Convert pixels to meters."""
        meters = pixels * self.ground_spacing(**kwargs)
        return meters

    def meters2pixels(self, meters, **kwargs):
        """Convert meters to pixels."""
        pixels = meters / self.ground_spacing(**kwargs)
        return pixels


def cdop_func(x):
    """
    """
    return 1. / (1. + np.exp(-x))


def cdop(u10, phi, inc, pol):
    """
    """
    # Check inputs
    sizes = np.array([np.size(inc), np.size(u10), np.size(phi)])
    size = sizes.max()
    if ((sizes != size) & (sizes != 1)).any():
        raise Exception('Inputs sizes do not agree.')
    if pol not in ['VV', 'HH']:
        raise Exception('Unknown polarisation : ' + pol)
    # NN coefficients (W=weights and B=biases)
    # (coefficient names in mouche2012 are given)
    if pol == 'VV':
        # lambda[0:2,1]
        B1 = np.array([-0.343935744939, 0.108823529412, 0.15],
                      dtype='float32')
        # lambda[0:2,0]
        W1 = np.array([0.028213254683, 0.0411764705882, .00388888888889],
                      dtype='float32')
        # omega[i,0]
        B2 = np.array([14.5077150927, -11.4312028555, 1.28692747109,
                       -1.19498666071, 1.778908726, 11.8880215573,
                       1.70176062351, 24.7941267067, -8.18756617111,
                       1.32555779345, -9.06560116738],
                      dtype='float32')
        # omega[i,[3,2,1]]
        W2 = np.array([[19.7873046673, 22.2237414308, 1.27887019276],
                       [2.910815875, -3.63395681095, 16.4242081101],
                       [1.03269004609, 0.403986575614, 0.325018607578],
                       [3.17100261168, 4.47461213024, 0.969975702316],
                       [-3.80611082432, -6.91334859293, -0.0162650756459],
                       [4.09854466913, -1.64290475596, -13.4031862615],
                       [0.484338480824, -1.30503436654, -6.04613303002],
                       [-11.1000239122, 15.993470129, 23.2186869807],
                       [-0.577883159569, 0.801977535733, 6.13874672206],
                       [0.61008842868, -0.5009830671, -4.42736737765],
                       [-1.94654022702, 1.31351068862, 8.94943709074]],
                      dtype='float32')
        # gamma[0]
        B3 = np.array(4.07777876994, dtype='float32')
        # gamma[1:11]
        W3 = np.array([7.34881153553, 0.487879873912, -22.167664703,
                       7.01176085914, 3.57021820094, -7.05653415486,
                       -8.82147148713, 5.35079872715, 93.627037987,
                       13.9420969201, -34.4032326496],
                      dtype='float32')
        # beta
        B4 = np.array(-52.2644487109, dtype='float32')
        # alpha
        W4 = np.array(111.528184073, dtype='float32')
    elif pol == 'HH':
        # lambda[0:2,1]
        B1 = np.array([-0.342097701547, 0.118181818182, 0.15],
                      dtype='float32')
        # lambda[0:2,0]
        W1 = np.array([0.0281843837385, 0.0318181818182, 0.00388888888889],
                      dtype='float32')
        # omega[i,0]
        B2 = np.array([1.30653883096, -2.77086154074, 10.6792861882,
                       -4.0429666906, -0.172201666743, 20.4895916824,
                       28.2856865516, -3.60143441597, -3.53935574111,
                       -2.11695768022, -2.57805898849],
                      dtype='float32')
        # omega[i,[3,2,1]]
        W2 = np.array([[-2.61087309812, -0.973599180956, -9.07176856257],
                       [-0.246776181361, 0.586523978839, -0.594867645776],
                       [17.9261562541, 12.9439063319, 16.9815377306],
                       [0.595882115891, 6.20098098757, -9.20238868219],
                       [-0.993509213443, 0.301856868548, -4.12397246171],
                       [15.0224985357, 17.643307099, 8.57886720397],
                       [13.1833641617, 20.6983195925, -15.1439734434],
                       [0.656338134446, 5.79854593024, -9.9811757434],
                       [0.122736690257, -5.67640781126, 11.9861607453],
                       [0.691577162612, 5.95289490539, -16.0530462],
                       [1.2664066483, 0.151056851685, 7.93435940581]],
                      dtype='float32')
        # gamma[0]
        B3 = np.array(2.68352095337, dtype='float32')
        # gamma[1:11]
        W3 = np.array([-8.21498722494, -94.9645431048, -17.7727420108,
                       -63.3536337981, 39.2450482271, -6.15275352542,
                       16.5337543167, 90.1967379935, -1.11346786284,
                       -17.57689699, 8.20219395141],
                      dtype='float32')
        # beta
        B4 = np.array(-66.9554922921, dtype='float32')
        # alpha
        W4 = np.array(136.216953823, dtype='float32')
    # Make inputs as a single matrix (and clip phi in [0,180])
    inputs = np.zeros((3, size), dtype='float32')
    for ivar, var in enumerate((inc, u10, phi)):
        if sizes[ivar] == 1:
            inputs[ivar, :] = np.repeat(var, size)
        else:
            inputs[ivar, :] = np.ravel(var)
        if ivar == 2:
            inputs[ivar, :] = np.abs(
                ((inputs[ivar, :] + 180) % 360) - 180)
        inputs[ivar, :] *= W1[ivar]
        inputs[ivar, :] += B1[ivar]
    # Compute CDOP
    B2 = np.tile(B2.reshape((11, 1)), (1, size))
    dop = W4 * \
        cdop_func(np.dot(W3, cdop_func(
            np.dot(W2, inputs) + B2)) + B3) + B4
    # Reshape output
    # (using the shape of input which have the maximum ndim)
    ndims = np.array([np.ndim(inc), np.ndim(u10), np.ndim(phi)])
    tmp = np.where(sizes == size)[0]
    ivar = tmp[ndims[tmp].argmax()]
    shp = np.shape((inc, u10, phi)[ivar])
    dop = dop.reshape(shp)
    return dop


def compute_roughness(sigma0, incidence, polarisation):
    """Compute sea surface roughness.

    Parameters
    ----------
    sigma0 : ndarray
        NRCS backscatter.
    incidence : ndarray
        Incidence angle in degrees.
    polarisation : str
        'VV' or 'HH' or 'VH' or 'HV'

    Returns
    -------
    ndarray
    """
    if polarisation == 'VV':  # Use cmod5
        sigma0_vv = cmod5(10, 45, incidence)
        return sigma0 / sigma0_vv
    elif polarisation == 'HH':  # Use cmod5 and Thompson polarisation ratio
        sigma0_vv = cmod5(10, 45, incidence)
        alpha = 0.7
        polrat = (1 + 2 * np.tan(incidence * np.pi / 180)**2)**2 / \
                 (1 + alpha * np.tan(incidence * np.pi / 180)**2)**2
        return sigma0 / sigma0_vv * polrat
    elif polarisation == 'VH' or polarisation == 'HV':  # Use simple model
        # nrcs_vh_db = 0.580*wsp - 35.652
        # nrcs_vh_lin = 10^(nrcs_vh_db/10.)
        sigma0_cross = 10**((0.58 * 10 - 35.652) / 10)
        return sigma0 / sigma0_cross
    else:
        raise Exception('Unknown polarisation : ' + polarisation)
