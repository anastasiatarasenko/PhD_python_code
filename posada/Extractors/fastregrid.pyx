"""

    Installation
        written with Cython. To compile:

    .. code-block:: python

       python setup.py build_ext --inplace
"""
import numpy as np
# "cimport" is used to import special compile-time information
# about the numpy module (this is stored in a file numpy.pxd which is
# currently part of the Cython distribution).
cimport numpy as np

# "ctypedef" assigns a corresponding compile-time type to DTYPE_t. For
# every type in the numpy module there's a corresponding compile-time
# type with a _t-suffix.
ctypedef np.float32_t DTYPEf32_t
ctypedef np.float64_t DTYPEf64_t


def fastregrid(lats, lons, values,
               latmin=-90., latmax=90., lonmin=-180., lonmax=180.,
               lonbinsize=1., latbinsize=1.):
    """
    Place unevenly spaced 2D data on a simple grid by 2D binning (mean).

    Args
        lat (ndarray (1D) of numpy float32): the independent data latitude-axis
            of the grid.
        lon (ndarray (1D) of numpy float32): the independent data longitude
            axis of the grid.
        values (ndarray (1D)): the dependent data in the form
            values = f(lat,lon).
        latmin, latmax, lonmin, lonmax (float, optional): the limits of the
            grid on which to remap the input data
        lonbinsize, latbinsize (float, optional): the full width and height of
             each bin on the grid.  If each bin is a cube, then this is the x
             and y dimension.  This is the step in both directions, x and y.
             Defaults to 1.

    Returns
        grid (ndarray (2D)): the evenly gridded data.  The value of each cell
            is the sum value of the contents of the bin.
        bins (ndarray (2D)): A grid the same shape as 'grid', except the value
            of each cell is the number of points in that bin.
    """
    if values.dtype == 'float64':
        result, nbsampl = __fastregrid_float64(
            lats.astype('float32'),
            lons.astype('float32'),
            values,
            latmin, latmax, lonmin, lonmax,
            lonbinsize, latbinsize
            )
    elif values.dtype == 'float32':
        result, nbsampl = __fastregrid_float32(
            lats.astype('float32'),
            lons.astype('float32'),
            values,
            latmin, latmax, lonmin, lonmax,
            lonbinsize, latbinsize
            )
    else:
        raise Exception('Array type no yet supported : %s', values.dtype)
    return result, nbsampl


def __fastregrid_float64(
        np.ndarray[DTYPEf32_t, ndim=1] lat,
        np.ndarray[DTYPEf32_t, ndim=1] lon,
        np.ndarray[DTYPEf64_t, ndim=1] values,
        float latmin=-90.,
        float latmax=90.,
        float lonmin=-180.,
        float lonmax=180.,
        float lonbinsize=1.,
        float latbinsize=1.):
    """
    Place unevenly spaced 2D float64 data on a simple grid by 2D binning
    (mean).

    Args
        lat (ndarray (1D) of numpy float32): the independent data latitude-axis
            of the grid.
        lon (ndarray (1D) of numpy float32): the independent data longitude
            axis of the grid.
        values (ndarray (1D)): the dependent data in the form
            values = f(lat,lon).
        latmin, latmax, lonmin, lonmax (float, optional): the limits of the
            grid on which to remap the input data
        lonbinsize, latbinsize (float, optional): the full width and height of
             each bin on the grid.  If each bin is a cube, then this is the x
             and y dimension.  This is the step in both directions, x and y.
             Defaults to 1.

    Returns
        grid (ndarray (2D)): the evenly gridded data.  The value of each cell
            is the sum value of the contents of the bin.
        bins (ndarray (2D)): A grid the same shape as 'grid', except the value
            of each cell is the number of points in that bin.
    """
    cdef int lonsize = np.round((lonmax - lonmin) / lonbinsize)
    cdef int latsize = np.round((latmax - latmin) / latbinsize)
    cdef np.ndarray[np.float64_t, ndim=2] grid\
        = np.zeros([latsize, lonsize], dtype=np.float64)
    cdef np.ndarray[np.int_t, ndim=2] bins\
        = np.zeros([latsize, lonsize], dtype=np.int)
    cdef unsigned int i
    cdef unsigned int x, y
    cdef double val
    for i, val in enumerate(values):
        if latmin <= lat[<unsigned int> i]\
                and lat[<unsigned int> i] < latmax\
                and lonmin <= lon[<unsigned int> i]\
                and lon[<unsigned int> i] < lonmax:
            y = <unsigned int> np.floor(
                (lat[<unsigned int> i] - latmin) / latbinsize)
            x = <unsigned int> np.floor(
                (lon[<unsigned int> i] - lonmin) / lonbinsize)
            grid[y, x] = grid[y, x] + values[<unsigned int> i]
            bins[y, x] += 1
    return grid, bins


def __fastregrid_float32(
        np.ndarray[DTYPEf32_t, ndim=1] lat,
        np.ndarray[DTYPEf32_t, ndim=1] lon,
        np.ndarray[DTYPEf32_t, ndim=1] values,
        float latmin=-90.,
        float latmax=90.,
        float lonmin=-180.,
        float lonmax=180.,
        float lonbinsize=1.,
        float latbinsize=1.):
    """
    Place unevenly spaced 2D float32 data on a simple grid by 2D binning
    (mean).

    Args
        lat (ndarray (1D) of numpy float32): the independent data latitude-axis
            of the grid.
        lon (ndarray (1D) of numpy float32): the independent data longitude
            axis of the grid.
        values (ndarray (1D)): the dependent data in the form
            values = f(lat,lon).
        latmin, latmax, lonmin, lonmax (float, optional): the limits of the
            grid on which to remap the input data
        lonbinsize, latbinsize (float, optional): the full width and height of
             each bin on the grid.  If each bin is a cube, then this is the x
             and y dimension.  This is the step in both directions, x and y.
             Defaults to 1.

    Returns
        grid (ndarray (2D)): the evenly gridded data.  The value of each cell
            is the sum value of the contents of the bin.
        bins (ndarray (2D)): A grid the same shape as 'grid', except the value
            of each cell is the number of points in that bin.
    """
    cdef int lonsize = np.round((lonmax - lonmin) / lonbinsize)
    cdef int latsize = np.round((latmax - latmin) / latbinsize)
    cdef np.ndarray[np.float32_t, ndim=2] grid\
        = np.zeros([latsize, lonsize], dtype=np.float32)
    cdef np.ndarray[np.int_t, ndim=2] bins\
        = np.zeros([latsize, lonsize], dtype=np.int)
    cdef unsigned int i
    cdef unsigned int x, y
    cdef double val
    for i, val in enumerate(values):
        if latmin <= lat[<unsigned int> i]\
                and lat[<unsigned int> i] < latmax\
                and lonmin <= lon[<unsigned int> i]\
                and lon[<unsigned int> i] < lonmax:
            y = <unsigned int> np.floor(
                (lat[<unsigned int>i] - latmin) / latbinsize)
            x = <unsigned int>np.floor(
                (lon[<unsigned int>i] - lonmin) / lonbinsize)
            grid[y, x] = grid[y, x] + values[<unsigned int>i]
            bins[y, x] += 1
    return grid, bins
