# coding: utf-8

"""
.. module::posada.Common.slices

Functions to handle slices when reading data from a :package:'Readers' or
:package:'Models' class.

:copyright: Copyright 2011-2016 SOLab.

:copyright: Copyright 2013 Ifremer / Cersat / SOLab.
:license: Released under GPL v3 license, see :ref:'license'.

.. sectionauthor:: Jeff Piolle <jfpiolle@ifremer.fr>, Alexander Myasoedov <mag@rshu.ru>
.. codeauthor:: Jeff Piolle <jfpiolle@ifremer.fr>, Alexander Myasoedov <mag@rshu.ru>
.. author:: SOLab Dev Team
"""


import numpy as np
from ..Common.logger import get_logger
logger = get_logger()


# TODO : method shape() from Slices object should be used instead !
def get_shape_from_slice(slices):
    """Returns the shape of the result from a slicing operation.

    Args:
        slices (list<slice>): a list of slices

    Returns:
        tuple<int>: the shape of the expected result
    """
    shape = []
    for slc in slices:
        shape.append((slc.stop - slc.start) / slc.step)
    return tuple(shape)


# TODO : Slices object creation should be used instead !
def get_nice_slices(slices, dimsizes):
    """Return filled slices with valid positive limits.
    
    Fill (replace None values with numbers). Trim the slices beyond the array
    dimensions. Transform negative slices into positive values.
    """
    newslices = []
    if slices is None:
        for d in dimsizes:
            newslices.append(slice(0, d, 1))
        return newslices
    for sli, dimsize in zip(slices, dimsizes):
        start, stop = sli.start, sli.stop
        if start is not None and start < 0:
            start += dimsize
        elif start is None:
            start = 0
        if stop is not None and stop < 0:
            stop += dimsize
        elif stop is None:
            stop = dimsize
        if start is not None and start > dimsize:
            start = dimsize
        if stop is not None and stop > dimsize:
            stop = dimsize
        step = sli.step
        if step is None:
            step = 1
        newslices.append(slice(start, stop, step))   
    return newslices    
    

# TODO : Slices object creation should be used instead !
def fill_slices(slices, dimsizes):
    """Fill (replace None values with numbers) and check slices with respect
    to the dimensions of an array.    
    """
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
        if (step > 0 and start >= stop) or (step < 0 and start <= stop):
            raise Exception('Unexpected slice start Vs stop : {}'
                            .format(sli))
        # Special case : negative step until the first element
        if stop == -1:
            stop = None
        filled_slices.append(slice(start, stop, step))
    return filled_slices


# TODO : method absolute_slices() from Slices object should be used instead !
def get_absolute_slices(view, slices, dimnames, dimsizes):
    """Returns absolute slices from the slices relative to a view."""
    # view slicing
    if view is not None:
        # make sure viewslice info is complete
        viewslices = format_slices(view, dimnames)
        viewslices = fill_slices(viewslices, dimsizes)
    # standard slicing
    if slices is not None:
        if view is None:
            finalslices = fill_slices(slices, dimsizes)
        else:
            # slices are relative to the opened view
            viewdimsizes = dimsizes
            finalslices = fill_slices(slices, viewdimsizes)
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
        finalslices = fill_slices(finalslices, dimsizes)
    return finalslices


# TODO : method shape() from Slices object should be used instead !
def adjust_dimsize(view, dimname, dimsize):
    """correct the size of a dimension with respect to a view"""
    if view is not None and dimname in view:
        viewslice = fill_slices([view[dimname]],
                                [dimsize])[0]
        start, stop, step = viewslice.start, viewslice.stop, viewslice.step
        if stop is None:
            stop = -1
        dimsize = 1 + (abs(stop - start) - 1) / abs(step)
    return dimsize


def format_slices(slices, fielddims):
    """format the slices to an explicit list matching the list and order of
    a field dimensions

    the slices as provided as a dictionary where keys are the dimensions
    over which subsetting needs to be performed. Implicit dimensions where
    no subsetting is requested don't need to be specified. The completion
    to a comprehensive ordered list of slices for all field dimensions is
    performed by this function.

    Args:
        slices (dict): keys are the names and values are the slice object
            of the dimensions in which to subset

        fielddims (list): full list of the field dimension names

    Returns:
        list: a comprehensive list of slice objects for all field
            dimensions.
    """
    if slices is None:
        return None
    formatted_slices = []
    dimensions = fielddims
    for _ in xrange(len(dimensions)):
        formatted_slices.append(slice(None, None, None))
    for k in slices.keys():
        try:
            idim = list(dimensions).index(k)
            formatted_slices[idim] = slices[k]
        except ValueError:
            pass
    return formatted_slices


def format_indices(indices, fielddims):
    """format the indices to an explicit list of 'slice' matching the list
    and order of a field dimensions

    the indices are provided as a dictionary where keys are the dimensions
    over which subsetting needs to be performed. Implicit dimensions where
    no subsetting is requested don't need to be specified. The completion
    to a comprehensive ordered list of slices for all field dimensions is
    performed by this function.

    Args:
        indices (dict): keys are the names and values are the indices
            of the dimensions in which to subset

        fielddims (list): full list of the field dimension names

    Returns:
        list: a comprehensive list of slice objects for all field
            dimensions.
    """
    if indices is None:
        return None
    formatted_indices = []
    dimensions = fielddims
    for i in xrange(len(dimensions)):
        formatted_indices.append(slice(None, None, None))
    for k in indices.keys():
        try:
            idim = list(dimensions).index(k)
            formatted_indices[idim] = indices[k]
        except ValueError:
            pass
    return tuple(formatted_indices)


class Slices(tuple):
    """Class defining a Slices object.

    A Slices object is a tuple of slice objects and therefore can be used
    to slice a numpy array or a netCDF4 variable for instance. During
    initialisation, input slice objects are checked and filled according to
    the dimension sizes of the variable to slice. It also provides some methods
    which aim at making the use of slices in posada.easier.

    Args:
        slices (None, slice, int, tuple<None,slice,int>): input slice object(s).

        dimsizes (None, int, list<int>, tuple<int>): dimension size(s) of the
            variable to slice.

    Example:
        >>> values = numpy.arange(10)
        >>> slices = Slices(slice(1, None, 4), values.shape)
        >>> print slices
        (slice(1, 10, 4),)
        >>> print values[1::4]
        [1 5 9]
        >>> print values[slices]
        [1 5 9]

    """

    def __new__(cls, slices=None, dimsizes=None):
        """
        """
        # Check dimsizes input type and copy it as a list
        if dimsizes is None:
            _dimsizes = []
        elif isinstance(dimsizes, int):
            _dimsizes = [dimsizes]
        elif isinstance(dimsizes, (list, tuple)):
            _dimsizes = list(dimsizes)
        else:
            raise Exception('Unexpected type for dimsizes.')
        # Check slices input type
        if slices is None:
            _slices = [slice(None)] * len(_dimsizes)
        elif isinstance(slices, (slice, int)):
            _slices = [slices]
        elif isinstance(slices, tuple):
            _slices = slices
        elif isinstance(slices, list):
            logger.warning('slices as list are deprecated, it will raise an '
                            'error in the future -> please use a tuple.')
            _slices = slices
        else:
            raise Exception('Unexpected type for slices.')
        # Check length of slices and dimsizes
        if len(_slices) != len(_dimsizes):
            raise Exception('Mismatch between slices/dimsizes lengths.')
        # Check and fill slices according to dimsizes
        # Warning : if the rules defining an acceptable slice are changed,
        # then some methods (eg shape() or absolute_slices()) may give false results.
        filled_slices = []
        _is_index = []
        for _slice, _dimsize in zip(_slices, _dimsizes):
            # Check _dimsize type
            if not isinstance(_dimsize, int):
                raise Exception('Unexpected type for dimsizes.')
            # Check _slice type and format it to a slice object if needed
            is_index = False
            if _slice is None:
                _slice = slice(None)
            elif isinstance(_slice, slice):
                pass
            elif isinstance(_slice, int):
                _slice = slice(_slice, _slice + 1, 1)
                is_index = True
            else:
                raise Exception('Unexpected type for slices.')
            _is_index.append(is_index)
            # Check _slice start/stop/step
            if _slice.start is not None and abs(_slice.start) > _dimsize:
                raise Exception('Unexpected slice start : {}'.format(_slice))
            if _slice.stop is not None and abs(_slice.stop) > _dimsize:
                raise Exception('Unexpected slice stop : {}'.format(_slice))
            if _slice.step is not None and _slice.step == 0:
                raise Exception('Unexpected slice step : {}'.format(_slice))
            # Fill slice
            start, stop, step = _slice.indices(_dimsize)
            if (step > 0 and start >= stop) or (step < 0 and start <= stop):
                raise Exception('Unexpected slice start Vs stop : {}'.format(_slice))
            # Special case : negative step until the first element
            # here : var[_slice] != var[slice(*_slice.indices(_dimsize))]
            if stop == -1:
                stop = None
            filled_slices.append(slice(start, stop, step))
        instance = super(Slices, cls).__new__(cls, tuple(filled_slices))
        instance._dimsizes = _dimsizes
        instance._is_index = _is_index
        return instance

    def __add__(self, slices):
        """
        """
        if not isinstance(slices, Slices):
            raise Exception('slices is not a Slices object.')
        added_slices = super(Slices, self).__add__(slices)
        added_dimsizes = self._dimsizes + slices._dimsizes
        new_slices = Slices(added_slices, added_dimsizes)
        new_slices._is_index = self._is_index + slices._is_index
        return new_slices

    def shape(self):
        """Returns the expected shape of the sliced variable.

        Returns:
            tuple: numpy array like shape.

        Example:
            >>> values = numpy.zeros((10, 5))
            >>> slices = Slices([slice(None, None, 3), slice(None, None, -2)], values.shape)
            >>> print values[::3, ::-2].shape
            (4, 3)
            >>> print values[slices].shape
            (4, 3)
            >>> print slices.shape()
            (4, 3)

        """
        shp = []
        for _slice in self:
            start, stop, step = _slice.start, _slice.stop, _slice.step
            if stop is None:
                stop = -1
            shp.append(1 + (abs(stop - start) - 1) / abs(step))
        return tuple(shp)

    def reduced_shape(self):
        """Returns the reduced shape of the sliced variable.

        During instantiation of this object, all int (integer indexing) are
        converted to simple slice because it is more convenient to deal only
        with slices objects. But as a consequence, after slicing, none of the
        dimensions are squeezed as it should occur with integer indexing.

        This method computes the shape by taking care of squeezing all the
        dimensions where integer indexing was initially specified. Reshaping
        the sliced variable with the result of this method allows to get back
        to classical indexing.

        Returns:
            tuple: numpy array like shape.

        Example:
            >>> values = numpy.arange(50).reshape((5, 10))
            >>> slices = Slices((2, slice(None, None, 3)), values.shape)
            >>> print slices.shape()
            (1, 4)
            >>> print slices.reduced_shape()
            (4,)
            >>> print values[slices]
            [[20 23 26 29]]
            >>> print values[slices].reshape(slices.reduced_shape())
            [20 23 26 29]

        """
        shp = []
        for _slice, _is_index in zip(self, self._is_index):
            if _is_index == True:
                continue
            start, stop, step = _slice.start, _slice.stop, _slice.step
            if stop is None:
                stop = -1
            shp.append(1 + (abs(stop - start) - 1) / abs(step))
        return tuple(shp)

    def absolute_slices(self, view_slices):
        """Returns absolute slices assuming this Slices object is relative to a view.

        It can be seen as slices composition : slicing with the absolute
        slices is equivalent to slice with the view slices then with the
        relative to view slices (i.e. absolute = view o relative).

        Args:
            view_slices (Slices): the view slices as a Slices object.

        Returns:
            Slices: the absolute slices as a new Slices object.

        Example:
            If a view is defined in a posada.Readers, the slices given by the user
            to the.Readers read_values() method are relative to this view. With this
            method, we get the right slices to apply to original data.

            >>> values = numpy.arange(100).reshape((2, 50))
            >>> view_slices = Slices([slice(1, 2), slice(None, None, -4)], values.shape)
            >>> print view_slices
            (slice(1, 2, 1), slice(49, None, -4))
            >>> rel_slices = Slices([slice(0, 1), slice(None, None, -2)], view_slices.shape())
            >>> print rel_slices
            (slice(0, 1, 1), slice(12, None, -2))
            >>> abs_slices = rel_slices.absolute_slices(view_slices)
            >>> print abs_slices
            (slice(1, 2, 1), slice(1, 50, 8))
            >>> print values[1:2, ::-4][0:1, ::-2]
            [[51 59 67 75 83 91 99]]
            >>> print values[view_slices][rel_slices]
            [[51 59 67 75 83 91 99]]
            >>> print values[abs_slices]
            [[51 59 67 75 83 91 99]]

        """
        if not isinstance(view_slices, Slices):
            raise Exception('view_slices is not a Slices object.')
        if len(self) != len(view_slices.reduced_shape()):
            raise Exception('Mismatch between slices/viewslices lengths.')
        abs_slices = []
        self_position = 0
        for view_sli, is_view_index in zip(view_slices, view_slices._is_index):
            if is_view_index:
                abs_slices.append(view_sli.start)
            else:
                # self slice (ie slice relative to view slice)
                sli = self[self_position]
                first = sli.start
                if sli.stop is None:
                    last = 0
                else:
                    last = sli.stop - sli.step / abs(sli.step)
                step = sli.step
                # view slice
                view_first = view_sli.start
                view_step = view_sli.step
                # absolute slice
                abs_first = view_first + first * view_step
                abs_last = view_first + last * view_step
                abs_start = abs_first
                abs_step = view_step * step
                abs_stop = abs_last + abs_step / abs(abs_step)
                if abs_stop == -1:
                    abs_stop = None
                if self._is_index[self_position]:
                    abs_slices.append(abs_start)
                else:
                    abs_slices.append(slice(abs_start, abs_stop, abs_step))
                self_position += 1
        return Slices(tuple(abs_slices), view_slices._dimsizes)

    def indices_array(self, dtype=None):
        """Returns arrays of indices corresponding to the slices.

        Args:
            dtype (dtype, optional): the type of the output indices arrays.

        Returns:
            list<ndarray>: indices array for each dimension.

        Example:
            >>> slices = Slices([slice(2, 12, 3), slice(9, 2, -2)], [15, 12])
            >>> print slices.indices_array()
            [array([ 2,  5,  8, 11]), array([9, 7, 5, 3])]

        """
        ind_arr = []
        for _slice in self:
            stop = _slice.stop
            if stop is None:
                stop = -1
            ind_arr.append(np.arange(_slice.start, stop, _slice.step,
                                     dtype=dtype))
        return ind_arr

    def indices_minmax(self):
        """Returns minimum and maximum of indices corresponding to the slices.

        Returns:
            list<list>: indices [min, max] for each dimension.

        Example:
            >>> slices = Slices([slice(2, 12, 3), slice(9, 2, -2)], [15, 12])
            >>> print slices.indices_minmax()
            [[2, 11], [3, 9]]

        """
        ind_minmax = []
        for _slice, size in zip(self, self.shape()):
            stopin = _slice.start + (size - 1) * _slice.step
            if _slice.step > 0:
                ind_minmax.append([_slice.start, stopin])
            else:
                ind_minmax.append([stopin, _slice.start])
        return ind_minmax

    def indices_offsetextent(self):
        """Returns offset and extent of indices corresponding to the slices.

        Returns:
            list<list>: indices [offset, extent] for each dimension.

        Example:
            >>> slices = Slices([slice(2, 12, 3), slice(9, 2, -2)], [15, 12])
            >>> print slices.indices_offsetextent()
            [[2, 10], [3, 7]]

        """
        ind_offext = []
        for minmax in self.indices_minmax():
            ind_offext.append([minmax[0], minmax[1] - minmax[0] + 1])
        return ind_offext

    def outsub_insub_slices(self, sub_offset, sub_size):
        """Returns new Slices objects to deal with a sub area of a variable.

        In 2D, for example, lets's assume an image is made up of subimages,
        each one located with an offset and a size (for each dimension). And
        let's assume also we want to build a sliced version of this image given
        this Slices object and given the subimages. Then, this method returns
        for a given subimage : an outsub Slices object defining where to fill
        the sliced image we want to build and an insub Slices object defining
        where to read in the subimage.

        Args:
            sub_offset (list<int>): offset of the sub area for each dimension.

            sub_size (list<int>): size of the sub area for each dimension.

        Returns:
            list<Slices>, list<None>: [outsub, insub] as new Slice objects or
                [None, None] if the sub area does not overlap with the slices
                area.

        Example:
            With a 4x4 image made up of 4 2x2 subimages.

            >>> # Subimages definition
            >>> sub_offsets = [[0, 0], [0, 2], [2, 0], [2, 2]]
            >>> sub_sizes = [[2, 2], [2, 2], [2, 2], [2, 2]]
            >>> sub_images = [numpy.ones(sub_sizes[i]) * i for i in range(4)]
            >>> # Image (that we don't want to build but shown for the example)
            >>> _image = numpy.array([[0., 0., 1., 1.], [0., 0., 1., 1.],
            ...                       [2., 2., 3., 3.], [2., 2., 3., 3.]])
            >>> print _image
            [[ 0.  0.  1.  1.]
             [ 0.  0.  1.  1.]
             [ 2.  2.  3.  3.]
             [ 2.  2.  3.  3.]]
            >>> # Sliced image
            >>> slices = Slices([slice(0, 2), slice(1, 3)], [4, 4])
            >>> sliced_image = numpy.zeros(slices.shape())
            >>> for i in range(4):
            ...     outsub, insub = slices.outsub_insub_slices(sub_offsets[i], sub_sizes[i])
            ...     if outsub is not None:
            ...         sliced_image[outsub] = sub_images[i][insub]
            ...     print outsub, insub
            ...
            (slice(0, 2, 1), slice(0, 1, 1)) (slice(0, 2, 1), slice(1, 2, 1))
            (slice(0, 2, 1), slice(1, 2, 1)) (slice(0, 2, 1), slice(0, 1, 1))
            None None
            None None
            >>> print sliced_image
            [[ 0.  1.]
             [ 0.  1.]]
            >>> print _image[slices]
            [[ 0.  1.]
             [ 0.  1.]]

        """
        if len(self) != len(sub_offset) or len(self) != len(sub_size):
            raise Exception('Mismatch between slices/offsets/sizes lengths.')
        outsub_slices = []
        insub_slices = []
        shp = self.shape()
        minmax = self.indices_minmax()
        for i in range(len(self)):
            # self slice (ie absolute slice)
            absmin = minmax[i][0]
            #absmax = minmax[i][1]
            abssize = shp[i]
            absstep = self[i].step
            # sub
            submin = sub_offset[i]
            subsize = sub_size[i]
            submax = submin + subsize - 1
            # Find outsub min and max
            outmin = np.ceil(float(submin - absmin) / abs(absstep))
            outmax = np.floor(float(submax - absmin) / abs(absstep))
            if outmax < 0 or outmin >= abssize:
                return None, None
            outmin = int(np.maximum(0, outmin))
            outmax = int(np.minimum(abssize - 1, outmax))
            # Check outmin / outmax
            assert outmin >= 0 and outmin < abssize and \
                outmax >= 0 and outmax < abssize
            # Find insub min and max
            inmin = int(absmin + outmin * abs(absstep) - submin)
            inmax = int(absmin + outmax * abs(absstep) - submin)
            # Check inmin / inmax
            assert inmin >= 0 and inmin < subsize and \
                inmax >= 0 and inmax < subsize
            # Make slices
            if self._is_index[i]:
                assert outmin == outmax and inmin == inmax
                outsub_slices.append(outmin)
                insub_slices.append(inmin)
            elif absstep > 0:
                outsub_slices.append(slice(outmin, outmax+1, 1))
                insub_slices.append(slice(inmin, inmax+1, absstep))
            else:
                outsub_slices.append(slice(abssize - 1 - outmax,
                                           abssize - outmin, 1))
                if inmin == 0:
                    instop = None
                else:
                    instop = inmin - 1
                insub_slices.append(slice(inmax, instop, absstep))
        res = (Slices(tuple(outsub_slices), shp),
               Slices(tuple(insub_slices), sub_size))
        # Check shapes
        assert res[0].shape() == res[1].shape()
        return res

