# coding: utf-8

"""
.. module::posada.Common.wind

This module contains the main functions for
WIND SPEED and WIND DIRECTION calculations.

:copyright: Copyright 2011-2016 SOLab.

.. sectionauthor:: Alexander Myasoedov <mag@rshu.ru>
.. codeauthor:: Alexander Myasoedov <mag@rshu.ru>
.. author:: SOLab Dev Team
"""

from .imports import zeros, atan, pi, reshape, size, sqrt
from numpy.ma import masked_array


def windDirection(u, v):
    U = u.ravel()
    V = v.ravel()
    direction = zeros(size(U))
    for i in range(0, len(U)):
        if U[i] >= 0 and V[i] > 0:
            direction[i] = ((180 / pi) * atan(abs(U[i] / V[i])) + 180)
        if U[i] < 0 and V[i] > 0:
            direction[i] = (-(180 / pi) * atan(abs(U[i] / V[i])) + 180)
        if U[i] >= 0 and V[i] < 0:
            direction[i] = (-(180 / pi) * atan(abs(U[i] / V[i])) + 360)
        if U[i] < 0 and V[i] < 0:
            direction[i] = ((180 / pi) * atan(abs(U[i] / V[i])))
        if V[i] == 0 and U[i] > 0:
            direction[i] = 270
        if V[i] == 0 and U[i] < 0:
            direction[i] = 90
        if V[i] == 0 and U[i] == 0:
            direction[i] = 0

    direction = reshape(direction, v.shape)
    return masked_array(data=direction, mask=u.mask)


def windSpeed(u, v):
    return sqrt(u**2 + v**2)
