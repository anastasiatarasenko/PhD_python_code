# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 12:08:27 2017

@author: atarasen
"""
import numpy as np

def compute_window_mean_and_var_strided(image, window_w, window_h):
    w, h = image.shape
    strided_image = np.lib.stride_tricks.as_strided(image, 
                                                    shape=[w - window_w + 1, h - window_h + 1, window_w, window_h],
                                                    strides=image.strides + image.strides)
    # important: trying to reshape image will create complete 4-dimensional compy
    means = strided_image.mean(axis=(2,3)) 
    mean_squares = (strided_image ** 2).mean(axis=(2, 3)) 
    maximums = strided_image.max(axis=(2,3))
    
    variations = mean_squares - means ** 2
    return means, maximums, variations 

def compute_window_mean_and_var_strided(image, window_w, window_h):
    w, h = image.shape
    strided_image = np.lib.stride_tricks.as_strided(image, 
                                                    shape=[w - window_w + 1, h - window_h + 1, window_w, window_h],
                                                    strides=image.strides + image.strides)
    # important: trying to reshape image will create complete 4-dimensional compy
    means = strided_image.mean(axis=(2,3)) 
    mean_squares = (strided_image ** 2).mean(axis=(2, 3)) 
    maximums = strided_image.max(axis=(2,3))
    
    variations = mean_squares - means ** 2
    return means, maximums, variations

image = np.random.rand((10, 5),)
