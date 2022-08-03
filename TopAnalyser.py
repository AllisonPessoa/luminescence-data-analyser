#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 31 16:28:18 2022

@author: allison
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from scipy import signal
import peakutils

class Topography():
    def __init__(self, filename, dimension):
        self.dim_x = dimension[0]
        self.dim_y = dimension[1]

        if filename:
            self.topData = np.load(filename, allow_pickle=True)
            self.topData = np.transpose(self.topData)

    @staticmethod
    def load_direct_array(array):
        dim_x = len(array[0])
        dim_y = len(np.transpose(array)[0])
        obj = Topography('', (dim_x, dim_y))
        obj.topData = array.tolist()
        return obj

    def correct_histeresis(self, pixels):
        # Hysteresis correction:
        mask = np.zeros((self.dim_y,self.dim_x))
        for i in range(self.dim_y):
            for j in range(self.dim_x-pixels):
                if i%2 == 0:
                    mask[i][j] = self.topData[i][j+pixels]
                else:
                    mask[i][j] = self.topData[i][j]
        self.topData = np.copy(mask)

    def plot_image(self, limits, fig = None, ax= None, **kwargs):
        if fig == None and ax == None:
            fig, ax =  plt.subplots()
        ax.imshow(self.topData, extent=limits, origin='lower', cmap=cm.get_cmap('gray'), **kwargs)
        return fig, ax

    def remove_line_background(self, start_pixel, stop_pixel):
        mask_bg_corr = np.zeros((self.dim_y,self.dim_x))
        for i in range(self.dim_y):
            start, stop = start_pixel, stop_pixel
            m,b = np.polyfit(range(stop-start), self.topData[i][start:stop], 1)
            mask_bg_corr[i] = [y - (b+m*x) for x,y in zip(range(self.dim_x), self.topData[i])]

        self.topData = np.copy(mask_bg_corr)

    def remove_baseline(self):
        mask_bg_corr = np.zeros((self.dim_y,self.dim_x))
        for i in range(self.dim_y):
            base = peakutils.baseline(self.topData[i], 1)
            mask_bg_corr[i] = np.subtract(self.topData[i], base)

        self.topData = np.copy(mask_bg_corr)
