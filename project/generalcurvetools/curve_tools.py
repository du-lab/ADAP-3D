#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2016-2017, Du-lab
# Author: Owen Myers
# Contact: omyers2@uncc.edu
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>.

import curves
from scipy.optimize import curve_fit
import pylab as pl
import scipy.integrate as integrate
import random

def asymmetric_gaussian_fit(peakLoc,peakIntArr,visualize):
    """
    Does an asymmetric gaussian fit and then uses the local similarity measure (see below) to measure the similarity
     between the input curve and the fit. Just one way to see how good the fit is.
    :param peakLoc: Index of peak maxima
    :param peakIntArr: Intensity array of cuve/signal/peak
    :param visualize: boolean. Plot the fit and the curve?
    :return: 
    """

    max_y = peakIntArr[peakLoc]
    x = pl.arange(len(peakIntArr))

    try:
        popt, pcov = curve_fit(curves.asym_gaussian,
                x,
                peakIntArr,
                bounds=([max_y - max_y*.15, peakLoc - 1, 1.0 / 2.35482,             0],
                          [max_y + max_y*.15, peakLoc + 1, 2*len(peakIntArr) / 2.35482,50]))

    except RuntimeError,e:
        return [0,0,0,0], 0



    cur_sim = get_similarity(peakIntArr,curves.asym_gaussian(x,popt[0],popt[1],popt[2],popt[3]))

    if visualize:
        print "cur similarity: " + str(cur_sim )
        pl.plot(x,curves.asym_gaussian(x,popt[0],popt[1],popt[2],popt[3]),c='r')
        pl.plot(x,peakIntArr,c='b')
        pl.show()

    return popt, cur_sim


def get_similarity(intensity1,intensity2):
    """
    This similarity measure looks at the area between two normalized curves.
    :param intensity1: Curve 1
    :param intensity2: Curve 2
    :return: Similarity will always be a value between 0(high similarity) and 1(low sim).
    """
    param_c=0.2
    param_b=pl.exp(-1.0/param_c)
    param_a=(1.0)/(1.0-param_b)

    x = pl.arange(0,len(intensity1))
    area = integrate.trapz(intensity1,x=x)
    normed1 = intensity1/area

    area2 = integrate.trapz(intensity2,x=x)
    normed2 = intensity2/area2

    diff_area = integrate.trapz(abs(normed1-normed2),x=x)

    cur_similarity = (pl.exp(-diff_area/param_c)-param_b)*param_a

    return cur_similarity

def estimate_fwhm_ms(int_matrix,numberOfScansForFWHMCalc):
    """
    This function is going to randomly select a variety of spectra and estimate a FWHM using a
    gaussian fit on some of the most intense peaks in the selected spectra.
    
    :param int_matrix: Matrix of all data points
    :param numberOfScansForFWHMCalc: how many scans will we use for estimate
    :return: estimate of full-width-half-max in number of points
    """

    # keep the values to be averaged over
    fwhm_arr = pl.array([])

    print "Estimating FWHM of profile MS"
    print "..."
    for i in range(numberOfScansForFWHMCalc):
        # randomly select a scan to investigate.
        spec_index = random.randint(0,len(pl.squeeze(pl.asarray(int_matrix[0,:].todense())))-1)
        cur_spec = pl.squeeze(pl.asarray(int_matrix[:,spec_index].todense()))
        # could be 0 intensity points in good spectra becuase of the way the int_matrix is bui
        # lets work in units of scan number
        x = pl.arange(0,len(cur_spec ))
        # Only going to be looking at the highest peak in each scan. Get the location of that peak.
        max_y = max(cur_spec )
        scan_max_y = pl.where(cur_spec == max_y)[0][0]
        popt, pcov = curve_fit(curves.gaussian,
                x,
                cur_spec,
                bounds = ([max_y * 0.8, scan_max_y - 1, 0.25 / 2.35482],
                          [max_y * 1.2, scan_max_y + 1, 5 / 2.35482]))
        tmp_x = x
        tmp_y = curves.gaussian(x,popt[0],popt[1],popt[2])

        #pl.plot(tmp_x[scan_max_y-20:scan_max_y+20],tmp_y[scan_max_y-20:scan_max_y+20],c='r')
        #pl.plot(x[scan_max_y-20:scan_max_y+20],cur_spec[scan_max_y-20:scan_max_y+20],c='b',marker='o')
        #pl.show()

        # with the right fit parameter you can determine the full width at half max.
        # 2.35... is a known value to take this parameter and get the FWHM with.
        fwhm = 2.35482*popt[2]
        fwhm_arr = pl.append(fwhm_arr, fwhm )

    avg_fwhm = fwhm_arr.mean()
    print "Done estimating FWHM of profile MS"
    print ("avg_fwhm: " + str(avg_fwhm ))

    # round up
    return int(avg_fwhm+0.5)