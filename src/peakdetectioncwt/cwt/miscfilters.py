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

import scipy.integrate as integrate
import pylab as pl

def numNonZeroGreaterThanThree(intensity,peakLeft,peakRight):
    numZeros = 0;
    numNotZero = 0;
    epsilon = 0.0001;
    #for (int alpha=peakLeft; alpha<=peakRight;alpha++){
    for alpha in range(peakLeft,peakRight):
        if (intensity[alpha]< epsilon):
            numZeros +=1
        else :
            numNotZero += 1

    if (numNotZero<=3):
        return False

    return True



# checks to see if the number of zero intensity points is greater than or equal to the number of
# non-zero intensity points.
# If there are too many zeros return false wich will mean do not continue with this peak.
def numberOfZeros(intensity,peakLeft,peakRight):
    numZeros = 0;
    numNotZero = 0;
    epsilon = 0.0001;
    #for (int alpha=peakLeft; alpha<=peakRight;alpha++){
    for alpha in range(peakLeft,peakRight):
        if (intensity[alpha]< epsilon):
            numZeros +=1
        else :
            numNotZero += 1

    #print "numZeros: " + str(numZeros)
    #print "numNotZero: " + str(numNotZero)

    if (numZeros>=numNotZero):
        return False

    return True

# compare the ratio of the mean boundry height to the mean of the signal
def peakBoundryPeakHeightComparison(intensity,peakLeft,peakRight):
    meanOfSignal = intensity[peakLeft:peakRight+1].mean()
    meanBoundary = (intensity[peakLeft]+intensity[peakRight])/2.0
    differenceSigBnd = meanOfSignal-meanBoundary
    
    
    if (differenceSigBnd/minimumFeatHeight < 0.2):
        return False

    return True


def coefficientOverAreaCheck(rt,intensity,peakLeft,peakRight,bestCoefficient,coefficientThreshold):

    area = integrate.trapz(intensity[peakLeft:peakRight],x=rt[peakLeft:peakRight])
    #print "area: " + str(area)
    #print "bestCoefficient: " + str(bestCoefficient)
    #print "coefficientThreshold: " + str(coefficientThreshold)

    normedCoef = bestCoefficient/area;
    #print "normedCoef: " + str(normedCoef)
    
    if (normedCoef <coefficientThreshold):
        return False, normedCoef

    return True, normedCoef

def minimumIntensityCheck(intensity,peakLeft,peakRight,minimumInt):
    if max(intensity[peakLeft:peakRight])<minimumInt:
        return False
    return True

def waveletPeakSimilarityCheck(cwtObject,peakLeft,peakRight,bestScale,bestCoeffLocation,waveletPeakSimilarityThresh):

    param_c=0.2
    param_b=pl.exp(-1.0/param_c)
    param_a=(1.0)/(1.0-param_b)
    rt = cwtObject.x

    intensity = cwtObject.signal

    # get the relavant part of the wavelet
    ricker = cwtObject.rickerWavelet(rt-rt[bestCoeffLocation],bestScale)

    # find the apropriate bounds of the wavelet
    maxIndex = pl.argmax(ricker)
    curInt = ricker[maxIndex ]
    rickerLeftIndex = maxIndex 
    i=0
    while curInt > 0.0:
        rickerLeftIndex = maxIndex-i
        i+=1
        if maxIndex-i < 0:
            break
        curInt = ricker[maxIndex-i]

    rickerRightIndex = maxIndex
    i=0
    curInt = ricker[maxIndex ]
    while curInt > 0.0:
        rickerRightIndex = maxIndex+i
        i+=1
        if maxIndex+i >= len(rt):
            break
        curInt = ricker[maxIndex+i]
    

    # add one because slicing ommits last 
    rickerRightIndex +=1

    # get point that is just below zero and adjust appropriatly -> gets you another point and
    # intensity should be only just under zero
    rickerRightIndex +=1
    rickerLeftIndex -=1
    if rickerRightIndex > len(rt):
        rickerRightIndex = len(rt)
    if rickerLeftIndex < 0:
        rickerLeftIndex = 0


    smallestVal = min(ricker[rickerLeftIndex:rickerRightIndex])
    ricker = ricker - smallestVal 

    areaRicker = integrate.trapz(ricker[rickerLeftIndex:rickerRightIndex],x=rt[rickerLeftIndex:rickerRightIndex])
    normedRicker = ricker[rickerLeftIndex:rickerRightIndex]/areaRicker

    area = integrate.trapz(intensity[rickerLeftIndex:rickerRightIndex],x=rt[rickerLeftIndex:rickerRightIndex])
    normedEIC=intensity[rickerLeftIndex:rickerRightIndex]/area

    pl.plot(rt,intensity/area,c='b',marker='o')
    pl.plot(rt[rickerLeftIndex:rickerRightIndex],intensity[rickerLeftIndex:rickerRightIndex]/area,c='r',marker='o')
    pl.plot(rt,ricker ,c='g',marker='o')
    pl.plot(rt[rickerLeftIndex:rickerRightIndex],ricker[rickerLeftIndex:rickerRightIndex]/areaRicker,c='r',marker='o')
    pl.plot(rt[rickerLeftIndex:rickerRightIndex],ricker[rickerLeftIndex:rickerRightIndex],c='r',marker='o')
    pl.show()

    diff_area = integrate.trapz(abs( normedRicker-normedEIC ),x=rt[rickerLeftIndex:rickerRightIndex])

    #cur_similarity = 1-diff_area 
    cur_similarity = (pl.exp(-diff_area/param_c)-param_b)*param_a
    #print "cur_similarity: " + str(cur_similarity)


    #print "normedCoef: " + str(normedCoef)
    
    if (cur_similarity < waveletPeakSimilarityThresh):
        return False, cur_similarity 

    return True, cur_similarity 

