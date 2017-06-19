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

import numpy as np

def findBestEstimateOfSignalToNoiseRatio(intensity,peakLeft,peakRight):

    curSN1 = filterBySNWindowInOutSweep(intensity,peakLeft,peakRight)
    curSN2 = filterBySNStaticWindowSweep(intensity,peakLeft,peakRight)
    curSN = max(curSN1,curSN2)

    return curSN


# This sweeps the window size ffrom the points furthest from the peak out and then from the points closest to the peak
# out twords the furthest.
def filterBySNWindowInOutSweep(intensities,peakLeft,peakRight):

    toFindMinNoPeakNoise = np.array([])
    toFindMinWithPeakNoise = np.array([])
    toFindLocalMean = np.array([])
    
    peakHeight = max(intensities[peakLeft:peakRight])
    meanOfSignal = intensities[peakLeft:peakRight].mean()
    
    peakWidth = peakRight-peakLeft
    
    # Initial size of expanding window will be twice the peak width on either side of the peak
    # (really half the total window size.
    initialWindowSize = 2*peakWidth
    # final size of expanding window.
    finalWindowSize = 8*peakWidth
    
    
    toFindMinNoPeakNoise = np.array([])
    toFindMinWithPeakNoise = np.array([]) 
    
    # loop over different window sizes
    #for (int i = 0; i<(finalWindowSize+1); i++){
    for i in range(finalWindowSize+1):
        stdDevStats = np.array([])
        curRight = peakRight + initialWindowSize + i+1
        curLeft = peakLeft - initialWindowSize - i-1
        
        this_time_num_added = 0
        
        # Make sure the these values don't go out of bounds
        if (curRight>=len(intensities)):
            curRight = len(intensities) - 1
    
        if (curLeft<0):
            curLeft = 0
        
        if (not(abs(curLeft-peakLeft)<initialWindowSize)):
            #for (int j = curLeft; j < peakLeft;j++){
            for j in range(curLeft,peakLeft):
                stdDevStats = np.append(stdDevStats,intensities[j])
                this_time_num_added+=1
    
        if (not (abs(curRight-peakRight)<initialWindowSize)):
            #for (int j = curRight; j > peakRight ;j--){
            for j in range(curRight,peakRight):
                stdDevStats = np.append(stdDevStats,intensities[j])
                this_time_num_added+=1
    
        if (this_time_num_added == 0):
            continue
        
        
        # first find the standard deviation from the windows EXCLUDING the peak
        curNoPeakNoise = stdDevStats.std()
        toFindMinNoPeakNoise = np.append(toFindMinNoPeakNoise,curNoPeakNoise)
        curNoPeakMean = stdDevStats.mean()
        toFindLocalMean = np.append(toFindLocalMean,curNoPeakMean)
        
        # Now add in the peak becausase the difference in std with the peak will tell us something about it height.
        #for (int j = peakLeft; j < (peakRight+1); j++){
        for j in range(peakLeft,(peakRight+1)):
            stdDevStats = np.append(stdDevStats,intensities[j])
    
        curWithPeakNoise = stdDevStats.std()
        toFindMinWithPeakNoise = np.append(toFindMinWithPeakNoise,curWithPeakNoise)
    
    # Now sweep the points closest to the peak out
    anchorRight = peakRight + initialWindowSize + finalWindowSize
    anchorLeft = peakRight - initialWindowSize - finalWindowSize
    # Make sure the these values don't go out of bounds
    if (anchorRight>=len(intensities)):
        anchorRight = len(intensities) - 1
    
    if (anchorLeft<0):
        anchorLeft = 0
        
    #for (int i = 0; i<(finalWindowSize-initialWindowSize); i++ ){
    for i in range(finalWindowSize-initialWindowSize):
        stdDevStats = np.array([])
        
        curRight = peakRight+1 + i
        curLeft = peakLeft-1 - i
        
        this_time_num_added = 0
        
        # Make sure the these values don't go out of bounds
        if (curRight>=len(intensities)):
            curRight = len(intensities) - 1
    
        if (curLeft<0):
            curLeft = 0
    
        
        if (not (abs(curLeft-anchorLeft)<initialWindowSize)):
            #for (int j = curLeft; j > anchorLeft;j--){
            for j in range(curLeft,anchorLeft,-1):
                stdDevStats = np.append(stdDevStats,intensities[j])
                this_time_num_added+=1
    
        if (not (abs(curRight-anchorRight)<initialWindowSize)):
            #for (int j = curRight; j < anchorRight ;j++){
            for j in range(curRight,anchorRight):
                stdDevStats = np.append(stdDevStats,intensities[j])
                this_time_num_added+=1
    
        if (this_time_num_added ==0):
            continue
        
        curNoPeakNoise = stdDevStats.std()
        toFindMinNoPeakNoise = np.append(toFindMinNoPeakNoise,curNoPeakNoise)
    
    # hadle the posibility that there is not enough data to get the SN ratio. With the
    # "if (not (abs(curRight-anchorRight)<initialWindowSize))" statements you will not get anything if
    # the peakwidth is comparable to the length of the entire EIC
    if len(toFindMinNoPeakNoise)==0:
        return 0.0

    # for a good looking peak we think the difference between the noise with and without the peak included should
    # be about a factor of 2. Specificaly with the peak the noise should be 2x higher than without if it is real.
    bestNoPeakNoise = min(toFindMinNoPeakNoise)
    
    #if (meanOfSignal>(1.5*smallestLocalMean)){
    #    
    #}
    #if (bestWithPeakNoise/bestNoPeakNoise <= 1.5){
    #    return false;
    #}
    
    # before calculating the signal to noise ratio we need to "normalize" the height.
    # What can happen is if a bad peak is found on a plateau the local stadard deviation will be 
    # small compared to the absolute height (or mean of peak) of peak.
    # 
    # Find the minimum intesity point starting at the boundry and going a boundry width out.
    # Then subtract this from the peak value (height or mean) before comparing to the standard
    # deviation.
    
    stdDevStats=np.array([])
    rightBound = peakRight+peakWidth
    leftBound = peakLeft-peakWidth
    
    # Make sure the these values don't go out of bounds
    if (rightBound>=len(intensities)):
        rightBound = len(intensities) - 1
    
    if (leftBound<0):
        leftBound = 0
    
    # out to right
    #for (int i = peakRight; i <= rightBound; i++):
    for i in range(peakRight,rightBound+1):
        stdDevStats = np.append(stdDevStats,intensities[i])
    
    smallIntensity1 = min(stdDevStats)
    stdDevStats = np.array([])
    #for (int i = leftBound; i <= peakLeft; i++){
    for i in range(leftBound,peakLeft+1):
        stdDevStats = np.append(stdDevStats,intensities[i])
    
    smallIntensity2 = min(stdDevStats)
    smallIntensityAvg = (smallIntensity1+smallIntensity2)/2.0
    
    
    if bestNoPeakNoise<1E-13:
        SNRatio = float("inf")
    else:
        SNRatio = (peakHeight-smallIntensityAvg)/bestNoPeakNoise
    #SNRatio = (meanOfSignal-smallIntensityAvg)/bestNoPeakNoise;
    
    return SNRatio;
                


# Uses two windows on eaiter side of the peak that are a constand width. These windows then both
# slide out to a set distance ithout changing side. The smallest standard deviation is taken to be 
# the noise.
def filterBySNStaticWindowSweep( intensities,peakLeft,peakRight):
    stdDevStats = np.array([])
    toFindMinNoPeakNoise = np.array([])
    
    # first get the value of the peak
    #for (int i=peakLeft; i<(peakRight+1); i++){
    for i in range(peakLeft,peakRight+1):
        stdDevStats = np.append(stdDevStats,intensities[i])
    
    peakHeight = max(stdDevStats)
    meanOfSignal = stdDevStats.mean()
    
    peakWidth = peakRight-peakLeft
    
    # Initial size of expanding window will be twice the peak width on either side of the peak
    # (really half the total window size.
    windowSize = 2*peakWidth
    # final size of expanding window.
    furthestPoint = 8*peakWidth
    
    finalCurLeft = 0
    finalCurRight = 0
    stop_sliding_right = False
    stop_sliding_left = False

    #for (int i = 0; i < (furthestPoint-windowSize);i++){
    for i in range(furthestPoint-windowSize):
        stdDevStats = np.array([])
        
        curRight = peakRight + i
        curRightRight = curRight + windowSize
        
        if ((curRightRight>=len(intensities))and(not stop_sliding_right)):
            stop_sliding_right=True
            finalCurRight = curRight-1

        if (stop_sliding_right):
            curRight = finalCurRight
            #curRightRight = curRightRight - windowSize
            curRightRight = len(intensities) - 1
        
        curLeft = peakLeft - i
        curLeftLeft = curLeft - windowSize
        
        if ((curLeftLeft<0) and (not stop_sliding_left)):
            stop_sliding_left=True
            # revert to last value
            finalCurLeft = curLeft+1

        if (stop_sliding_left):
            curLeft = finalCurLeft
            #curLeftLeft = curLeftLeft + windowSize
            curLeftLeft = 0
            
        
        #for (int j = curRight; j<=curRightRight; j++){
        for j in range(curRight,curRightRight+1):
            stdDevStats = np.append(stdDevStats,intensities[j])

        #for (int j = curLeft; j>=curLeftLeft; j--){
        for j in range(curLeft,curLeftLeft-1,-1):
            stdDevStats = np.append(stdDevStats,intensities[j])

        curNoPeakNoise = stdDevStats.std()
        
        toFindMinNoPeakNoise = np.append(toFindMinNoPeakNoise,curNoPeakNoise)

    bestNoPeakNoise = min(toFindMinNoPeakNoise)
    # before calculating the signal to noise ratio we need to "normalize" the height.
    # What can happen is if a bad peak is found on a plateau the local stadard deviation will be 
    # small compared to the absolute height (or mean of peak) of peak.
    # 
    # Find the minimum intesity point starting at the boundry and going a boundry width out.
    # Then subtract this from the peak value (height or mean) before comparing to the standard
    # deviation.
    
    stdDevStats = np.array([])
    rightBound = peakRight+peakWidth
    leftBound = peakLeft-peakWidth
    
    # Make sure the these values don't go out of bounds
    if (rightBound>=len(intensities)):
        rightBound = len(intensities) - 1

    if (leftBound<0):
        leftBound = 0

    # out to right
    #for (int i = peakRight; i <= rightBound; i++){
    for i in range(peakRight,rightBound+1):
        stdDevStats = np.append(stdDevStats,intensities[i])
        
    smallIntensity1 = min(stdDevStats)
    stdDevStats = np.array([])
    #for (int i = leftBound; i <= peakLeft; i++){
    for i in range(leftBound,peakLeft+1):
        stdDevStats = np.append(stdDevStats,intensities[i])

    smallIntensity2 = min(stdDevStats)
    smallIntensityAvg = (smallIntensity1+smallIntensity2)/2.0;
    
    if bestNoPeakNoise<1E-13:
        SNRatio = float("inf")
    else:
        SNRatio = (peakHeight-smallIntensityAvg)/bestNoPeakNoise
    #double SNRatio = (meanOfSignal-smallIntensityAvg)/bestNoPeakNoise;
    
    return SNRatio
