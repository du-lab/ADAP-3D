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

PERCENT_PEAK_HEIGHT_STOP = 0.05

def fixLeftBoundry(intensities, peakLeft, peakPosition):
    """
    Correct the boundaries that were determined through the CWT initial boundary estimation method. 
    The correction really just descends along the slope of the curve until it hits a local minima, 
    OR until the intensity is less than PERCENT_PEAK_HEIGHT_STOP*<the peaks max intensity>. 
    :param intensities: The y-axis values of the input signal.
    :type intensities: numpy array
    :param peakLeft: The current index of the left boundary of the peak.
    :type peakLeft: int
    :param peakPosition: The index of the peaks maximum intensity point.
    :type peakPosition: int
    :return: index of the corrected left boundary.
    """

    foundLocalMin = False
    curLeft = peakLeft
    bestLeft = peakLeft
    
    peakIntensity = intensities[peakPosition]

    while (not foundLocalMin):
        checkLeftPoint = curLeft-1
       
        checkRightPoint = curLeft+1

        if (checkRightPoint>=len(intensities)):
            return peakLeft
        
        if (checkLeftPoint < 0):
            bestLeft = curLeft
            foundLocalMin = True

        #elif (intensities[curLeft]<1.0):
        #    bestLeft = curLeft
        #    foundLocalMin = True

        elif ((intensities[checkRightPoint]>=intensities[curLeft])
                and (intensities[checkLeftPoint]>=intensities[curLeft])):
            bestLeft = curLeft
            foundLocalMin = True

        elif(intensities[checkLeftPoint] <= intensities[curLeft]):
            curLeft -=1

        elif (intensities[checkRightPoint]<intensities[curLeft]):
            curLeft +=1

        else:
            return -1

        # the condition regarding the intesities is a little add hock but should prevent boundries from
        # being too wide for very intense peaks
        if (intensities[curLeft]<(PERCENT_PEAK_HEIGHT_STOP*peakIntensity)):
            bestLeft = curLeft
            foundLocalMin = True


    return bestLeft;


def fixRightBoundry(intensities, peakRight, peakPosition):
    """
    Correct the boundaries that were determined through the CWT initial boundary estimation method. 
    The correction really just descends along the slope of the curve until it hits a local minima, 
    OR until the intensity is less than PERCENT_PEAK_HEIGHT_STOP*<the peaks max intensity>. 
    :param intensities: The y-axis values of the input signal.
    :type intensities: numpy array
    :param peakRight: The current index of the right boundary of the peak.
    :type peakRight: int
    :param peakPosition: The index of the peaks maximum intensity point.
    :type peakPosition: int
    :return: index of the corrected right boundary.
    """

    foundLocalMin = False
    curRight = peakRight
    bestRight = peakRight
    
    peakIntensity = intensities[peakPosition]
    
    while (not foundLocalMin):
        checkRightPoint = curRight+1
       
        checkLeftPoint = curRight-1
        if (checkLeftPoint<=0):
            return peakRight
        
        if (checkRightPoint >= len(intensities)):
            bestRight = curRight
            foundLocalMin = True
    
        #elif (intensities[curRight]<1.0):
        #    bestRight = curRight
        #    foundLocalMin = True

        elif ((intensities[checkRightPoint]>=intensities[curRight]) and
                (intensities[checkLeftPoint]>=intensities[curRight])):
            bestRight = curRight
            foundLocalMin = True
    
        elif(intensities[checkRightPoint] <= intensities[curRight]):
            curRight += 1
    
        elif (intensities[checkLeftPoint]<intensities[curRight]):
            curRight -= 1

        else:
            return -1

        # the condition regarding the intesities is a little add hock but should prevent boundries from
        # being too wide for very intense peaks
        if (intensities[curRight]<(PERCENT_PEAK_HEIGHT_STOP*peakIntensity)):
            bestRight = curRight
            foundLocalMin = True

        
    return bestRight



def cropZerosFromEdges(intensity,peakLeft,peakRight):
    """
    It is possible for the bounds to contain null or 0 intensity points. For example if 
    for some reason the CWT finds a peak with only a single point in it but it thinks the 
    bounds are a couple of points to the left and right then it would pass the above if statement. 
    To make sure we get rid of these points we need to do one more check which is below.
    
    :param intensity: The y-axis values of the input signal.
    :type intensity: numpy array
    :param peakLeft: The current index of the left boundary of the peak.
    :type peakLeft: int
    :param peakRight: The current index of the right boundary of the peak.
    :type peakRight: int
    :return: Returns the new boundaries after crop ([new left,new right]). 
    if there is something wrong with the peak it will return boundaries of [0,0].
    """
    croppedPeakLeft=-1
    croppedPeakRight=-1
    
    allZero=True
    #for (int alpha=peakLeft; alpha<peakRight; alpha++):
    for alpha in range(peakLeft,peakRight):
        curInt = intensity[alpha]
        if (curInt!=0.0):
            
            allZero = False
            break

    if (allZero):
        return 0,0

    #for (int alpha=peakLeft; alpha<peakRight; alpha++){
    for alpha in range(peakLeft,peakRight):
        curInt = intensity[alpha]
        if (curInt!=0.0):
            croppedPeakLeft = alpha
            break

    #for (int alpha=peakRight; alpha>peakLeft; alpha--){
    for alpha in range(peakRight,peakLeft,-1):
        curInt = intensity[alpha]
        if (curInt!=0.0):
            croppedPeakRight = alpha
            break

    # Everything could be zero
    #if (croppedPeakRight==-1){
    #    croppedPeakRight = peakRight;
    #}
    #if (croppedPeakLeft==-1){
    #    croppedPeakLeft = peakLeft;
    #}
    
    # the most left and right points could/should be zero so by adding/subtracting from alpha we can make sure that  remains the case
    # Otherwise the peak width is not acuretly being represented.
    if (croppedPeakLeft!=peakLeft):
        croppedPeakLeft-=1

    if (croppedPeakRight!=peakRight):
        croppedPeakRight+=1

    if(croppedPeakLeft==-2):
        print("bug in correctbounds.py")
        exit()

    return croppedPeakLeft,croppedPeakRight
    
    #croppedRetentionTimeRight = retentionTimes[croppedPeakRight];
    #double croppedRetentionTimeLeft = retentionTimes[croppedPeakLeft];
    #if(! peakWidth.contains(croppedRetentionTimeRight- croppedRetentionTimeLeft))
    #{
    #    continue;
    #}
