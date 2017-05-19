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

import peakdetector as pd
import pylab as pl

def main():

    #f = open('TestEICs/testing_eic.txt','r')
    f = open('TestEICs/testing_eic2.txt','r')
    d = pl.genfromtxt(f)
    intensities = d[0,:]
    rt = d[1,:]
    
    #print 'intensities: ' + str(type(intensities))
    #exit()
    
    peakdetector = pd.PeakDetector()

    # intensities and rt (inputs) must be numpy arrays
    peakdetector.setSignal(intensities)
    peakdetector.setX(rt)
    # These are just the default values. Just showing how to set
    peakdetector.setMinIntensityThreshold(100)
    peakdetector.setCoefOverAreaThreshold(10)
    peakdetector.setWaveletPeakHighestSimilarity("off")
    # This can be a number or it can be "off" to tell the program not to calculate the SN threshold.
    peakdetector.setSignalToNoiseThreshold(5)
    peakdetector.setVisualize(False)
    allLeftBounds,allRightBounds , allPeakPositions = peakdetector.findPeaks()

    allSNValues = peakdetector.getallSNVals()

    peakdetector.showSNAnnotatedEIC()

    for i in range(len(allRightBounds)):
        print "--------- peak #"+str(i)+"------------"
        print "Right bound: " + str(allRightBounds[i]) + " Left bound: " + str(allLeftBounds[i])
        print "Peak position: " + str(allPeakPositions[i])
        print "Peak SN: "+str(allSNValues[i])
        print "------------------------------"

if __name__=="__main__":
    main()
