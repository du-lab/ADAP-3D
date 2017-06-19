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

import cwt.cwt as cwt
import cwt.correctbounds as correctbounds 
import cwt.signaltonoise as csn
import cwt.miscfilters as cmf
import pylab as pl

class PeakDetector:
    """
    The PeakDetector class is the interface for the user to the peak detection tools.
    The class essential takes a signal from the user and using a continuous wavelet transform method,
    detects peaks within that signal. The user can set parameters with set() like methods and get
    some additional information using get() like methods. After the peak detector has been configured
    the findPeaks() function will return the results in the form of a list of all the left boundaries
    of the peaks, a list of all the right boundaries of the peaks, and a list of all the positions
    of the peak maxima. The returned values are the indices of the boundaries and maxima respectively.
    
    An EXAMPLE can be found in the example_use_peakdetector.py script
    """
    def __init__(self):
        self.minimumIntensity = 1000.0
        self.coefOverAreaThreshold = 100.0
        self.SNThresh = 10.0
        self.waveletPeakSimilarity = .3

        self.signal = pl.array([])
        self.x = pl.array([])

        self.visualize  = False

        self.wavelet_small_scale = 1
        self.wavelet_large_scale = 10
        self.wavelet_scale_increment = 1

        self.peak_duration_range = [0.0,pl.inf]

        self.cur_peak_pos_list = []
        self.cur_left_bound_list = []
        self.cur_right_bound_list = []

        self.cur_coef_over_area_values = []
        self.cur_SN_values = []
        self.cur_wav_sim_vals = []

        self.verbose = False
            
        self.min_ridgeline_length = 5

    def getallSNVals(self):
        """ 
        Get a list of signal to noise ratios for all detected peaks.
        Call only after findPeaks() function.
        """
        return self.cur_SN_values

    def getallCoefOverAreaVals(self):
        """
        Get a list of coefficient over area values for all detected peaks.
        Call only after findPeaks() function.
        """
        return self.cur_coef_over_area_values

    def getallWaveletSimilarityVals(self):
        """
        Get a list of wavelet/peak similarity values for all detected peaks.
        Call only after findPeaks() function.
        """
        return self.cur_wav_sim_vals

    def setPeakDurationRange(self,rangeIn):
        """
        Set the peak duration range which is used to filter peaks by their widths.
        
        :param rangeIn: Two element list [lowest peak width, highest peak width].
        :type rangeIn: list [float, float]
        """
        self.peak_duration_range = rangeIn

    def setWaveletSmallScale(self,wssIn):
        """
        Set the smallest wavelet scale to be used in the CWT.
        
        :param wssIn: number->smallest scale.
        :type wssIn: float
        """
        self.wavelet_small_scale = wssIn

    def setWaveletLargeScale(self,wlsIn):
        """
        Set the largest wavelet scale to be used in the CWT.
        
        :param wssIn: number->largest scale.
        :type wlsIn: float
        """
        self.wavelet_large_scale = wlsIn

    def setWaveletPeakHighestSimilarity(self,hwpsIn):
        """
        The peak can be compared to the wavelet itself. The user can specify a similarity value
         below which all peaks will be discarded
         
        :param hwpsIn: Similarity threshold between 0 (completely different) and 1 (completely the same)
        :type hwpsIn: float
        """
        self.waveletPeakSimilarity = hwpsIn

    def setWaveletScaleIncrement(self,wciIn):
        """
        From the lowest scale to the highest -> what is the desired step size.
        
        :param wciIn: Scale step size.
        :type wciIn: float
        """
        self.wavelet_scale_increment = wciIn

    def setMinIntensityThreshold(self,minIntIn):
        """
        Set the intensity value below which all peaks will be discarded.
        
        :param minIntIn: Minimum intensity threshold value.
        :type minIntIn: float
        """
        self.minimumIntensity = minIntIn

    def setCoefOverAreaThreshold(self,minCoefOverA):
        """
        The coefficient over area is a way to filter out bad peaks based on their shape (some
        important subtleties not discussed here). The coefficient is the inner product
        of the wavelet with the peak at the best scale (largest inner product). The area is
        the area under the curve of the peak.
        
        :param minCoefOverA: Threshold value below which peaks will be discarded.
        :type minCoefOverA: float
        """
        self.coefOverAreaThreshold = minCoefOverA

    def setSignalToNoiseThreshold(self,minSN):
        """
        Signal to noise estimated from intensities. Choose a value near 10 and greater
        than 5 for best results. See cwt signal_to_noise for more information.
        
        :param minSN: minimum acceptable signal to noise ratio of peak.
        :type minSN: float
        """
        self.SNThresh = minSN

    def setSignal(self,signalIn):
        """
        Set the signal that peak detection will be preformed on.
        
        :param signalIn: Array of intensity values.
        :type signalIn: numpy array.
        """
        self.signal = signalIn

    def setX(self,xIn):
        """
        Set the values of the horizontal axis for the signal (retention time for XCMS data)
        
        :param xIn: Horizontal axis values (should be nearly evenly spaced).
        :type xIn: numpy array.
        """
        self.x = xIn

    def setVisualize(self,v):
        """
        Does the user want to visualize the peak picking procedure? Such visualization is not
        comprehensive but certainly useful and informative. Displays many aspects of procedure as
        it occurs.
        
        :param v: True or False 
        :type v: bool
        """
        self.visualize = v

    def setVerbose(self,v):
        """
        If true the peak picking procedure will print out some useful information about the peak picking
        process. This output does not accurately represent everything occurring in peak picking.
        
        :param v: True or False
        :type v: bool
        """
        self.verbose = v

    def setMinRidgelineLength(self,rlin):
        """
        This minimum length is very important and it can be dificult to choose a length threshold that
        minimizes false positives and maximizes false negatives. This length is the minimum number of
        connected local maxima in the wavelet coefficients required for a ridgeline to be considered indicative
        of an actual peak.
        
        :param rlin: Minimum ridgeline length
        :type rlin: int
        """
        self.min_ridgeline_length = rlin

    def findPeaks(self):
        """
        This is the function that detects the peaks in self.signal. It uses a continuous wavelet transform method
        to find their locations and boundaries as well as several filtering mechanism to reduce false positives.
        This function should only be called after the user has properly configured the PeakDetector object with 
        the desired options (S/N ratio, peak duration range, ... etc.)
        
        :return: List of all peaks left boundaries,
                list of all peaks right boundaries,
                list of all peak positions (position is position of maxima between boundaries).
        """
        
        finalLeftBoundList = []
        finalRightBoundList = []
        finalPeakPositionList = []
        self.cur_coef_over_area_values = []
        self.cur_SN_values = []

        if self.visualize:
            pl.plot(self.x,self.signal)
            pl.show()

        cwtObject = cwt.ContinuousWaveletTransform(
                        self.wavelet_small_scale,
                        self.wavelet_large_scale,
                        self.wavelet_scale_increment)
        # This line is probably a bad idea -> Should probably have this be a hardcoded value.
        #cwtObject.setMinRidgelineLength(self.peak_duration_range[0])
        cwtObject.setMinRidgelineLength(self.min_ridgeline_length)
        cwtObject.setVisualize(self.visualize)
        cwtObject.setSignal(self.signal)
        cwtObject.setX(self.x)
        cwtObject.buildRidgelines()
        cwtObject.filterRidgelines()

        # each column is one peak. 
        # row 0: left B
        # row 1: right B
        # row 2: best coefficient
        newPeaks = cwtObject.findBoundries();

        for i in range(len(newPeaks[0,:])):

            leftBound = int(newPeaks[0,i])
            rightBound = int(newPeaks[1,i])
            bestCoefficient = newPeaks[2,i]
            bestCoefficientLocation = int(newPeaks[3,i])
            bestScale = int(newPeaks[4,i])

            peakPosition = pl.argmax(self.signal[leftBound:rightBound])

            if self.visualize:
                pl.plot(self.x,self.signal)
                pl.fill_between(self.x[leftBound:rightBound+1],self.signal[leftBound:rightBound+1],facecolor='r')
                pl.title('peak before boundary correction')
                pl.show()

            # right bound has to go first before the leftBound value changes
            rightBound = correctbounds.fixRightBoundry(self.signal,
                                                     rightBound,
                                                     peakPosition+leftBound)
            leftBound = correctbounds.fixLeftBoundry(self.signal,
                                                     leftBound,
                                                     peakPosition +leftBound)
            leftBound,rightBound \
                = correctbounds.cropZerosFromEdges(self.signal,
                                                   leftBound,
                                                   rightBound)

            if self.visualize:
                pl.plot(self.x,self.signal)
                pl.fill_between(self.x[leftBound:rightBound+1],
                                self.signal[leftBound:rightBound+1],facecolor='r')
                pl.title('peak after boundary correction')
                pl.show()

            if (leftBound==0) and (rightBound == 0):
                if self.verbose:
                    print "PEAK NOT FOUND REASON: Left and right bound are = 0. " \
                          "It is porbably not a real peak becasue the boundry correction could not be done"
                continue

            rightBound+=1

            # get the peak position again because it needs to be with respect to the corrected
            # bounds
            peakPosition = pl.argmax(self.signal[leftBound:rightBound])

            # peak position should not also be the boundary
            if (leftBound==(peakPosition+leftBound))or(rightBound==(peakPosition+leftBound)):
                if self.verbose:
                    print "PEAK NOT FOUND REASON: Peak position is the same as boundry"
                    print "     leftBound: " + str(leftBound)
                    print "     rightBound: " + str(rightBound)
                    print "     peakPosition: " + str(peakPosition)
                continue

            # make sure the peak width is reasonable
            cur_peak_width = abs(rightBound - leftBound)
            if cur_peak_width<self.peak_duration_range[0]:
                if self.verbose:
                    print "PEAK NOT FOUND REASON: peakwidth too small."
                    print "     peak width: " + str(cur_peak_width)
                    print "     lower bound: " + str(self.peak_duration_range[0])
                    #print "     leftBound: " + str(leftBound)
                    #print "     rightBound: " + str(rightBound)
                    #print "     peakPosition: " + str(peakPosition)
                    #print "     self.peak_duration_range[0]: " + str(self.peak_duration_range[0])
                continue
            elif cur_peak_width>self.peak_duration_range[1]:
                if self.verbose:
                    print "PEAK NOT FOUND REASON: peakwidth too large."
                    print "     peak width: " + str(cur_peak_width)
                    print "     upper bound: " + str(self.peak_duration_range[1])
                    #print "     leftBound: " + str(leftBound)
                    #print "     rightBound: " + str(rightBound)
                    #print "     peakPosition: " + str(peakPosition)
                continue

            if self.waveletPeakSimilarity!="off":
                waveletPeakShapeBool, waveletPeakSimilarity = cmf.waveletPeakSimilarityCheck(cwtObject,leftBound,rightBound,bestScale,bestCoefficientLocation,self.waveletPeakSimilarity)
                if not waveletPeakShapeBool:
                    if self.verbose:
                        print "PEAK NOT FOUND REASON: wavelet similarity not met."
                        print "     Similarity = "+str(waveletPeakSimilarity)
                    continue
            else:
                waveletPeakSimilarity = -1


            if self.SNThresh!="off":
                curSN =csn.findBestEstimateOfSignalToNoiseRatio(self.signal,leftBound,rightBound)
                if curSN<self.SNThresh:
                    if self.verbose:
                        print "PEAK NOT FOUND REASON: S/N not met."
                    continue
            else:
                curSN = -1.0

            if not cmf.minimumIntensityCheck(self.signal,leftBound,rightBound,self.minimumIntensity):
                if self.verbose:
                    print "PEAK NOT FOUND REASON: Minimum intensity not met."
                continue

            if self.coefOverAreaThreshold!="off":
                coefOverAreaBool, curCoefOverArea = cmf.coefficientOverAreaCheck(self.x,self.signal,leftBound,rightBound,bestCoefficient,self.coefOverAreaThreshold)
                if not coefOverAreaBool:
                    if self.verbose:
                        print "PEAK NOT FOUND REASON: C/A not met."
                        print "     C/A = "+str(curCoefOverArea)
                    continue
            else:
                curCoefOverArea = -1

            if not cmf.numberOfZeros(self.signal,leftBound,rightBound):
                if self.verbose:
                    print "PEAK NOT FOUND REASON: Number of zero points greater than number of non-zero points."
                #print "h7"
                continue
            if not cmf.numNonZeroGreaterThanThree(self.signal,leftBound,rightBound):
                if self.verbose:
                    print "PEAK NOT FOUND REASON: Only has 3 non-zero points -> should not consider this a real peak."
                continue

            finalLeftBoundList.append(leftBound)
            finalRightBoundList.append(rightBound)

            finalPeakPositionList.append(peakPosition+leftBound)

            self.cur_coef_over_area_values.append(curCoefOverArea)
            self.cur_SN_values.append(curSN)
            self.cur_wav_sim_vals.append(waveletPeakSimilarity)

            self.cur_peak_pos_list = finalPeakPositionList 
            self.cur_left_bound_list = finalLeftBoundList
            self.cur_right_bound_list = finalRightBoundList

        return finalLeftBoundList,finalRightBoundList,finalPeakPositionList 

    def showSNAnnotatedEIC(self):
        """
        Plot the self.signal with all the peaks detected shaded and each annotated with
        the estimated signal to noise ratio.
        """
        pl.plot(self.x,self.signal)
        mostMaxInten = 0.0
        for i in range(len(self.cur_peak_pos_list)):
            pl.fill_between(self.x[self.cur_left_bound_list[i]:self.cur_right_bound_list[i]+1],self.signal[self.cur_left_bound_list[i]:self.cur_right_bound_list[i]+1],facecolor='r')
            peakIndex = self.cur_peak_pos_list[i]
            curSN = self.cur_SN_values[i]
            if curSN > 1E12:
                strSN = "inf"
            else:
                strSN = "%0.1f"%(curSN)
            pl.annotate(strSN,xy=(self.x[peakIndex],self.signal[peakIndex]),
                    xytext = (self.x[peakIndex],self.signal[peakIndex]+self.signal[peakIndex]*0.2),
                    arrowprops=dict(facecolor='black',shrink = 0.05))
            if self.signal[peakIndex]>mostMaxInten:
                mostMaxInten = self.signal[peakIndex]

        pl.ylim([0,mostMaxInten+mostMaxInten*0.4])
        pl.show()

