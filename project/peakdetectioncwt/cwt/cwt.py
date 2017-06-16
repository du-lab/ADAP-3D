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

import pylab as pl
import ridgeline 

class ContinuousWaveletTransform:
    """
    This class contains the parameters and implementation of the continuous wavelet transform. It only
    requires the Ridgeline objects and is limited to the most basic elements of a CWT. We chose to 
    separate estimates of signal to noise ratio and other filtering mechanisms from this implementation
    to keep it as simple as possible.
    """
    def __init__(self,smallScaleIn,largeScaleIn,incrementScaleIn):
        # visualize things as you go?
        self.visualize = False

        self.smallScale= smallScaleIn
        self.largeScale= largeScaleIn
        self.incrementScale= incrementScaleIn

        self.signal = []
        self.x = []
        self.avgXSpace = []
        self.allCoefficients = pl.array([])
        self.minRidgelineLength =6

        # how far in each direction from the current point do we need to grab data for a succesful wavelet transform?
        # This number is the factor we multiply by the scale. 5 should be good because this is the estimated compact support
        self.scaleCoefHowFarOut = 5;

        self.arrScales = pl.array([])
        self.mapScaleToIndex = {}

        self.ridgeLineArr = []

        index = 0

        self.arrScales = pl.arange(self.smallScale,self.largeScale+self.incrementScale,self.incrementScale)
        for curScale in self.arrScales:
            self.mapScaleToIndex[curScale] = index
            index += 1

    def setMinRidgelineLength(self,minIn):
        """
        This minimum length is very important and it can be dificult to choose a length threshold that
        minimizes false positives and maximizes false negatives. This length is the minimum number of
        connected local maxima in the wavelet coefficients required for a ridgeline to be considered indicative
        of an actual peak.
        
        :param minIn: Minimum ridgeline length 
        :return: int
        """
        self.minRidgelineLength = minIn

    def setVisualize(self,visIn):
        """
        Does the user want to visualize the peak picking procedure? Such visualization is not
        comprehensive but certainly useful and informative. Displays many aspects of procedure as
        it occurs.
        
        :param visIn: True or False 
        :type visIn: bool
        """
        self.visualize = visIn

    def findBoundries(self):
        """
        Using the best scale (the smallest scale at which there is a local maxima in the coefficients of 
        the ridgeline) of each peak estimate the boundaries of the peaks.
        :return: 2D list of 
                list of left boundaries, 
                list of right boundaries, 
                list of best coefficient of each peak,
                list of the index (location in signal) of best coefficient,
                list of the best scales.
        """
        boundsAndBestCoef = pl.zeros([5,len(self.ridgeLineArr)])

        count = 0
        for curRL in self.ridgeLineArr:
            bestIndex = curRL.getBestIndex()
            # this is the actuale scale, not the index of the best scale.
            bestScale = curRL.getBestScale()

            bestCoefficient = curRL.getMaxCor();
            
            curRightBound = bestIndex+int(bestScale+0.5);
            if (curRightBound>=len(self.x)):
                curRightBound = len(self.x)-1;
             
            curLeftBound = bestIndex-int(bestScale+0.5);
            if (curLeftBound<0):
                curLeftBound = 0;
            
            boundsAndBestCoef[0,count] = curLeftBound;
            boundsAndBestCoef[1,count] = curRightBound;
            boundsAndBestCoef[2,count] = bestCoefficient;
            boundsAndBestCoef[3,count] = bestIndex;
            boundsAndBestCoef[4,count] = bestScale;
            
            count+=1;

        return boundsAndBestCoef

    def filterRidgelines(self):
        """
        Get rid of ridgelines that are too short i.e. do not mean the minimum ridgeline length parameter. 
        --> Modifies self.ridgeLineArr by removing bad RLs.
        """

        filteredRidgelines = []
        
        if (self.visualize):
            self.plotRidgelinesOverCoefs('before length filter')
        
        for i in range(len(self.ridgeLineArr)):
            #Ridgeline curRL = ridgeLineArr[i];
            curRL = self.ridgeLineArr[i]
            ridgeLength = curRL.getRidgeLength()
            NScales = curRL.totalNumberOfScales
            
            # When we make this CWT more general this check should be in terms of some precentage of the total number of scales.
            #if (ridgeLength<(NScales-4)):
            #    continue
            if (ridgeLength<self.minRidgelineLength):
                continue

            filteredRidgelines.append(curRL)
            
        self.ridgeLineArr = filteredRidgelines;

        if (self.visualize):
            self.plotRidgelinesOverCoefs('after length filter')

    def buildRidgelines(self):
        """
        By looking for local maxima in the coefficients at different scales but are close in the
        retention time domain the ridgelines are built. Ridgeline objects contain the sets of 
        connected local maxima and perform some simple parts of the building algorithm as well.
        Adds ridgelines to self.ridgeLineArr as they are built.
        """

        self.getCoefficientsForAllScales()

        # start from the largest scale and go to the smallest
        #for (int i=arrScales.size()-1; i>=0; i--){
        for i in range(len(self.arrScales)-1,0,-1):
            curScale = self.arrScales[i]
            indexOfThisWaveletScale = self.mapScaleToIndex[curScale]
            curCoefficients = self.allCoefficients[indexOfThisWaveletScale]

            #if self.visualize:
            #    print "curScale: " + str(curScale)
            #    print "indexOfThisWaveletScale: " + str(indexOfThisWaveletScale)
            
            
            thisScaleBestMaxima = self.findMaximaForThisScale(curScale);

            #if self.visualize:
            #    print "thisScaleBestMaxima: " + str(thisScaleBestMaxima)
            
            
            #for (int j=0; j<thisScaleBestMaxima.length;j++){
            for j in range(0,len(thisScaleBestMaxima)):
                wasMatched = False;
                curBestMaxLoc = thisScaleBestMaxima[j];


                #for (int alpha=0; alpha<ridgeLineArr.size(); alpha++){
                for alpha in range(0,len(self.ridgeLineArr)):
                    
                    wasAdded = self.ridgeLineArr[alpha].tryAddPoint(
                                        curScale,
                                        curBestMaxLoc,
                                        curCoefficients[thisScaleBestMaxima[j]])
                    
                    if (wasAdded):
                        wasMatched=True
                    

                #if self.visualize:
                #    print "wasMatched: " + str(wasMatched)
                #    print "curBestMaxLoc: " + str(curBestMaxLoc)

                # if it was not added to at least one then make a new redge line
                if (not wasMatched):
                    
                    curStartRidge = ridgeline.Ridgeline(curScale,
                                              thisScaleBestMaxima[j],
                                              curCoefficients[thisScaleBestMaxima[j]],
                                              len(self.arrScales),
                                              self.visualize)
                                              
                    self.ridgeLineArr.append(curStartRidge)

                #if self.visualize:
                #    self.plotRidgelinesOverCoefs('watch build')

        #writeRidgelines()


    # returns the indecies of the location of the maxima
    def findMaximaForThisScale(self,waveletScale):
        """
        Starting with the most intense local maxima and then removing points near ("near" based on
        particular wavelet scale under investigation) that maxima the first is detected. After
        the removal the next most intense local maxima should be the next real peak in the 
        coefficients and the process is repeated until all points have been removed.
        
        :param waveletScale: The scale for which the maxima are desired.
        :return: A list of indices corresponding to the locations of the maxima.
        """

        # when we are removing points adjacent to the current maxima this is the 
        # number of points to go in either direction before stopping.
        #removeCutOff = int(round(waveletScale*2.5))
        #removeCutOff = int(round(waveletScale*2.5))
        removeCutOff = int(round(waveletScale*2.5))
        
        maximaLocations = []
        
        # sort and keep track of the original idecies
        indexOfThisWaveletScale = self.mapScaleToIndex[waveletScale]

        curCoefficients = pl.copy(self.allCoefficients[indexOfThisWaveletScale])

        #if self.visualize:
        #    print "waveletScale: " + str(waveletScale)
        #    pl.plot(curCoefficients,marker='o')
        #    pl.show()

        # sorted stuff will be from smallest to greates
        indecies = sorted(range(len(curCoefficients )), key=lambda k: curCoefficients[k])
        #curCoefficients.sort()

       
        mapIndexToBoolRemain = {}
        
        for i in range(0,len(indecies)):
            mapIndexToBoolRemain[indecies[i]]=True


        #for (int i = indecies.length-1; i>=0; i--){
        for i in range(len(indecies)-1,0-1,-1):
            if (mapIndexToBoolRemain[indecies[i]]):

                curLargestIndex = indecies[i]

                #if self.visualize:
                #    print "curLargestIndex: " + str(curLargestIndex)

                if (curCoefficients[curLargestIndex] <= 0.0):
                    continue

                maximaLocations.append(curLargestIndex)
                # remove points. num points to right and left equal to current scale
                mapIndexToBoolRemain[curLargestIndex]=False

                #for (int j=1; j<removeCutOff; j++){
                for j in range(1,removeCutOff):
                  
                    curRemoveIndexRight = curLargestIndex+j
                    curRemoveIndexLeft = curLargestIndex-j

                    #if self.visualize and (j==(removeCutOff-1)):
                    #    print "curRemoveIndexRight: " + str(curRemoveIndexRight)
                    #    print "curRemoveIndexLeft: " + str(curRemoveIndexLeft)

                    if(curRemoveIndexLeft>=0):
                        mapIndexToBoolRemain[curRemoveIndexLeft]=False

                    if(curRemoveIndexRight<len(self.x)):
                        mapIndexToBoolRemain[curRemoveIndexRight]=False


        #if self.visualize and (j==(removeCutOff-1)):
        #    print "maximaLocations: " + str(maximaLocations)
        
        return maximaLocations
        

    def getCoefficientsForAllScales(self):
        """
        Fill the matrix of self.allCoefficents by performing the CWT at each scale on the signal.
        """

        NScales = len(self.arrScales)
        self.allCoefficients = pl.zeros([NScales,len(self.x)])
        count = 0
        #for (Integer curScale: arrScales){
        for curScale in self.arrScales:
            self.allCoefficients[count] = self.getCoefficientsForThisScale(curScale);
            count += 1;
        #writeAllCoeffs();

        if (self.visualize):
            pl.pcolor(self.allCoefficients)
            pl.title('coefficients')
            pl.show()

    def getCoefficientsForThisScale(self,waveletScale):
        """
        For the given scale take the inner product of the wavelet with the signal for all
        translations of the wavelet.
        
        :param waveletScale: Input wavelet scale.
        :return: Return numpy array of all the coefficients found for each translation of the wavelet.
        """

        coefficientsForThisScale = pl.zeros([len(self.x)])
        #for (int i=0; i<x.length; i++){
        for i in range(len(self.x)):
            currentCoefficient = self.signalWaveletInnerProductOnePoint(i,waveletScale);
            coefficientsForThisScale[i] = currentCoefficient;
            
        return coefficientsForThisScale


    

    def signalWaveletInnerProductOnePoint(self,xIndexOfWaveletMax, waveletScale):
        """
        This function needs to be careful with two things: 
        1) the boundaries of the signal were it needs to either pad or pretend to 
        pad values below off the boundary. 
        2) The location of the wavelet has to be set correctly.
        
        Note: waveletScale is in units of indices NOT rt or anything else
        
        :param xIndexOfWaveletMax: Index of wavelet maxima
        :type xIndexOfWaveletMax: int
        :param waveletScale: Current scale of wavelet.
        :type waveletScale: float
        :return: (float) the value of the inner product between the wavelet and the signal
            for this particular translation and scale of the wavelet.
        """

        leftBoundIntegrate = xIndexOfWaveletMax-self.scaleCoefHowFarOut*int(round(waveletScale))
        rightBoundIntegrate = xIndexOfWaveletMax+self.scaleCoefHowFarOut*int(round(waveletScale))
        if (leftBoundIntegrate<0):
            leftBoundIntegrate=0
        if (rightBoundIntegrate>=len(self.x)):
            rightBoundIntegrate=len(self.x)-1
        
        curX = pl.zeros([rightBoundIntegrate-leftBoundIntegrate +1])
        curY = pl.zeros([rightBoundIntegrate-leftBoundIntegrate +1])
        waveletY = pl.zeros([rightBoundIntegrate-leftBoundIntegrate +1])
        
        curIndex = 0;
        #for (int i=leftBoundIntegrate;i<=rightBoundIntegrate; i++){
        for i in range(leftBoundIntegrate,rightBoundIntegrate+1):
            curX[curIndex] = self.x[i]
            curY[curIndex] = self.signal[i]
            # for the wavelt work in units of indecies because wavelt for  numbers smaller than one is not approriate.
            waveletY[curIndex] = self.rickerWavelet(self.x[i]-self.x[xIndexOfWaveletMax], waveletScale)
            curIndex+=1

        #double[] doublePtsCurX = doubleTheNumberOfPtsX(curX);
        #double[] doublePtsDataY = doubleTheNumberOfPtsDataY(curY);
        #double[] doublePtsWavelet = doubleTheNumberOfPtsWavelet(waveletY, (double) waveletScale);
        
        #writeDataAndWavelet(curX,curY,waveletY);
        
        innerProd = self.innerProduct(curX,curY,waveletY)
        return innerProd


    # This just takes an x value and the parameters of the wavelet and retuns the y value for that x
    def rickerWavelet(self, x, scalParam):
        """
        The mathmatical form of the ricker (Mexican hat) wavelet.
        :param x: horizontal axis input
        :type x: float
        :param scalParam: wavelet scale
        :type scalParam: float
        :return: Returns the function evaluated at x for the given scale.
        """

        scalParam = scalParam*self.avgXSpace
        A = 2.0/pl.sqrt(3.0 * scalParam*pl.sqrt(pl.pi)) * (1.0-(x**2.0)/(scalParam**2.0))
        return pl.exp(-(x**2.0)/(2.0*(scalParam**2)))*A

    def innerProduct(self,x, arr1, arr2):
        """
        This function can only take two arrays of equivalent length.
        in the msconvert code they just add the wavelet * the intensity... Lets just do this 
        for now to see if we can get the same results.
        
        :param x: horizontal axis values
        :param arr1: array of first curves y axis values
        :param arr2: array of second curves y axis values
        :return: The value of the inner product
        """

        multArr = arr1*arr2
        return multArr.sum()
        ## This is the way it should be done but becase msconvert only did the above we are ommiting
        ## this for now to keep the inner product results consistent with those of msconvert.
       # // Because EICs can be messy best to just use trapezoidal rule
       # double area = 0.0;
       # for (int i=0; i < l-1; i++){
       #     double curXSpace = x[i+1]-x[i];
       #     // lowest height of the two adjacent points
       #     double curYLow = Math.min(multArr[i+1],multArr[i]);
       #     double curYHigh= Math.max(multArr[i+1],multArr[i]);
       #     double triangleArea = 0.5*curXSpace*(curYHigh-curYLow);
       #     // triangle area needs to be set as negative if the points are below zero
       #     if (){
       #         
       #     }
       #     double rectangleArea = curXSpace*curYLow;
       #     
       # }

    def setSignal(self, signalIn):
        """
        Set the signal that the CWT will be performed on.
        
        :param signalIn: intensities (y-axis values).
        :type signalIn: numpy array
        """
        self.signal = signalIn
        if (self.visualize):
            pl.plot(self.signal)
            pl.title('signal input')
            pl.show()

    def setX(self, xIn):
        """
        Set the horizontal axis values of the curve that the CWT will be perfomed on.
        
        :param xIn: input horizontal values.
        :type xIn: numpy array.
        """

        self.x = xIn

        spaceArr = xIn[1:]-xIn[:-1]

        self.avgXSpace = spaceArr.mean() 

        if (self.visualize):
            pl.plot(self.x,self.signal)
            pl.title('signal with rt')
            pl.show()

    def plotRidgelinesOverCoefs(self,titleString):
        """ 
        Useful method to visualize the coefficients matrix as well as the
        ridgelines that were built from it.
        
        :param titleString: Title of the plot
        :type titleString: String
        """

        fig = pl.figure()
        ax = fig.add_subplot(211)
        ax2 = fig.add_subplot(212)

        
        rls = pl.zeros(pl.shape(self.allCoefficients))

        ax.pcolor(self.allCoefficients)
        
        for i in self.ridgeLineArr:
            for j in range(len(i.scales_)):
                curScaleIndex = self.mapScaleToIndex[i.scales_[j]]
                curIndex = i.indecies_[j]
                rls[curScaleIndex,curIndex] = 100

        ax2.pcolor(rls)
        ax.set_title('ridgelines '+titleString)
        pl.show()

