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

NEGATIVE_INFINITY = float("-inf")

class Ridgeline:
    def __init__(self,
                firstScale,
                firstIndex,
                corValue,
                NScales,
                visualize):

        self.scales_ = []
        self.indecies_ = []
        self.corValues_ = []

        self.scales_.append(firstScale)
        self.indecies_.append(firstIndex)
        self.corValues_.append(corValue)
        self.curRunningGap_ = 0
        self.totalNumberOfScales = NScales

        self.visualize=visualize

    def getRunningGapNum(self):
        return self.curRunningGap_
    

    def getRidgeLength(self):
        return len(self.scales_)

    def getRidgeStartScale(self):
        return self.scales_[0]

    def getRidgeEndScale(self):
        return self.scales_[-1]

    def getBestIndex(self):
        # find local maxima
        locMaxCorList =[]
        locMaxIndexList =[]
        for i in range(1,len(self.indecies_)-1):
            curCor = self.corValues_[i]
            curIndex = self.indecies_[i]
            if (curCor>self.corValues_[i-1])and(curCor>self.corValues_[i+1]):
                locMaxCorList.append(curCor)
                locMaxIndexList.append(curIndex)

        #check very edges
        if self.corValues_[0]>self.corValues_[1]:
            locMaxCorList  =[self.corValues_[0]]+locMaxCorList
            locMaxIndexList=[self.indecies_[0]]+locMaxIndexList
        if self.corValues_[-1]>self.corValues_[-2]:
            locMaxCorList  . append(self.corValues_[-1])
            locMaxIndexList. append(self.indecies_[-1])

        curBestIndex = locMaxIndexList[-1]


        # original way
        #curBestInd=-1
        #maxCorVal=NEGATIVE_INFINITY
        ##for (int i=0; i< self.indecies_.size();i++){
        #for i in range(0,len(self.indecies_)):
        #    curCor = self.corValues_[i]
        #    if (curCor>maxCorVal):
        #        maxCorVal = curCor
        #        curBestInd = self.indecies_[i]

        return curBestIndex;

    # we want this to find how many local maxima and then choose the local max coresponding to the
    # smallest scale to avoid lumping close peaks together.
    def getBestScale(self):

        # find local maxima
        locMaxCorList =[]
        locMaxScaleList =[]
        #for (int i=0; i< self.scales_.size();i++){
        for i in range(1,len(self.scales_)-1):
            curCor = self.corValues_[i]
            curScale = self.scales_[i]
            if (curCor>self.corValues_[i-1])and(curCor>self.corValues_[i+1]):
                locMaxCorList.append(curCor)
                locMaxScaleList.append(curScale)

        #check very edges
        if self.corValues_[0]>self.corValues_[1]:
            locMaxCorList  =[self.corValues_[0]]+locMaxCorList
            locMaxScaleList=[self.scales_[0]]+locMaxScaleList
        if self.corValues_[-1]>self.corValues_[-2]:
            locMaxCorList  . append(self.corValues_[-1])
            locMaxScaleList. append(self.scales_[-1])

        curBestScale = locMaxScaleList[-1]

        if self.visualize:
            print "curBestScale: " + str(curBestScale) 
            print "locMaxScaleList: " + str(locMaxScaleList) 
            pl.plot(self.scales_,self.corValues_)
            pl.scatter(locMaxScaleList,locMaxCorList,c='r')
            pl.title("all local max RL coefs")
            pl.show()
            pl.plot(self.scales_,self.corValues_)
            pl.scatter(locMaxScaleList[-1],locMaxCorList[-1],c='r')
            pl.title("smallest scale local max RL coefs")
            pl.show()



        #Original way
        #curBestScale=-1
        #maxCorVal=NEGATIVE_INFINITY
        ##for (int i=0; i< self.scales_.size();i++){
        #for i in range(0,len(self.scales_)):
        #    curCor = self.corValues_[i]
        #    if (curCor>maxCorVal):
        #        maxCorVal = curCor
        #        curBestScale = self.scales_[i]

        return curBestScale

    def getMaxCor(self):
        # find local maxima
        locMaxCorList =[]
        locMaxScaleList =[]
        #for (int i=0; i< self.scales_.size();i++){
        for i in range(1,len(self.scales_)-1):
            curCor = self.corValues_[i]
            curScale = self.scales_[i]
            if (curCor>self.corValues_[i-1])and(curCor>self.corValues_[i+1]):
                locMaxCorList.append(curCor)
                locMaxScaleList.append(curScale)

        #check very edges
        if self.corValues_[0]>self.corValues_[1]:
            locMaxCorList  =[self.corValues_[0]]+locMaxCorList
            locMaxScaleList=[self.scales_[0]]+locMaxScaleList
        if self.corValues_[-1]>self.corValues_[-2]:
            locMaxCorList  . append(self.corValues_[-1])
            locMaxScaleList. append(self.scales_[-1])

        maxCorVal = locMaxCorList[-1]


        # original way
        #maxCorVal=0.0
        ##for (int i=0; i< self.corValues_.size();i++){
        #for i in range(0,len(self.corValues_)):
        #    curCor = self.corValues_[i]
        #    if (curCor>maxCorVal):
        #        maxCorVal = curCor

        return maxCorVal

    def tryAddPoint(self,scale, index, corValue):

        # see if where this index is in relation to the last added
        lastAddedInd = self.indecies_[-1]
        indexDiff = abs(lastAddedInd - index)
        indexTol = self.findIndexTolFromScale(scale)

        #print "lastAddedInd: " + str(lastAddedInd )
        #print "indexDiff: " + str(indexDiff)
        #print "indexTol: " + str(indexTol)

        # Need to see if something has already been added for this scale
        haveThisScaleAlready = False

        if (self.scales_[-1]==scale):
            haveThisScaleAlready = True;

        #print "haveThisScaleAlready: " + str(haveThisScaleAlready)
        if (not haveThisScaleAlready ):

            # times 2 for pluss minus tollerance
            if (indexDiff<(2*indexTol)):

                self.scales_.append(scale)
                self.indecies_.append(index)
                self.corValues_.append(corValue)
                self.curRunningGap_ = 0
                return True;
            else:
                self.curRunningGap_+=1
                return False;

        else:
            # two things to check 
            # 1) is it closer in endex to previous?
            # 2) is it larger or smaller correlation value
            # For now lets just take the closest point unless this the first scale still.
            # If it is the first scale then lets pick the largest value.
            if (len(self.scales_) > 1):
                # Lets try taking the largest one instead
                # prevCor = self.corValues_[-1]
                # curCor = corValue

                #if (curCor>prevCor):
                #    self.indecies_[-1]=index
                #    self.corValues_[-1]=corValue
                #    return True

                prevIndexDiff = abs(self.indecies_[-2]-self.indecies_[-1])
                # this is -2 comparaed to index because the other one to compare has already been
                # added want to compare index to the -2th one also
                curIndexDiff = abs(self.indecies_[-2]-index)

                if (prevIndexDiff>curIndexDiff ):
                    self.indecies_[-1]=index
                    self.corValues_[-1]=corValue
                    return True

            else :
                # only compare magnitued if they are close points
                if (indexDiff<(2*indexTol)):
                    prevCor = self.corValues_[0]

                    if (corValue>prevCor):
                        self.indecies_[-1]=index
                        self.corValues_[-1]=corValue
                        return True

        return False
    
    def findIndexTolFromScale(self,scale):
        # This is probably too simple but it apears to work well enough. Use for now ans then look into
        # correct value

        # window size for scale = 1 is [-5,5]
        #return scale*5;
        #
        #I think above is much to big. Going to try this
        return int(scale+0.5);
