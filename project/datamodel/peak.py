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

class Peak:
    def __init__(self):
        self.rt = -1.0
        self.mz = -1.0
        self.height = -1.0
        self.mzMin = -1.0
        self.mzMax = -1.0
        self.rtMin = -1.0
        self.rtMax = -1.0
        self.mzIndex = -1
        self.rtIndex = -1
        self.mzMinIndex = -1
        self.mzMaxIndex = -1
        self.rtMinIndex = -1
        self.rtMaxIndex = -1
        self.hasIsotope = False
        self.isotopeList = []
        self.coefOverArea = 0.0
        self.asymGausFit  = 0.0
        self.similarityAdjMZ = []
    def setSimilarityAdjMZ(self,tobe):
        self.similarityAdjMZ = tobe
    def setAsymGausFit(self,tobe):
        self.asymGausFit = tobe
    def setCoefOverArea(self,tobe):
        self.coefOverArea = tobe
    def setRTMinIndex(self,tobe):
        self.rtMinIndex = tobe
    def setRTMaxIndex(self,tobe):
        self.rtMaxIndex = tobe
    def setMZMinIndex(self,tobe):
        self.mzMinIndex = tobe
    def setMZMaxIndex(self,tobe):
        self.mzMaxIndex = tobe
    def setRTIndex(self,tobe):
        self.rtIndex = tobe
    def setMZIndex(self,tobe):
        self.mzIndex = tobe
    def setRT(self,tobe):
        self.rt = tobe
    def setMZ(self,tobe):
        self.mz = tobe
    def setHeight(self,tobe):
        self.height = tobe
    def setMZmin(self,tobe):
        self.mzMin = tobe
    def setMZmax(self,tobe):
        self.mzMax = tobe
    def setRTMin(self,tobe):
        self.rtMin = tobe
    def setRTMax(self,tobe):
        self.rtMax = tobe
    def addIsotope(self,isotopePeakObject):
        self.isotopeList.append(isotopePeakObject)
        self.hasIsotope = True
