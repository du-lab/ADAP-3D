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

class Result:
    def __init__(self):
        self.final_peak_list_rt_vals = []
        self.final_peak_list_mz_vals = []
        self.final_peak_list_peak_height_vals = []
        self.final_peak_list_mz_max_vals = []
        self.final_peak_list_mz_min_vals = []
        self.final_peak_list_peak_rt_start_vals = []
        self.final_peak_list_peak_rt_end_vals = []
        self.mz_index = []
        self.rt_index = []
        self.mz_min_index = []
        self.mz_max_index = []
        self.rt_min_index = []
        self.rt_max_index = []
        self.coef_over_areas = []
        self.asym_gaus_fit_vals = []
        self.similarity_adj_mz_vals = []
    # adjacent similarity in mz domain. List of lists.
    def addSimilarityAdjMZVals(self,siIn):
        self.similarity_adj_mz_vals.append(siIn)
    def addAsymGausFitVal(self,agIn):
        self.asym_gaus_fit_vals.append(agIn)
    def addCoefOverArea(self,caIn):
        self.coef_over_areas.append(caIn)
    def addMZMinIndex(self,mzI):
        self.mz_min_index.append(mzI)
    def addMZMaxIndex(self,mzI):
        self.mz_max_index.append(mzI)
    def addRTMinIndex(self,rtI):
        self.rt_min_index.append(rtI)
    def addRTMaxIndex(self,rtI):
        self.rt_max_index.append(rtI)
    def addMZIndex(self,mzI):
        self.mz_index.append(mzI)
    def addRTIndex(self,rtI):
        self.rt_index.append(rtI)
    def addRT(self,rt):
        self.final_peak_list_rt_vals.append(rt)
    def addMZ(self,mz):
        self.final_peak_list_mz_vals.append(mz)
    def addHeight(self,height):
        self.final_peak_list_peak_height_vals.append(height)
    def addMaxMZ(self,maxMZ):
        self.final_peak_list_mz_max_vals.append(maxMZ)
    def addMinMZ(self,minMZ):
        self.final_peak_list_mz_min_vals.append(minMZ)
    def addMinRT(self,minRT):
        self.final_peak_list_peak_rt_start_vals.append(minRT)
    def addMaxRT(self,maxRT):
        self.final_peak_list_peak_rt_end_vals.append(maxRT)
