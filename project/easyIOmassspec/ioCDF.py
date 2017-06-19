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

from scipy.io import netcdf
import numpy as np


class netCDFReadHandler:
    """ 
    This class contains methods that handle mass spectra in netCDF format.
    """

    def __init__(self, data_file_str_w_path):
        self.data_file_full_name = data_file_str_w_path
        self.ncid = netcdf.netcdf_file(self.data_file_full_name, 'r')
                
        
        total_intensity_variable = self.ncid.variables['total_intensity']
        self.total_intensity = total_intensity_variable[:].copy()
        
        scan_acquisition_time_variable = self.ncid.variables['scan_acquisition_time']
        self.scan_acquisition_time = scan_acquisition_time_variable[:].copy()
        
        point_count_variable = self.ncid.variables['point_count']
        self.point_count = point_count_variable[:].copy()
        
        mass_values_variable = self.ncid.variables['mass_values']
        self.mass_values = mass_values_variable[:].copy()
        
        intensity_values_variable = self.ncid.variables['intensity_values']
        self.intensity_values = intensity_values_variable[:].copy()

        index_of_intensity_values_variable = self.ncid.variables['scan_index']
        self.index_of_intensity_values = index_of_intensity_values_variable[:].copy()
        
        self.total_number_of_scans = len(self.scan_acquisition_time)
        self.total_number_of_data_points = len(self.mass_values)

        #self.ncid.close()


        self.cur_index = 0

    def get_rt_from_scan_num(self,scan_num):
        return self.scan_acquisition_time[scan_num]



    def get_next_scan_mzvals_intensities(self):
        """ Returns the arrays of the mz values and the intensities (in that order) of the next
        scan. If this is the first time this function is called it just gets that information for
        the first scan."""
        mzandints = self.get_one_scan_by_scan_index(self.cur_index)
        self.cur_index+=1

        mzs = mzandints['mz']

        ints = mzandints['intensity']

        if len(mzs)==0:
            return None, None

        return mzs,ints
        
    def get_file_name(self):
        return self.data_file_full_name
        
        
    def get_TIC(self):
        """
        This method extracts the TIC.
        
        Return a dictionary consisting of the scan_acquisition_time and the total_intensity.
        """
        
        tic = {}
        tic['scan_acquisition_time'] = self.scan_acquisition_time
        tic['total_intensity'] = self.total_intensity
        
        return tic
        
        
    def get_one_scan_by_scan_index(self, scan_index):
        """
        This method extracts one scan from the raw data.
        """
        

        if scan_index == 0:
            start_index = 0
            end_index = self.point_count[0]
        elif scan_index == 1:
            start_index = self.point_count[0]
            end_index = start_index + self.point_count[1]
        else:
            start_index = np.sum(self.point_count[0:scan_index])
            end_index = start_index + self.point_count[scan_index]

        if (start_index>len(self.mass_values)-1) or (end_index>len(self.mass_values)-1):
            scan = {}
            scan['mz'] = []
            scan['intensity'] = []
            return scan
        
        scan = {}
        scan['mz'] = self.mass_values[start_index:end_index]
        scan['intensity'] = self.intensity_values[start_index:end_index]
        
        return scan

    def get_act_scan_num(self, scan_number):
        """
        In CDF files the scans are just incrementally numbered so this is just a trick
        so you can call this function correctly from inside easyio.py regardless of
        the file type
        
        :param scan_number: integer referring to the index
        :return: scan_number
        """

        return scan_number
        

