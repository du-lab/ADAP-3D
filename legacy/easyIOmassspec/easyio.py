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

import ioMZXML
import ioCDF
import numpy as np
import warnings
try:
    import netCDF4
except ImportError:
    warnings.warn("Warning: netCDF4 is not available on this machine. You will not be able to write to CDF (you should still be able to read CDF)")


class DataFileReader():
    def __init__(self,file_and_path_str,is_profile_bool):
        self.is_profile_bool = is_profile_bool
        self.file_str = file_and_path_str

        # look at the extension to find the file type. 
        # this is just a string. can either be:
        # 1) mzxml
        # 2) mzml
        # 3) cdf
        location_last_period = self.file_str[::-1].find('.')
        self.type_of_file = self.file_str[-location_last_period : ]

        if self.type_of_file == 'mzXML':
            self.file_handler = ioMZXML.EasyHandler(self.file_str,self.is_profile_bool)
        elif self.type_of_file == 'mzml':
            print ("Do not have mzml reader implemented")
            exit()
        elif (self.type_of_file == 'cdf') or (self.type_of_file == 'CDF'):
            self.file_handler = ioCDF.netCDFReadHandler(self.file_str)

    def get_rt_from_scan_num(self,scan_num):
        return self.file_handler.get_rt_from_scan_num(scan_num)
        
    def get_num_scans(self):
        return self.file_handler.get_num_scans()

    def get_next_scan_mzvals_intensities(self):
        """ Returns the arrays of the mz values and the intensities (in that order) of the next
        scan. If this is the first time this function is called it just gets that information for
        the first scan."""

        return self.file_handler.get_next_scan_mzvals_intensities()

    def set_cur_mzvalues_and_intensities(self, act_scan_number):
        self.file_handler.set_cur_mzvalues_and_intensities(act_scan_number)

    def get_mzvalues_and_intensities_from_scan_number(self, scan_number):
        return self.file_handler.get_mzvalues_and_intensities_from_scan_number(scan_number)

    def set_cur_scan(self, scan_number):
        self.file_handler.set_cur_scan(scan_number)

    def get_act_scan_num(self,scan_number):
        return self.file_handler.get_act_scan_num(scan_number)

    def get_cur_scan_num_from_act(self,act_scan_number):
        """ essentially this returns the index of teh actual scan number if you have an ordered list of
        the scans"""
        return self.file_handler.get_cur_scan_num_from_act(act_scan_number)

    def get_mzvals(self):
        return self.file_handler.get_mzvals()
    def get_intensities(self):
        return self.file_handler.get_intensities()
    def plot_cur_scan(self):
        self.file_handler.plot_cur_scan()
    def get_mzvals_intensities_by_scan_num(self,scan_number):
        return self.file_handler.get__mzvals_intensities_by_scan_num(scan_number)

    def get_TIC(self):
        num_scan = self.get_num_scans()
        TIC = np.array([])
        rt = np.array([])
        for i in range(num_scan):
            cur_rt = self.get_rt_from_scan_num(i)
            print ("cur_rt" + str(cur_rt))
            rt = np.append(rt,cur_rt)

            cur_mz, cur_intensities = self.get_next_scan_mzvals_intensities()
            TIC = np.append(TIC,cur_intensities.sum() )

        return TIC,rt


# initialize
# set the  number of points in every scan arr
# With the above we can set the dimension of

class DataFileWriter():
    def __init__(self,file_and_path,num_scans):
        # number of scans to be written 
        self.num_scans = num_scans
        # The file string that is passed in will not have the extension
        self.file_str = file_and_path 

        # the writer file object
        #self.f_obj = netcdf.netcdf_file(self.file_str+'.CDF', 'w')
        #self.f_obj = netCDF4.Dataset(self.file_str+'.CDF', mode='w',format='NETCDF4_CLASSIC')
        self.f_obj = netCDF4.Dataset(self.file_str+'.CDF', mode='w',format='NETCDF3_64BIT')

        # the points_in_scan array contains the number of data points per scan. Each element of the
        # aray corresponds to a particular scan
        self.points_in_scan_arr = np.array([])

        self.scan_index_arr = []

        self.how_many_scans = 0

        # some booleans to make sure that everything has been done before closing the file
        self.have_set_write_points_in_scan = False
        self.have_set_mz_values = False
        self.have_set_intensities = False
        self.have_set_scan_index = False
        self.have_set_rt_values = False
        self.have_set_tot_intensities = False

    def set_write_points_in_scan_arr(self,to_set_arr):
        """
        Also sets the scan index array which is just a simple incremental array
        """
        self.points_in_scan_arr = to_set_arr
        self.how_many_scans = len(self.points_in_scan_arr)

        print ("self.how_many_scans")
        print (self.how_many_scans)
        print ("len(to_set_arr)")
        print (len(to_set_arr))

        self.f_obj.point_count = 'How many points in scan'
        self.f_obj.createDimension('point_count',len(to_set_arr))
        point_count = self.f_obj.createVariable('point_count','i',('point_count',))
        point_count[:] = self.points_in_scan_arr
        point_count.units = 'count'

        #self.scan_index_arr = range(1,self.how_many_scans + 1,1)

        self.have_set_write_points_in_scan = True

    def set_mz_values(self,to_set_arr):
        self.f_obj.mass_values = 'all mz values'
        self.f_obj.createDimension('mass_values',len(to_set_arr))
        mass_values = self.f_obj.createVariable('mass_values','d',('mass_values',))
        mass_values[:] = to_set_arr
        mass_values.units = 'MZ'

        self.have_set_mz_values = True


    def set_intensity_values(self,to_set_arr):
        self.f_obj.intensity_values = 'all intensities'
        self.f_obj.createDimension('intensity_values',len(to_set_arr))
        intensity_values = self.f_obj.createVariable('intensity_values','d',('intensity_values',))
        #intensity_values[:] = [000.0,110.0,1000.0,2000.0,2500.0,3000.0,2500.0,2000.0,1000.0,000.0]
        print "len(to_set_arr): "+str(len(to_set_arr))
        intensity_values[:] = to_set_arr
        intensity_values.units = 'intensity'

        self.have_set_intensities = True

    #def set_scan_index(self,scan_index_arr):
    def set_scan_index(self):
        """
        The "scan_index" is where the magic happens in MZmine. this array is a set of numbers
        corresponding to the index of the mz_values array telling you where the beginning of each scan
        is.
        """
        
        print ("setting the scan index")
        #scan_index_arr = [long(0)]
        scan_index_arr = np.array([0])

        cz = 0
        for i,j in enumerate(self.points_in_scan_arr):
            #print "i"
            #print i
            #print "j"
            #print j

            if j == 0:
                cz+=1

            
            if i == 0: 
                continue
            #scan_index_arr.append(long(j+scan_index_arr[-1]))
            scan_index_arr = np.append(scan_index_arr ,j+scan_index_arr[-1])

        print ("scan_index_arr")
        print (scan_index_arr)
        print ("zeros")
        print (cz)



        self.f_obj.scan_index= 'idecies of full mass list where scan starts??'
        self.f_obj.createDimension('scan_index',len(scan_index_arr))
        scan_index = self.f_obj.createVariable('scan_index','i',('scan_index',))
        scan_index[:] = scan_index_arr.astype(long)
        scan_index.units = 'scan_number'

        self.have_set_scan_index = True


    def set_rt_values(self,to_set_arr):
        self.f_obj.scan_acquisition_time= 'time of scan beginings'
        self.f_obj.createDimension('scan_acquisition_time',len(to_set_arr))
        scan_acquisition_time = self.f_obj.createVariable('scan_acquisition_time','d',('scan_acquisition_time',))
        scan_acquisition_time[:] = to_set_arr
        scan_acquisition_time.units = 'seconds'

        self.have_set_rt_values = True

    
    # XCMS needs this one
    def set_tot_intensity_arr(self,to_set_arr):

        self.f_obj.total_intensity= 'Total intensity of each scan'
        self.f_obj.createDimension('total_intensity',len(to_set_arr))
        total_intensity = self.f_obj.createVariable('total_intensity','d',('total_intensity',))
        total_intensity[:] = to_set_arr

        total_intensity.units = 'intensity'

        self.have_set_tot_intensities = True

    #def set_act_scan_num(self,to_set_arr):
    #    self.f_obj.actual_scan_number = 'The actual scan number'
    #    self.f_obj.createDimension('actual_scan_number',len(to_set_arr))
    #    total_intensity = self.f_obj.createVariable('actual_scan_number','i',('actual_scan_number',))
    #    total_intensity[:] = to_set_arr

    #    total_intensity.units = 'intensity'

    def closeWriter(self):

        self.f_obj.close()











        
