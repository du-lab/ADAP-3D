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

import sys
import mspy.parser_mzxml as pmzxml
import scipy as sp 
from matplotlib import pyplot as plt


class EasyHandler:
    """ This class just serves as a different interface for the parsing program writen by Martin
    Strohalm. Because of this it requires the mass package and rides entirely off of the hard work
    done in producing that package. 

    Important note is that the "scan_number" will refer to an incremented list of numbers from 0 to
    the number of scans. The "actual" scan number may not be incremented and will be referred to as
    the "act_scan_number"
    
    """

    def __init__(self,file_and_path_str, is_profile_bool):
        self.is_profile_bool = is_profile_bool

        self.file_str = file_and_path_str
        self.file_parser = pmzxml.parseMZXML(self.file_str)
        self.file_parser.load()
        self.scan_list = self.file_parser.scanlist()
        
        self.total_number_of_scans = len(self.scan_list.keys())
        
        
        self.scan_number_array = sp.array([])
        self.act_scan_number_array = sp.array([])

        self.cur_mzvals = sp.array([])
        self.cur_intensities = sp.array([])
        self.rt_arr = sp.array([])

        self.get_scan_nums_in_orderd_array()
        # lets sort them so we can go through them croologicaly
        self.act_scan_number_array.sort()
        self.act_scan_number_array = self.act_scan_number_array.astype(int)

        self.cur_scan = {}
        self.cur_num_scan_list = 0
        
    def get_one_scan_by_scan_index(self, scan_index):
        """
        This method is to replace the old one "get_mzvalues_and_intensities_from_scan_number()".
        """
        self.set_cur_mzvalues_and_intensities(self.act_scan_number_array[scan_index])
        
        scan = {}
        scan['mz'] = self.get_mzvals()
        scan['intensity'] = self.get_intensities()
        return scan

    def get_actual_scan_number_by_scan_index(self,scan_index):
        """
        This method gets the actual scan number from scan index. This method is to replace the method "get_act_scan_num()".
        """
        return self.act_scan_number_array[scan_index]    

    def get_scan_list(self):
        return self.scan_list

    def get_rt_from_scan_num(self,scan_num):
        scan = self.file_parser.scan(self.act_scan_number_array[scan_num])
        return scan.retentionTime

    def get_rt_arr(self):
        return self.rt_arr
        
    def get_num_scans(self):
        return len(self.scan_number_array)


    def get_next_scan_mzvals_intensities(self):
        """ Returns the arrays of the mz values and the intensities (in that order) of the next
        scan. If this is the first time this function is called it just gets that information for
        the first scan."""
        if (self.cur_num_scan_list >= len(self.act_scan_number_array)):
            return None,None

        act_cur_scan_num = self.act_scan_number_array[self.cur_num_scan_list]

        ms_level = self.file_parser._scanlist[act_cur_scan_num]['msLevel']

        while ms_level >= 2:
            self.cur_num_scan_list += 1
            if (self.cur_num_scan_list >= len(self.act_scan_number_array)):
                return None,None
            act_cur_scan_num = self.act_scan_number_array[self.cur_num_scan_list]
            ms_level = self.file_parser._scanlist[act_cur_scan_num]['msLevel']

        self.set_cur_scan(act_cur_scan_num )
        self.set_cur_mzvalues_and_intensities(act_cur_scan_num)

        self.cur_num_scan_list += 1

        return self.get_mzvals(), self.get_intensities()

    def set_cur_mzvalues_and_intensities(self, act_scan_number):
        self.set_cur_scan(act_scan_number)
        if self.is_profile_bool:
            if len(self.cur_scan.profile) != 0:
                self.cur_mzvals = self.cur_scan.profile[:,0]
                self.cur_intensities = self.cur_scan.profile[:,1]
        else: 
            self.cur_mzvals = sp.array([])
            self.cur_intensities = sp.array([])
            for i in self.cur_scan.peaklist:
                self.cur_mzvals = sp.append(self.cur_mzvals, i.mz)
                self.cur_intensities = sp.append(self.cur_intensities, i.intensity)

    def get_mzvalues_and_intensities_from_scan_number(self, scan_number):
        self.set_cur_mzvalues_and_intensities(self.act_scan_number_array[scan_number])
        return self.get_mzvals(), self.get_intensities()
        
        
    #def set_cur_mzvalues_and_intensities(self):
    #    self.mzvals = cur_scan.profile[:,0]
    #    self.intensities = cur_scan.profile[:,1]

    def set_cur_scan(self, scan_number):
        self.cur_scan = self.file_parser.scan(scan_number)

    def get_act_scan_num(self,scan_number):
        return self.act_scan_number_array[scan_number]

        

    def get_scan_index_from_actual_scan_number(self,act_scan_number):
        """ 
        This method retuns the scan index of an actual scan number. Note that the scan index starts at 0.
        """
        for i in range(len(self.act_scan_number_array)):
            if act_scan_number == self.act_scan_number_array[i]:
                return i
                
                
                

    def get_mzvals(self):
        return self.cur_mzvals
        
        
        
        
    def get_intensities(self):
        return self.cur_intensities
        
        
        
        
        
    def plot_cur_scan(self):
        plt.plot(self.mzvals,intensities)
        
    def get_scan_nums_in_orderd_array(self):
        if len(self.scan_number_array) == 0:
            count = 0
            #print self.scan_list
            for i in self.scan_list:
                self.act_scan_number_array = sp.append(self.act_scan_number_array,i)
                self.scan_number_array = sp.append(self.scan_number_array,count)
                #print("i['retentionTime']")
                #print(self.scan_list[i]['retentionTime'])
                self.rt_arr = sp.append(self.rt_arr,self.scan_list[i]['retentionTime'])
                count +=1

            self.rt_arr.sort()


    def get_mzvals_intensities_by_scan_num(self,scan_number):
        self.cur_scan = self.file_parser.scan(scan_number)

        self.set_cur_scan(cur_scan_num )
        self.set_cur_mzvalues_and_intensities(cur_scan_num)

        return self.get_mzvals(), self.get_intensities()


    def get_file_str(self):
        return self.file_str
