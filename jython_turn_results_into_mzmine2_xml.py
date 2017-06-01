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

import xml.etree.ElementTree as ET
from xml.dom import minidom
import time
from java.io import ByteArrayOutputStream
from java.io import DataOutputStream
import base64
import argparse

def clean_up_output(elem):
    """Return an XML string for the Element with proper indents and returns. """

    rough_string = ET.tostring(elem, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    #return reparsed.toprettyxml(indent="   ")
    return reparsed.toxml()


def playWithByteStream():
    baos = ByteArrayOutputStream()
    ds = DataOutputStream(baos)
    ds.write(1)
    ds.flush()
    ds.write(2)
    ds.flush()

    print baos.toByteArray()

    print "type(baos.toByteArray()): " + str(type(baos.toByteArray()))

    print "base 64 encoded: "
    print base64.b64encode(baos.toByteArray())


def getProperBase64EncodingOfIntArr(cur_scan_arr):
    baos = ByteArrayOutputStream()
    ds = DataOutputStream(baos)
    for i in cur_scan_arr:
        ds.writeInt(i)
        ds.flush()

    return base64.b64encode(baos.toByteArray())


def getProperBase64EncodingOfFloatArr(cur_scan_arr):
    baos = ByteArrayOutputStream()
    ds = DataOutputStream(baos)
    for i in cur_scan_arr:
        ds.writeFloat(i)
        ds.flush()

    return base64.b64encode(baos.toByteArray())

def get_data_file_str(lines):
    dfs = ''
    for i in lines:
        if "Data file" in i:
            dfs = i[i.find(":")+2:]
            break

    reversed_dfs = dfs[::-1]
    sliced_rev_dfs = reversed_dfs[:reversed_dfs.find("/")]
    return sliced_rev_dfs[::-1].strip('\n')

def get_index_of_mz_values(s):
    splited = s.split(',')
    for i in range(len(splited)):
        if 'row m/z' in splited[i]:
            return i

def get_index_of_rt_values(s):
    splited = s.split(',')
    for i in range(len(splited)):
        if 'row retention time' in splited[i]:
            return i

def get_index_of_height_values(s):
    splited = s.split(',')
    for i in range(len(splited)):
        if 'Peak height' in splited[i]:
            return i

def find_area(rt,inten):
    area = 0.0
    for i in range(len(rt)-1):
        left_h = inten[i]
        right_h = inten[i+1]
        diff = abs(right_h-left_h)
        triangle_area = 0.5*diff*(rt[i+1]-rt[i])
        rectangle_area = min(left_h,right_h)*(rt[i+1]-rt[i])
        area+= triangle_area + rectangle_area

    return area

def main():
    parser = argparse.ArgumentParser()
    
    # Specify the directory containing results. This program needs several things so
    # it cant be a specific file has to the the un modified results directory.
    parser.add_argument('-d',action='store',dest = 'd',type=str, required = True)

    inargs = parser.parse_args()
    working_directory = inargs.d
    if working_directory[-1] != '/':
        working_directory += '/'

    misc_info_file = open(working_directory + "run_info.txt","r")
    misc_info_lines = misc_info_file.readlines()
    data_file_str = get_data_file_str(misc_info_lines)

    peak_list_file = open(working_directory + "combined_no_repeat_peak_list.csv","r")
    peak_list_lines = peak_list_file.readlines()
    number_of_peaks = len(peak_list_lines)-1

    index_of_mz_values = get_index_of_mz_values(peak_list_lines[0])
    index_of_rt_values = get_index_of_rt_values(peak_list_lines[0])
    index_of_height_values = get_index_of_height_values(peak_list_lines[0])

    file_name = data_file_str
    print "file name"
    print file_name

    number_of_peaks = number_of_peaks;
    method_name = 'Peak Detection: ADAP 3D profile'
    # This is in the structure:
    # parameter1: value1, parameter2: value2: ...
    method_parameters = 'param1: 0.0'

    top = ET.Element('peaklist')

    file_name_elem = ET.SubElement(top, 'pl_name')
    file_name_elem.text = file_name + ': Peaks detected'
    time_created_elem = ET.SubElement(top, 'created')
    time_created_elem.text = time.strftime("%d/%m/%Y") + ' ' + time.strftime("%H:%M:%S")
    number_of_peaks_elem = ET.SubElement(top, 'quantity')
    number_of_peaks_elem.text = str(number_of_peaks)

    applied_method_elem = ET.SubElement(top, 'applied_method')
    method_name_elem = ET.SubElement(applied_method_elem, 'method_name')
    method_name_elem.text = method_name
    method_parameters_elem = ET.SubElement(applied_method_elem, 'method_parameters')
    method_parameters_elem.text = method_parameters

    raw_file_elem = ET.SubElement(top, 'raw_file')
    raw_file_elem.text = file_name

    for k in range(len(peak_list_lines)-1):
        #print peak_list_lines[k]
        # get the file containing the peak shape data
        peak_shape_file = open(working_directory + 'PeakAllData/'+str(k)+'_peak.dat')
        peak_shape_file.readline()
        peak_shape_lines = peak_shape_file.readlines()
        cur_peak_shape_rt = []
        cur_peak_shape_intensity = []
        cur_peak_shape_scans = []
        cur_peak_best_scan = 0
        cur_peak_best_inten = 0.0
        for i in peak_shape_lines:
            cur_peak_shape_rt.append(float(i.split(',')[0]))
            cur_peak_shape_intensity.append(float(i.split(',')[1]))
            cur_peak_shape_scans.append(int(i.split(',')[2]))
            if float(i.split(',')[1])>cur_peak_best_inten:
                cur_peak_best_scan = int(i.split(',')[2])
                cur_peak_best_inten = float(i.split(',')[1])

        #######################################################################
        # begin rows -> this will need to be in a loop when you actually use it
        #######################################################################
        cur_line = peak_list_lines[k+1]
        splited_line = cur_line.split(',')

        cur_id_num = k
        cur_peak_mz = float(splited_line[index_of_mz_values])
        cur_peak_rt = float(splited_line[index_of_rt_values])*60.0
        cur_peak_height = float(splited_line[index_of_height_values])
        cur_peak_area = find_area(cur_peak_shape_rt,cur_peak_shape_intensity)
        cur_peak_status = "DETECTED"
        cur_peak_charge = 0

        cur_peak_shape_mz = []
        for i in range(len(peak_shape_lines)):
            cur_peak_shape_mz.append(cur_peak_mz)
        # we have to do these next two lines otherwise mzmine can not plot correctly
        cur_peak_shape_mz[0]-=0.0001
        cur_peak_shape_mz[-1]+=0.0001

        # This is for the best scan element which looks like it is just the scan number of the highest
        # intensity in the peak.
        cur_best_scan = cur_peak_best_scan
        # This will not be necessary for now -> should be -1 if no meaning
        cur_fragment_scan = -1
        # I don't know what this is. look into it. Might just be random label
        cur_quantity = len(cur_peak_shape_scans)

        cur_scan_arr = cur_peak_shape_scans
        cur_mz_arr = cur_peak_shape_mz
        cur_height_arr = cur_peak_shape_intensity
        for i in range(cur_quantity):
            cur_scan_arr64 = getProperBase64EncodingOfIntArr(cur_scan_arr)
            cur_mz_arr64 = getProperBase64EncodingOfFloatArr(cur_mz_arr)
            cur_height_arr64 = getProperBase64EncodingOfFloatArr(cur_height_arr)

        # Figure out correct format for these
        cur_scan_id = cur_scan_arr64
        cur_mz = cur_mz_arr64
        cur_height = cur_height_arr64

        row_elem = ET.SubElement(top, 'row', {'id': str(cur_id_num)})
        peak_elem = ET.SubElement(row_elem, 'peak', {
            'column_id': file_name,
            'mz': str(cur_peak_mz),
            'rt': str(cur_peak_rt),
            'height': str(cur_peak_height),
            'area': str(cur_peak_area),
            'status': cur_peak_status,
            'charge': str(cur_peak_charge)
        })

        best_scan_elem = ET.SubElement(peak_elem, 'best_scan')
        best_scan_elem.text = str(cur_best_scan)
        fragment_scan_elem = ET.SubElement(peak_elem, 'fragment_scan')
        fragment_scan_elem.text = str(cur_fragment_scan)

        mz_peaks_elem = ET.SubElement(peak_elem, 'mzpeaks', {'quantity': str(cur_quantity)})
        scan_id_elem = ET.SubElement(mz_peaks_elem, 'scan_id')
        scan_id_elem.text = str(cur_scan_id)
        mz_elem = ET.SubElement(mz_peaks_elem, 'mz')
        mz_elem.text = cur_mz
        height_elem = ET.SubElement(mz_peaks_elem, 'height')
        height_elem.text = cur_height

    # f = open('test_out.xml','w')
    # rough_string = ET.tostring(top , 'utf-8')
    # f.write(rough_string)
    # f.close()

    f = open(working_directory + 'full_peaklist.xml.mpl', 'w')
    f.write(clean_up_output(top))
    f.close()

if __name__ == '__main__':
    main()

# jython jython_turn_results_into_mzmine2_xml.py -d /Users/xdu4/Documents/Duxiuxia/Analysis/my_projects/adap-3d/results/    
