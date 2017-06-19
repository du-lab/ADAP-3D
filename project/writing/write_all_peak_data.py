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

# This is the main body of the program which detects peaks using similarity measurements between
# adjacent (in m/z) slices of peaks.

import os
import shutil
import sys
from os import path

sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
from plotting import peak_plotting

def write(result_dir_str,all_peak_objects,int_matrix,rt,scans):
    """
    Write the intensities, retention times, and scan numbers of each peak in a seperate file.
    Might be useful for a variaty of things but this is done primaraly so the script 
    jython_turn_results_into_mzmine2_xml.py can be run. This script makes a file that
    MZmine 2 can read so the results can be imported into MZmine 2.
    
    :param result_dir_str: Where will the data be writen (directory)?
    :param all_peak_objects: List of peak objects (see datamodel peak.py)
    :param int_matrix: Matrix of intensities. Rows are m/z, columns are RT.
    :param rt: Array of retention times.
    :param scans: Array of the scan numbers.
    """

    list_of_all_peak_mz_and_rt = []

    peakDataDirStr = result_dir_str+"PeakAllData"
    if os.path.exists(peakDataDirStr):
        shutil.rmtree(peakDataDirStr)
    os.mkdir(peakDataDirStr)

    count = 0

    for k in range(len(all_peak_objects)):
        cur_peak = all_peak_objects[k]

        if not peak_plotting.didWePlotThisPeakAlready(list_of_all_peak_mz_and_rt, cur_peak.mz, cur_peak.rt):

            f = open(peakDataDirStr+"/"+str(count)+"_peak.dat", "w")
            f.write("rt, intensity, scan\n")

            # get intensities
            intensities = int_matrix[cur_peak.mzIndex, cur_peak.rtMinIndex:cur_peak.rtMaxIndex]
            cur_rt = rt[cur_peak.rtMinIndex:cur_peak.rtMaxIndex]
            cur_scans = scans[cur_peak.rtMinIndex:cur_peak.rtMaxIndex]

            for i in range(len(cur_rt)):
                f.write(str(cur_rt[i]) + ", " + str(intensities[i]) + ", " + str(cur_scans[i]) + "\n")

            f.close()

            list_of_all_peak_mz_and_rt.append((cur_peak.mz, cur_peak.rt))

            count+=1

        if cur_peak.hasIsotope:
            for kappa in range(len(cur_peak.isotopeList)):
                cur_iso_peak = cur_peak.isotopeList[kappa]

                if peak_plotting.didWePlotThisPeakAlready(list_of_all_peak_mz_and_rt, cur_iso_peak.mz, cur_iso_peak.rt):
                    continue

                f = open(peakDataDirStr+"/"+str(count)+"_peak.dat", "w")
                f.write("rt, intensity\n")

                # get intensities
                intensities = int_matrix[cur_iso_peak.mzIndex, cur_iso_peak.rtMinIndex:cur_iso_peak.rtMaxIndex]
                cur_rt = rt[cur_iso_peak.rtMinIndex:cur_iso_peak.rtMaxIndex]
                cur_scans = scans[cur_iso_peak.rtMinIndex:cur_iso_peak.rtMaxIndex]

                for i in range(len(cur_rt)):
                    f.write(str(cur_rt[i]) + ", " + str(intensities[i]) + ", " + str(cur_scans[i]) + "\n")

                f.close()

                list_of_all_peak_mz_and_rt.append((cur_iso_peak.mz, cur_iso_peak.rt))

                count+=1
