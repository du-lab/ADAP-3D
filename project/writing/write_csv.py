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

import sys
from os import path
sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
from plotting import peak_plotting

def write(result_dir_str, all_peak_objects):

    list_of_all_peak_mz_and_rt = []

    f = open(result_dir_str+"peaklist.csv", "w")
    f.write("row m/z,row retention time,Peak height,Peak m/z max,Peak m/z min,Peak RT start,Peak RT end\n")
    iso_file = open(result_dir_str+"isotope_peaklist.csv", "w")
    iso_file.write("row m/z,row retention time,Peak height,Peak m/z max,Peak m/z min,Peak RT start,Peak RT end\n")
    combined_no_repeat = open(result_dir_str+"combined_no_repeat_peak_list.csv","w")
    combined_no_repeat.write("row m/z,row retention time,Peak height,Peak m/z max,Peak m/z min,Peak RT start,Peak RT end\n")

    for k in range(len(all_peak_objects)):
        cur_peak = all_peak_objects[k]

        f.write(str(cur_peak.mz))
        f.write(',')
        f.write(str(cur_peak.rt))
        f.write(',')
        f.write(str(cur_peak.height))
        f.write(',')
        f.write(str(cur_peak.mzMax))
        f.write(',')
        f.write(str(cur_peak.mzMin))
        f.write(',')
        f.write(str(cur_peak.rtMin))
        f.write(',')
        f.write(str(cur_peak.rtMax))
        f.write("\n")

        if not peak_plotting.didWePlotThisPeakAlready(list_of_all_peak_mz_and_rt, cur_peak.mz, cur_peak.rt):
            combined_no_repeat.write(str(cur_peak.mz))
            combined_no_repeat.write(',')
            combined_no_repeat.write(str(cur_peak.rt))
            combined_no_repeat.write(',')
            combined_no_repeat.write(str(cur_peak.height))
            combined_no_repeat.write(',')
            combined_no_repeat.write(str(cur_peak.mzMax))
            combined_no_repeat.write(',')
            combined_no_repeat.write(str(cur_peak.mzMin))
            combined_no_repeat.write(',')
            combined_no_repeat.write(str(cur_peak.rtMin))
            combined_no_repeat.write(',')
            combined_no_repeat.write(str(cur_peak.rtMax))
            combined_no_repeat.write("\n")
            list_of_all_peak_mz_and_rt.append((cur_peak.mz, cur_peak.rt))

        if cur_peak.hasIsotope:
            for kappa in range(len(cur_peak.isotopeList)):
                cur_iso_peak = cur_peak.isotopeList[kappa]

                iso_file.write(str(cur_iso_peak.mz))
                iso_file.write(',')
                iso_file.write(str(cur_iso_peak.rt))
                iso_file.write(',')
                iso_file.write(str(cur_iso_peak.height))
                iso_file.write(',')
                iso_file.write(str(cur_iso_peak.mzMax))
                iso_file.write(',')
                iso_file.write(str(cur_iso_peak.mzMin))
                iso_file.write(',')
                iso_file.write(str(cur_iso_peak.rtMin))
                iso_file.write(',')
                iso_file.write(str(cur_iso_peak.rtMax))

                iso_file.write("\n")

                if not peak_plotting.didWePlotThisPeakAlready(list_of_all_peak_mz_and_rt, cur_iso_peak.mz, cur_iso_peak.rt):
                    combined_no_repeat.write(str(cur_iso_peak.mz))
                    combined_no_repeat.write(',')
                    combined_no_repeat.write(str(cur_iso_peak.rt))
                    combined_no_repeat.write(',')
                    combined_no_repeat.write(str(cur_iso_peak.height))
                    combined_no_repeat.write(',')
                    combined_no_repeat.write(str(cur_iso_peak.mzMax))
                    combined_no_repeat.write(',')
                    combined_no_repeat.write(str(cur_iso_peak.mzMin))
                    combined_no_repeat.write(',')
                    combined_no_repeat.write(str(cur_iso_peak.rtMin))
                    combined_no_repeat.write(',')
                    combined_no_repeat.write(str(cur_iso_peak.rtMax))
                    combined_no_repeat.write("\n")
                    list_of_all_peak_mz_and_rt.append((cur_iso_peak.mz, cur_iso_peak.rt))

