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
from os import path
sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
import easyIOmassspec.easyio as eio
from prototypenewpeakdetect.datamodel import data_point
import tqdm



def load_data_points_in_lists(dfr, absolute_intensity_thresh):
    """
    Fills lists of mz values, RT values, and intensities as well as list of all data points

    :param dfr: Data file reader object
    :param absolute_intensity_thresh: Intensity below which everything is thrown out.
    :return: Lists in this order mz_by_scan, inten_by_scan, rt, list_all_data_points
    """

    RT_MIN = 0.5
    RT_MAX = 1.5
    MZ_MIN = 136.0474 - 5
    MZ_MAX = 136.0474 + 5

    mz_by_scan = []
    inten_by_scan = []
    rt = []
    list_all_data_points = []

    count = 0
    mz, inten = dfr.get_next_scan_mzvals_intensities()

    scan_index_count = 0
    while mz != None:
        # while (count<100):
        # line below is to skip chemical noise in data Change back to aboveline HERE
        # if (count>200)and(count<1200):
        cur_rt = dfr.get_rt_from_scan_num(count)
        if (cur_rt < (RT_MAX * 60.0)) and (cur_rt > (RT_MIN * 60.0)):
            sys.stdout.write("\r" + str(count))
            sys.stdout.flush()
            rt.append(cur_rt)
            cur_mz_arr = []
            cur_inten_arr  = []
            for i in range(len(mz)):
                # For testing HERE
                if (mz[i] > MZ_MIN) and (mz[i] < MZ_MAX):
                    # if (mz[i]>0.0)and(mz[i]<1000.0):
                    if inten[i] < absolute_intensity_thresh:
                        continue
                    cur_dp = data_point.DataPoint(scan_index_count, i, mz[i], inten[i])
                    list_all_data_points.append(cur_dp)
                    cur_mz_arr.append(mz[i])
                    cur_inten_arr.append(inten[i])
            scan_index_count += 1
            mz_by_scan.append(cur_mz_arr)
            inten_by_scan.append(cur_inten_arr)

        mz, inten = dfr.get_next_scan_mzvals_intensities()

        count += 1

    return mz_by_scan, inten_by_scan, rt, list_all_data_points

def main():
    fstr = "/Users/owenmyers/Desktop/Data/DC_010814_StandardsMixTest1_34StandardMix_01.mzXML"
    dfr = eio.DataFileReader(fstr, True)
    absolute_intensity_thresh = 2000
    mz_by_scan, inten_by_scan, rt_arr, list_all_data_points = load_data_points_in_lists(dfr, absolute_intensity_thresh)

    mz_arr = []
    int_arr = []
    tot_int_arr = []
    points_in_scan_arr = []
    for i in tqdm.tqdm(range(len(mz_by_scan))):
        points_in_scan_arr.append(len(mz_by_scan[i]))
        tot_int_arr.append(sum(inten_by_scan[i]))
        for j in range(len(mz_by_scan[i])):
            mz_arr.append(mz_by_scan[i][j])
            int_arr.append(inten_by_scan[i][j])

    #rt_arr = [1.0,10.0,20.0]
    #scan_index = [0,2,4]
    #act_scan_num =[1,3,5]
    #tot_int_arr = [int_arr[0]+int_arr[1],int_arr[2]+int_arr[3],int_arr[4]]
    #points_in_scan_arr = [2,2,2]

    num_scans = len(rt_arr)

    dfwrite = eio.DataFileWriter('test_out',num_scans)

    dfwrite.set_write_points_in_scan_arr(points_in_scan_arr )
    dfwrite.set_mz_values(mz_arr)
    dfwrite.set_intensity_values(int_arr)
    #dfwrite.set_scan_index(scan_index)
    dfwrite.set_scan_index()
    dfwrite.set_rt_values(rt_arr)
    dfwrite.set_tot_intensity_arr(tot_int_arr)
    #dfwrite.set_act_scan_num(act_scan_num)

    dfwrite.closeWriter()

if __name__ == "__main__":
    main()
