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
import argparse
from easyIOmassspec import easyio


def get_start_point(mz1,division):

    mult = int(round(mz1/division))

    return mult*division


def get_interpolated_intensity(mz1, mz2, mz_in, int1, int2):
    slope = (int2-int1)/(mz2-mz1)
    b = int2-slope*mz2

    new_int = mz_in*slope+b

    return new_int

def main():

    parser = argparse.ArgumentParser()
    # data file to be resampled
    parser.add_argument('-f',action='store',dest = 'f',type=str, required = True)
    # Name of resampled data file output.
    parser.add_argument('-o',action='store',dest = 'o',type=str, required = False,
                        default='resampled_data_file.CDF')
    smallest_diff = 0.0004

    inargs = parser.parse_args()

    data_file_str = inargs.f
    output_file_str = inargs.o

    data_file_reader = easyio.DataFileReader(data_file_str, True)

    mz_by_scan = []
    inten_by_scan = []
    rt = []
    scan_numbers = []
    #list_all_data_points = []

    count = 0
    mz,inten = data_file_reader.get_next_scan_mzvals_intensities()
    while mz != None:
        # while (count<100):
        # line below is to skip chemical noise in data Change back to above line HERE
        # if (count>200)and(count<1200):
        print("\r" + str(count))
        #sys.stdout.write("\r" + str(count))
        #sys.stdout.flush()
        mz_by_scan.append(mz)
        inten_by_scan.append(inten)
        rt.append(data_file_reader.get_rt_from_scan_num(count))
        scan_numbers.append(data_file_reader.get_act_scan_num(count))
        #for i in range(len(mz)):
        # # For testing HERE
        # if (mz[i] > 0.0) and (mz[i] < 1000.0):
        #     # if (mz[i]>0.0)and(mz[i]<1000.0):
        #     if inten[i] < absolute_intensity_thresh:
        #         continue
        #     cur_dp = data_point.DataPoint(count, i, mz[i], inten[i])
        #     list_all_data_points.append(cur_dp)

        mz, inten = data_file_reader.get_next_scan_mzvals_intensities()

        count += 1

#    print "Finding smallest difference: "
#
#    # search scans for the smallest mz spacing
#    smallest_diff = 1.0
#    for i in range(len(mz_by_scan)):
#        cur_scan = mz_by_scan[i]
#        sys.stdout.write("\r" + str(i) + " out of " + str(len(mz_by_scan)))
#        sys.stdout.flush()
#
#        for j in range(len(cur_scan)-1):
#            cur_mz_diff = abs(cur_scan[j+1] - cur_scan[j])
#            if cur_mz_diff<smallest_diff:
#                smallest_diff = cur_mz_diff

    print "smallest mz difference is: " +str(smallest_diff)

    new_mz_by_scan = []
    new_inten_by_scan = []

    print "Making resampled data"
    for i in range(len(mz_by_scan)):

        print( str(i) + " out of " + str(len(mz_by_scan)))
        print "\n"
        #sys.stdout.write("\r" + str(i) + " out of " + str(len(mz_by_scan)))
        #sys.stdout.flush()

        cur_scan_mz = mz_by_scan[i]
        cur_scan_int = inten_by_scan[i]

        # the point in mz we are currently working with
        working_mz = 0.0
        cur_resampled_mz = []
        cur_resampled_inten = []
        for j in range(len(cur_scan_mz)-1):
            # check gap between points 1
            if abs(cur_scan_mz[j+1]-cur_scan_mz[j])>10.0*smallest_diff:
                continue
            # check gap between points 2
            if abs(working_mz-cur_scan_mz[j])>10.0*smallest_diff:
                working_mz = get_start_point(cur_scan_mz[j],smallest_diff)
                cur_resampled_mz.append(working_mz)
                working_intensity = get_interpolated_intensity(cur_scan_mz[j],
                                                               cur_scan_mz[j+1],
                                                               working_mz,
                                                               cur_scan_int[j],
                                                               cur_scan_int[j+1])

                cur_resampled_inten.append(working_intensity)

            while working_mz+smallest_diff <= cur_scan_mz[j+1]:
                working_mz+=smallest_diff
                cur_resampled_mz.append(working_mz)
                working_intensity = get_interpolated_intensity(cur_scan_mz[j],
                                                           cur_scan_mz[j+1],
                                                           working_mz,
                                                           cur_scan_int[j],
                                                           cur_scan_int[j+1])

                cur_resampled_inten.append(working_intensity)

        new_inten_by_scan.append(cur_resampled_inten)
        new_mz_by_scan.append(cur_resampled_mz)

    print "resampling compleate"

    print "type(new_inten_by_scan): " + str(type(new_inten_by_scan))

    print "preparing data for writing"

    points_in_scan_arr = []
    all_inten_arr = []
    all_mz_arr = []
    tic = []
    for i in range(len(new_mz_by_scan)):
        points_in_scan_arr.append(len(new_mz_by_scan[i]))
        all_inten_arr+=new_inten_by_scan[i]
        all_mz_arr+=new_mz_by_scan[i]
        tic.append(sum(new_inten_by_scan[i]))

    print "done preparing data for writing"
    print "writing"

    dfwrite = easyio.DataFileWriter('resampled_profile_data', len(points_in_scan_arr))

    dfwrite.set_write_points_in_scan_arr(points_in_scan_arr )
    dfwrite.set_mz_values(all_mz_arr)
    dfwrite.set_intensity_values(all_inten_arr)
    #dfwrite.set_scan_index(scan_index)
    dfwrite.set_scan_index()
    dfwrite.set_rt_values(rt)
    dfwrite.set_tot_intensity_arr(tic)
    #dfwrite.set_act_scan_num(act_scan_num)

    dfwrite.closeWriter()

    print "done writing"


if __name__ == '__main__':
    main()

