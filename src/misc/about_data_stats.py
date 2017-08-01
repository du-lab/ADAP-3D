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

import pylab as pl
import sys
sys.path.append('../easyIOmassspec')
import easyio

class dataPoint():
    def __init__(self,scan_index,mz_index,mz,intensity):
        self.scan_index = scan_index
        self.mz_index = mz_index
        self.intensity = intensity
        self.mz = mz

def main():
    df_str = "/Users/owenmyers/Desktop/Data/DC_010814_StandardsMixTest1_34StandardMix_01.mzXML"
    dfr = easyio.DataFileReader(df_str,True)

    absoluteIntensityThresh = 0.0

    mz,inten = dfr.get_next_scan_mzvals_intensities()
    count = 0
    mz_by_scan = []
    inten_by_scan = []
    rt = []

    listAllDP = []
    all_inten = []
    while mz!=None:
    #while (count<100):
        # line below is to skip chemical noise in data Change back to aboveline HERE
        #if (count>200)and(count<1200):
        sys.stdout.write("\r"+str(count))
        sys.stdout.flush()
        mz_by_scan.append(mz)
        inten_by_scan.append(inten)
        rt.append(dfr.get_rt_from_scan_num(count))
        for i in range(len(mz)):
            # For testing HERE
            if (mz[i]>0.0)and(mz[i]<1000.0):
            #if (mz[i]>0.0)and(mz[i]<1000.0):
                if inten[i]<absoluteIntensityThresh:
                    continue
                cur_dp = dataPoint(count,i,mz[i],inten[i])
                all_inten.append(inten[i])
                listAllDP.append(cur_dp)

        mz,inten = dfr.get_next_scan_mzvals_intensities()

        count+=1

    pl.hist(all_inten,bins = 1000)
    pl.show()

if  __name__ == "__main__":
    main()
