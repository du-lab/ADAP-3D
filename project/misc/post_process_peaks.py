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

def main():

    mz_tol = 0.000001
    rt_tol = 0.00001

    fin = open("peaklist.csv","r")
    fout = open("peaklist_after_postprocessing.csv","w")

    headerLine = fin.readline()
    fout.write(headerLine)

    lines = fin.readlines()

    mzList = []
    rtList = []

    print "line 0:"
    print lines[0]

    countDuplicates = 0
    for i in lines:
        curMZ = float(i[:i.find(",")])
        curRT = float(i[i.find(",")+1:][:i[i.find(",")+1:].find(",")])

        foundDuplicate = False
        for j in range(len(mzList)):
            if (curMZ<(mzList[j]+mz_tol)) and (curMZ>(mzList[j]-mz_tol)):
                if (curRT<(rtList[j]+rt_tol)) and (curRT>(rtList[j]-rt_tol)):

                    countDuplicates+=1
                    foundDuplicate = True
        
        if not foundDuplicate:
            fout.write(i)

            mzList.append(curMZ)
            rtList.append(curRT)
        
    print "number duplicates: " + str(countDuplicates)

    fin.close()
    fout.close()
if __name__ == "__main__":
    main()
