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

import shutil
import os
import matplotlib
matplotlib.use('Agg')
import pylab as pl

##### Plot styles#########
markerfacecolor = 'k'
markeredgecolor = 'k'
linestyle = 'None'
markersize = 3
alpha = 1.0
capsize = 2
elinewidth = 0.4
markeredgewidth = 0.25
linewidth = 0.5
#####################
font = {'family': 'sans',
        'weight': 'normal',
        'size': 9}

matplotlib.rc('font', **font)


def didWePlotThisPeakAlready(listOfAllPeakMZAndRT,mz,rt):
    """
    Checks to see if the peak at mz and rt has already been plotted. If it has not then it is added to the list
    and that list is returned.
    :param listOfAllPeakMZAndRT: List of tuples containing (mz,rt) of peaks that have already been plotted
    :param mz: m/z value of peak being checked
    :param rt: rt value of the peak being checked
    :return: returns listOfAllPeakMZAndRT with the current peak added if it was not already present. 
    """
    mzTol = 0.001
    rtTol = 0.01
    inList = False
    for i in range(len(listOfAllPeakMZAndRT)):
        cur_mz = listOfAllPeakMZAndRT[i][0]
        cur_rt = listOfAllPeakMZAndRT[i][1]

        if ((mz-mzTol)<cur_mz) and ((mz+mzTol)>cur_mz):
            if ((rt-rtTol)<cur_rt) and ((rt+rtTol)>cur_rt):
                inList = True

    return inList

def cm2inch(value):
    return value/2.54

def plot_these_peaks(allPeakObjects,int_matrix,rt,plotDirStr,listOfAllPeakMZAndRT,countRepeats):

    for k in range(len(allPeakObjects)):
        curPeak = allPeakObjects[k]
        mzIndex = allPeakObjects[k].mzIndex
        rtIndex = allPeakObjects[k].rtIndex
        mzMaxIndex = allPeakObjects[k].mzMaxIndex
        mzMinIndex = allPeakObjects[k].mzMinIndex
        rtMaxIndex = allPeakObjects[k].rtMaxIndex
        rtMinIndex = allPeakObjects[k].rtMinIndex

        cur_peak_mz = allPeakObjects[k].mz
        cur_peak_rt = allPeakObjects[k].rt

        if didWePlotThisPeakAlready(listOfAllPeakMZAndRT, cur_peak_mz, cur_peak_rt):
            countRepeats += 1
            continue
        listOfAllPeakMZAndRT.append((cur_peak_mz, cur_peak_rt))

        # show beyond the peak
        eicMin = rtMinIndex - 40
        eicMax = rtMaxIndex + 40

        # print "eicMin before " + str(eicMin )
        # print "eicMax before " + str(eicMax )
        # print "len(pl.squeeze(pl.asarray(int_matrix[0,:].todense()))) " + str(len(pl.squeeze(pl.asarray(int_matrix[0,:].todense()))))

        lenIntMatrix = len(pl.squeeze(pl.asarray(int_matrix[0, :].todense())))

        if eicMin < 0: eicMin = 0
        if eicMax > lenIntMatrix: eicMax = lenIntMatrix
        # print "eicMin after " + str(eicMin )
        # print "eicMax after " + str(eicMax )

        cur_eic_int = pl.squeeze(pl.asarray(int_matrix[mzIndex, eicMin:eicMax].todense()))
        peak_only_eic_int = pl.squeeze(pl.asarray(int_matrix[mzIndex, rtMinIndex:rtMaxIndex].todense()))

        fig = matplotlib.pyplot.gcf()
        fig.set_size_inches(cm2inch(8.9), cm2inch(5.5))
        ax = fig.add_subplot(111)

        fig.subplots_adjust(bottom=0.14, left=0.12, right=0.98, top=0.92)
        # Set axis labels and locations
        ax.xaxis.set_label_coords(0.5, -0.1)
        ax.yaxis.set_label_coords(-0.085, 0.5)

        ax.yaxis.get_major_formatter().set_powerlimits((0, 1))
        ax.ticklabel_format(style='sci', useOffset=False)

        ax.set_xlabel("rt")
        ax.set_ylabel("intensity")
        # ax.set_ylim([most_min,most_max])
        # ax.set_xlim([min_rt_arr[0],min_rt_arr[-1]])

        # ax.plot(min_rt_arr,cur_eic_int,c='k',linewidth=linewidth,marker='o',markersize=markersize)
        # ax.plot(pl.arange(eicMin,eicMax),cur_eic_int,c='k',linewidth=linewidth,marker='o',markersize=markersize)
        ax.plot(rt[eicMin:eicMax], cur_eic_int, c='k', linewidth=linewidth, marker='o', markersize=markersize)
        # ax.fill_between(peak_only_rt_arr,peak_only_eic_int,color='red',interpolate=True)
        # ax.fill_between(pl.arange(rtMinIndex,rtMaxIndex),peak_only_eic_int,color='red',interpolate=True)
        ax.fill_between(rt[rtMinIndex:rtMaxIndex], peak_only_eic_int, color='red', interpolate=True)

        ax.annotate("C/A=%.2f" % (curPeak.coefOverArea), xy=(0.130, .870), xycoords='figure fraction', fontsize=8)
        ax.annotate("AsymFit=%.2f" % (curPeak.asymGausFit), xy=(0.130, .80), xycoords='figure fraction', fontsize=8)
        # ax.annotate("RT=%.2f"%(curPeak.rt),xy=(0.130,.73),xycoords='figure fraction',fontsize=8)
        yloc = 0.73
        for simIdx in range(len(curPeak.similarityAdjMZ)):
            ax.annotate("Sim=%.2f" % (curPeak.similarityAdjMZ[simIdx]), xy=(0.130, yloc), xycoords='figure fraction',
                        fontsize=8)
            yloc -= 0.07

        pl.savefig(plotDirStr + "/" + str(k) + "_mz" + str(cur_peak_mz)+"_rt"+ str(cur_peak_rt) + ".pdf", dpi=300)
        pl.close(fig)

    return listOfAllPeakMZAndRT, countRepeats

def plot_all_the_peaks(result_dir_str,inital_peaks_objects,allPeakObjects,int_matrix,rt):
    """
    Takes a list of peak objects and plots all of them including their isotopes. Isotopes in separate folder.
    :param allPeakObjects: List of peak objects
    :param int_matrix: matrix containing all data. Rows -> mz*10000, Cols -> RT in scan number, Values -> intensities
    :param rt: Retention time arr for plots 
    :return: nothing
    """



    plotDirStr = result_dir_str+"PeakPlots"
    isoPlotDirStr = result_dir_str+"IsoPeakPlots"
    initialPeaksPlotDirStr = result_dir_str+"ParamDetectionPeakPlots"

    if os.path.exists(initialPeaksPlotDirStr):
        shutil.rmtree(initialPeaksPlotDirStr)
    os.mkdir(initialPeaksPlotDirStr)
    if os.path.exists(plotDirStr):
        shutil.rmtree(plotDirStr)
    os.mkdir(plotDirStr)
    if os.path.exists(isoPlotDirStr):
        shutil.rmtree(isoPlotDirStr)
    os.mkdir(isoPlotDirStr)

    listOfAllPeakMZAndRT = []
    countRepeats = 0

    ############################
    ##### plot all peaks #######
    ############################
    listOfAllPeakMZAndRT, countRepeats = plot_these_peaks(allPeakObjects,int_matrix,rt,plotDirStr,listOfAllPeakMZAndRT,countRepeats)

    ############################
    ##### plot isotopes ########
    ############################
    for k in range(len(allPeakObjects)):
        curPeak = allPeakObjects[k]
        if curPeak.hasIsotope:
            listOfAllPeakMZAndRT, countRepeats = plot_these_peaks(curPeak.isotopeList,int_matrix,rt,isoPlotDirStr,listOfAllPeakMZAndRT,countRepeats)


    ############################
    ##### plot initial peaks ###
    ############################
    # These are separate. They should be contained in the allPeakObjects list but we
    # also want to plot them separately so it is easy to see what was used to find
    # the parameters
    plot_these_peaks(inital_peaks_objects,int_matrix,rt,initialPeaksPlotDirStr,[],0)


    print "Number of repeated peaks avoided while doing plots:"
    print "-----> " + str(countRepeats)