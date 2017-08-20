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

# import matplotlib
# matplotlib.use('Agg')


import sys
import os
import shutil
import pylab as pl
import scipy.integrate as integrate
from scipy.sparse import dok_matrix
import argparse

import generalcurvetools.curve_tools as cts
from plotting import peak_plotting
from datamodel import data_point
from datamodel import peak
from datamodel import result
from writing import write_csv
from writing import write_all_peak_data
from datamodel import hard_parameter
from peakdetectioncwt import peakdetector
from easyIOmassspec import easyio
from datamodel import efficient_find_next_max


PLOT_ALL_PEAKS = True

VERBOSE = True

USE_HARD_CODED_DETECTION_PARAMETERS = False
USE_ISOTOPE_PARAMETERS_FOR_ALL = False
USE_SMALL_TEST_WINDOW = False

# isotopes will be shown for the tmp_mz_of_peak value set
ONLY_VISUALIZE_ISOTOPES_FOR_NORMAL_DETECTED_PEAK = False

if USE_SMALL_TEST_WINDOW:
    tmp_mz_of_peak = 212.006378174  # 20.25 21.75 rt
    # RT_MIN = 17.25
    # RT_MAX = 21.75
    # MZ_MIN = tmp_mz_of_peak - 10
    # MZ_MAX = tmp_mz_of_peak + 10
    # RT_MIN = 19.0
    # RT_MAX = 21.0
    # MZ_MIN = 115
    # MZ_MAX = 117

    RT_MIN = 21.0
    RT_MAX = 24.0
    MZ_MIN = 114.0
    MZ_MAX = 119.0

##########################################################################
########### Important numbers that need to be set ########################
##########################################################################
HP = hard_parameter.RequiredParameters()
if USE_ISOTOPE_PARAMETERS_FOR_ALL:
    HP.use_isotope_parameters_for_all()

##########################################################################
############ Done setting important numbers ##############################
##########################################################################

if ONLY_VISUALIZE_ISOTOPES_FOR_NORMAL_DETECTED_PEAK:
    investigate_specific_mz = [-1]
else:
    investigate_specific_mz = [-1]
    # investigate_specific_mz = [tmp_mz_of_peak]

RESULT_LIST = []
ISOTOPE_RESULT_PEAK_LIST = []


def store_result(result):
    """ This is for multiprocessing to store the results in the global RESULT_LIST

    :param result: result to append to list
    :type result: result objects (see data_model)
    """

    RESULT_LIST.append(result)


def isotope_store_result(result):
    """ So multiprocessing can store the results in global ISOTOPE_RESULT_PEAK_LIST

    :param result: list of results to append to list
    :type result: list of peak objects (see data_model)
    """

    for i in result:
        ISOTOPE_RESULT_PEAK_LIST.append(i)


def peak_similarity_test(int_matrix, rightB, leftB, estimate_of_ms_fwhm, a, visualize):
    """
    For now lest use the normalized area similarity I came up with. It will be more sensitive.
    Starting with the EIC sliced from the maxima, move up and down in mz and compare the profiles of
    the slices with the maxima slice. 0-> bad similarity 1-> good similarity

    :param int_matrix: matrix of all intensities
    :param rightB: index of right bound of the peak.
    :param leftB: index of left bound of the peak.
    :param estimate_of_ms_fwhm: estimate of mass spec full width half max in number of points.
    :param a: row (mz) index of peak
    :param visualize: Boolean-> do you want to see what is going on?
    :return: good peak boolean, mz upper bound (determined from similarity measures), mz lower bound,
            list of all similarity values.
    """

    # Parameters for exponential transform of similarity value.
    # maps [0,1] to [0,1] but should separate the good values from
    # the bad a little more clearly.
    param_c = 0.2
    param_b = pl.exp(-1.0 / param_c)
    param_a = 1.0 / (1.0 - param_b)

    # EIC of highest point. leftB and rightB are the boundaries of the peak found
    # from the CWT and boundary correction
    just_peak_eic = pl.squeeze(pl.asarray(int_matrix[a, leftB:rightB].todense()))
    x = pl.arange(0, len(just_peak_eic))
    # trapezoidal integration.
    area = integrate.trapz(just_peak_eic, x=x)

    normalized_eic = just_peak_eic / area

    if visualize:
        pl.plot(normalized_eic)

    # initialize some values before while loop.
    cur_similarity = 1.0
    # count the number of similar EICs for increasing m/z
    cur_inc = 0
    # if we find zero area EICs we skip them so we need this to correct for that
    subtract_from_cur_inc = 0

    # keep track of the similarity values that are measured so that
    # We can include the values in the output plots.
    sim_vals = []

    while cur_similarity > HP.get_peak_similarity_threshold():
        cur_inc += 1
        # Make sure we don't go too far out. 2 times the full width half max of a
        # mass spec peak should be plenty.
        if ((a + cur_inc) >= int_matrix.shape[0]) or (cur_inc >= 2 * estimate_of_ms_fwhm):
            break
        cur_eic = pl.squeeze(pl.asarray(int_matrix[a + cur_inc, leftB:rightB].todense()))
        # because of imperfections in building the the intensity matrix (int_matrix) there can be
        # gaps between real slices of data filled with zeros. This is because somewhere
        # there are slightly different mz sampling intervals for different m/z that mess up the
        # assumed perfect grid of the matrix.
        area2 = integrate.trapz(cur_eic, x=x)
        # epsilon is just a "very small" number.
        if area2 < HP.get_epsilon():
            # If we are in here we are skipping one of the empty slices mentioned above.
            subtract_from_cur_inc += 1
            continue
        normed_cur_eic = cur_eic / area2
        diff_area = integrate.trapz(abs(normalized_eic - normed_cur_eic), x=x)

        if VERBOSE: print ("diff_area: " + str(diff_area))

        # cur_similarity = 1-diff_area
        cur_similarity = (pl.exp(-diff_area / param_c) - param_b) * param_a

        if visualize and (cur_similarity > HP.get_peak_similarity_threshold()):
            print ("cur_similarity: " + str(cur_similarity))
            pl.plot(normed_cur_eic)
        if cur_similarity > HP.get_peak_similarity_threshold():
            sim_vals.append(cur_similarity)

    cur_inc -= 1
    cur_inc -= subtract_from_cur_inc

    if VERBOSE or visualize:
        print("cur_inc: " + str(cur_inc))
        print("subtract_from_cur_inc: " + str(subtract_from_cur_inc))

    # count the number of similar EICs for decreasing m/z. For details
    # see comments above.
    cur_inc2 = 0
    cur_similarity = 1.0
    subtract_from_cur_inc = 0
    while cur_similarity > HP.get_peak_similarity_threshold():
        cur_inc2 += 1
        if ((a - cur_inc2) < 0) or (cur_inc >= 2 * estimate_of_ms_fwhm):
            break
        cur_eic = pl.squeeze(pl.asarray(int_matrix[a - cur_inc2, leftB:rightB].todense()))
        area2 = integrate.trapz(cur_eic, x=x)
        if area2 < HP.get_epsilon():
            subtract_from_cur_inc += 1
            continue
        normed_cur_eic = cur_eic / area2
        diff_area = integrate.trapz(abs(normalized_eic - normed_cur_eic), x=x)

        if VERBOSE: print("diff_area: " + str(diff_area))

        # cur_similarity = 1-diff_area
        cur_similarity = (pl.exp(-diff_area / param_c) - param_b) * param_a

        if visualize and (cur_similarity > HP.get_peak_similarity_threshold()):
            print("cur_similarity: " + str(cur_similarity))
            pl.plot(normed_cur_eic)
        if cur_similarity > HP.get_peak_similarity_threshold():
            sim_vals.append(cur_similarity)

    cur_inc2 -= 1
    cur_inc2 -= subtract_from_cur_inc

    if VERBOSE or visualize:
        print ("cur_inc2: " + str(cur_inc2))

    if visualize:
        pl.title("similar normalized curves")
        pl.show()

    # HP.get_peak_num_sim_either_side_threshold() = int(round(estimate_of_ms_fwhm/2.0))

    if (cur_inc >= HP.get_peak_num_sim_either_side_threshold()) and (
        cur_inc2 >= HP.get_peak_num_sim_either_side_threshold()) and (
        (cur_inc + cur_inc2) >= HP.get_peak_num_sim_total_threshold()):
        good_peak_bool = True
    else:
        good_peak_bool = False

    if (a + cur_inc + 1) >= int_matrix.shape[0]:
        to_return_up = int_matrix.shape[0] - 1
    else:
        to_return_up = a + cur_inc + 1
    if (a - cur_inc2 - 1) < 0:
        to_return_down = 0
    else:
        to_return_down = a - cur_inc2 - 1

    if VERBOSE or visualize:
        print "good_peak_bool: " + str(good_peak_bool)

    return good_peak_bool, to_return_up, to_return_down, sim_vals


def remove_bad_peak_info_best_we_can_sparse(to_mod_int_matrix,
                                            estimate_of_ms_fwhm,
                                            eic_left_bound,
                                            eic_right_bound,
                                            a,
                                            visualize,
                                            efficient_next_max):
    """
    This function will remove, from to_mod_int_matrix, the data that should be bad. It uses the
    full width at half max of the mass spec data (estimated value) to select a strip of data that will
    be set to zero. Everything in the retention time range of the EIC is removed. There is a danger that
    the edge of the EIC will be a good peak but its ok because when a new EIC is made around a maxima it
    uses the int_matrix from which no data has been removed.

    :param to_mod_int_matrix: (sparce matrix) The matrix of intensities that will have values set to zero if they are "bad" (not part of peak)
    :param estimate_of_ms_fwhm: (int) Estimated value of the full width at half max of a mass spec in number of data points. scan peak.
    :param eic_left_bound: (int) Scan number left boundary of EIC investigated.
    :param eic_right_bound: (int) Scan number right boundary of peak investigated.
    :param a: (int) Row (mz) index of peak (even though it is bad peak)
    :param visualize: (boolean) Visualize whats happening?
    :return: return the modified intensity matrix
    """

    rm_mz_down = a - estimate_of_ms_fwhm
    rm_mz_up = a + estimate_of_ms_fwhm
    if rm_mz_down < 0:
        rm_mz_down = 0
    if rm_mz_up > to_mod_int_matrix.shape[0]:
        rm_mz_up = to_mod_int_matrix.shape[0]

    mz_low_bound = rm_mz_down - 15
    mz_high_bound = rm_mz_down + 15
    if mz_low_bound < 0:
        mz_low_bound = 0
    if mz_high_bound > to_mod_int_matrix.shape[0]:
        mz_high_bound = to_mod_int_matrix.shape[0]

    if visualize:
        pl.pcolor(pl.asarray(
            to_mod_int_matrix[mz_low_bound:mz_high_bound, max(eic_left_bound - 10, 0):eic_right_bound + 10].todense()))
        pl.colorbar()
        pl.title('before removal')
        pl.show()

    to_mod_int_matrix[rm_mz_down:rm_mz_up, eic_left_bound:eic_right_bound + 1] = 0
    efficient_next_max.done_with_rows_cols(rm_mz_down, rm_mz_up, eic_left_bound, eic_right_bound + 1)

    if visualize:
        pl.pcolor(pl.asarray(
            to_mod_int_matrix[mz_low_bound:mz_high_bound, max(eic_left_bound - 10, 0):eic_right_bound + 10].todense()))
        pl.colorbar()
        pl.title('after removal')
        pl.show()


def get_peakwidth_and_coef_over_area(int_matrix,
                                     all_peak_positions,
                                     all_right_bounds,
                                     all_left_bounds,
                                     eic_left_bound,
                                     a,
                                     b,
                                     visualize,
                                     all_coef_over_area,
                                     all_wavelet_similarity_vals):
    """
    Select the correct peak information to return based off of the maxima being investigated.
    Maxima under investigation is at int_matrix[a,b].
    :param int_matrix: Matrix of intensity values.
    :param all_peak_positions: (list int) Indices of the peak positions.
    :param all_right_bounds: (list int) Indices of peak right most boundaries.
    :param all_left_bounds: (list int) Indices of peak left most boundaries.
    :param eic_left_bound: (int) Index of EIC left bound.
    :param a: (int) Row (mz) of maxima investigating.
    :param b: (int) Col (rt) of maxima investigating.
    :param visualize: (boolean) See the peaks?
    :param all_coef_over_area: (list float) all coefficient over area values for each peak.
    :param all_wavelet_similarity_vals: (list float) all wavelet similarity values.
    :return: (int) scan number (index) of right side of final peak,
                (int) scan number (index) of left side of final peak,
                coef/area of final peak,
                wavelet similarity of final peak.
    """

    scan_num_eic_right = 0
    scan_num_eic_left = 0
    cur_coef_over_area = 0
    cur_wav_peak_sim = 0

    for i, j in enumerate(all_peak_positions):

        # make sure you are getting the peak that contains the relevant maxima
        eic_max_index = j
        scan_num_index_max = j + eic_left_bound
        if VERBOSE:
            print ("------ THESE NEED TO BE THE SAME FOR A GOOD PEAK -------")
            print ("scan_num_index_max: " + str(scan_num_index_max))
            print ("b: " + str(b))
            print ("--------------------------------------------------------")

        # b is the index of the maximum intensity point of the data file that is currently under
        # investigation. Want to get the peak that is found (if there is one) at this point.
        if scan_num_index_max != b:
            continue

        # check similarity between adjacent mz EICs, both to determine the boundary in m/z of
        # removal and to make sure this is a good peak. Lets use cts.estimate_fwhm_ms to determine
        # an appropriate number of similar EICs for it to be a good peak.
        # Lets say on either side if must have FWHM/2 + 1 (in units of num. intervals) similar
        # EICs to be good peak

        # first just grab the intensities of the actual peak
        scan_num_eic_left = all_left_bounds[i] + eic_left_bound
        scan_num_eic_right = all_right_bounds[i] + eic_left_bound
        if VERBOSE:
            print ("scan_num_eic_left : " + str(scan_num_eic_left))
            print ("scan_num_eic_right : " + str(scan_num_eic_right))
        eic_peak = int_matrix[a, scan_num_eic_left: scan_num_eic_right]

        if visualize:
            pl.plot(pl.squeeze(pl.asarray(eic_peak.todense())))
            pl.title("does this look like the right peak???")
            pl.show()

        cur_coef_over_area = all_coef_over_area[i]
        cur_wav_peak_sim = all_wavelet_similarity_vals[i]

    return scan_num_eic_right, scan_num_eic_left, cur_coef_over_area, cur_wav_peak_sim


def plot_color_map_peak_and_box_bounds(int_matrix, rightB, leftB, mz_up_bound, mz_low_bound):
    """
    Plot a color map of the peak a a box around it showing it's detected boundaries.

    :param int_matrix: Matrix of intensities.
    :param rightB: (int) index of right bound.
    :param leftB:  (int) index of left bound.
    :param mz_up_bound: (int) Index of top bound.
    :param mz_low_bound: (int) Index of low bound.
    :return: nothing
    """
    mzl = mz_low_bound - 15
    mzr = mz_up_bound + 15
    if mzl < 0:
        mzl = 0
    if mzr > int_matrix.shape[0]:
        mzr = int_matrix.shape[0]
    rtl = leftB - 10
    rtr = rightB + 10
    if rtl < 0:
        rtl = 0
    if rtr > int_matrix.shape[1]:
        rtl = int_matrix.shape[1]
    print "pl.shape(int_matrix): " + str(pl.shape(int_matrix))

    shiftMZ = mz_low_bound - mzl
    shiftRT = leftB - rtl

    pl.pcolor(pl.asarray(int_matrix[mzl:mzr, rtl:rtr].todense()))
    pl.plot([shiftRT, shiftRT + rightB - leftB + 1], [shiftMZ, shiftMZ], c='r')
    pl.plot([shiftRT, shiftRT + rightB - leftB + 1],
            [shiftMZ + mz_up_bound - mz_low_bound + 1, shiftMZ + mz_up_bound - mz_low_bound + 1], c='r')
    pl.plot([shiftRT, shiftRT], [shiftMZ, shiftMZ + mz_up_bound - mz_low_bound + 1], c='r')
    pl.plot([shiftRT + rightB - leftB + 1, shiftRT + rightB - leftB + 1],
            [shiftMZ, shiftMZ + mz_up_bound - mz_low_bound + 1], c='r')
    pl.title("color plot of peak and detected 2D bounds")
    pl.show()


def get_peak_list(peakDetector,
                  int_matrix,
                  to_mod_int_matrix,
                  eic_num_scans_leftright,
                  visualize_post_initial_peaks,
                  estimate_of_ms_fwhm,
                  min_intensity_thresh,
                  rt,
                  unique_mz_list,
                  efficient_next_max,
                  mz_index_low=0,
                  rt_index_low=0
                  ):
    """
    This takes a small section of the entire data set and looks for peaks within it.
    Starting with the most intense data point it take a slice in RT and performs a wavelet
    transform to look for peaks. If it finds a peak at the maxima under investigation it looks at
    adjacent slices in the m/z domain to determine if the peak is good by comparing the similarity
    of those curves with the curve at the maxima. If not good peaks (at all) are found it sets all that
    data to zero in to_mod_int_matrix before performing a search for the next local maxima. If good peaks are
    found that are not at the most intense location, all the data EXCEPT the data corresponding to those good
    peaks is set to zero. A good peak is added to the result list before that data is set to zero. int_matrix
    is never modified so all EICs contain accurate description of data.

    :param peakDetector: Peak detection object configured for finding peaks in EICs.
    :param int_matrix: Matrix of all intensity values (row->m/z, col->rt).
    :param to_mod_int_matrix: Modified int_matrix. Detected peaks and bad peaks have been set to zero.
    :param eic_num_scans_leftright: (int) How many scans in RT domain to build EIC around maxima.
    :param visualize_post_initial_peaks: (boolean) see the process?
    :param estimate_of_ms_fwhm: (int) Estimated value of the full width at half max of a mass spec in number of data points. scan peak.
    :param min_intensity_thresh: (float) Peaks below this intensity value are not considered.
    :param rt: (list float) list of retention times.
    :param unique_mz_list: (list int (floats rounded to 4th decimal place then times 10000 and converted to int))
        list of all unique mz values as integers
    :param mz_low_edge: (int) Because the matrices are only portions of the data this variable tells you where the
        lower m/z (row) edge of that data begins in the full data set.
    :param rt_low_edge: (int) Same as above but reference for RT (col).
    :return: returns result object with all the detected peaks.
    """

    result_obj = result.Result()

    a, b = efficient_next_max.find_max()
    cur_most_max = int_matrix[a, b]
    while min_intensity_thresh < cur_most_max:

        # If both are negative 1 it means there are no non-removed points left.
        if (a == -1) and (b == -1):
            break

        for alpha in investigate_specific_mz:
            if ((unique_mz_list[a] / 10000.0) < (alpha + 0.005)) and ((unique_mz_list[a] / 10000.0) > (alpha - 0.005)):
                global VERBOSE
                VERBOSE = True
                print "\n**************** investigate this mz *********************"
                print "**************** " + str(alpha) + "*********************"
                print "Current actual mz: " + str(unique_mz_list[a] / 10000.0)
                print "Current rt: " + str(rt[b])
                print "\n**********************************************************"
                peakDetector.setVerbose(True)
                peakDetector.setVisualize(True)
                visualize_post_initial_peaks = True

        if VERBOSE:
            print "Current actual mz: " + str(unique_mz_list[a] / 10000.0)
            print "Current rt: " + str(rt[b])

        # get the "EIC" associated with this maxima. Use the original intensity matrix so the eic does
        # not contain a discontinuity from the removal of data
        eic_left_bound = b - eic_num_scans_leftright
        eic_right_bound = b + eic_num_scans_leftright
        if eic_left_bound < 0:
            eic_left_bound = 0
        if eic_right_bound >= int_matrix.shape[1]:
            eic_right_bound = int_matrix.shape[1] - 1
        eic = int_matrix[a, eic_left_bound:eic_right_bound]
        eic = pl.squeeze(pl.asarray(eic.todense()))

        # Now do peak picking. It is possible that several peaks will be found so you will have to check
        # each for the one with boundaries that contain the index of the maximum found in the
        # "pl.unravel_index..." line.
        peakDetector.setSignal(eic)
        peakDetector.setX(rt[eic_left_bound:eic_right_bound])
        all_left_bounds, all_right_bounds, all_peak_positions = peakDetector.findPeaks()
        all_coef_over_area_vals = peakDetector.getallCoefOverAreaVals()
        all_wavelet_similarity_vals = peakDetector.getallWaveletSimilarityVals()

        if VERBOSE:
            print ("len(all_right_bounds): " + str(len(all_right_bounds)))
            print ("eic_left_bound: " + str(eic_left_bound))
            print ("eic_right_bound: " + str(eic_right_bound))
            print ("len(rt): " + str(len(rt)))
            print ("a: " + str(a))
            print ("b: " + str(b))
            print "len(int_matrix[a,:]): " + str(int_matrix.shape[1])
            print "all peaks information:"
            for beta in range(len(all_coef_over_area_vals)):
                print "     all_coef_over_area_vals[beta]: " + str(all_coef_over_area_vals[beta])
                print "     all_left_bounds[beta]: " + str(all_left_bounds[beta])
                print "     all_right_bounds[beta]: " + str(all_right_bounds[beta])
                print "     all_peak_positions[beta]: " + str(all_peak_positions[beta])
                print "     all_wavelet_similarity_vals[beta]: " + str(all_wavelet_similarity_vals[beta])

        # If there are no good peaks then you can not use similarity between adjacent mz EICs to
        # determine boundary in m/z (for point removal). We will have to use the cts.estimate_fwhm_ms
        if len(all_right_bounds) == 0:
            if VERBOSE:
                print (
                "remove_bad_peak_info_best_we_can\n either because b not in (all_peak_positions+eic_left_bound)\n or no peak detected using original cwt peak detector.")
                if len(all_right_bounds) == 0: print ("reason -----> len(all_right_bounds)==0)")

            remove_bad_peak_info_best_we_can_sparse(to_mod_int_matrix, estimate_of_ms_fwhm, eic_left_bound,
                                                    eic_right_bound, a, visualize_post_initial_peaks,
                                                    efficient_next_max)

        # if we pass the above if and are here it means that there is a good peak somewhere in
        # here -> just not the one associated with the current maxima under investigation. Be
        # careful not the remove the good peak(s)
        elif b not in (pl.array(all_peak_positions) + eic_left_bound):
            if VERBOSE:
                print ("IN -----> b not in (pl.array(all_peak_positions)+eic_left_bound)")
                print "     pl.array(all_peak_positions) = " + str(pl.array(all_peak_positions))
                print "     pl.array(all_peak_positions)+eic_left_bound = " + str(
                    pl.array(all_peak_positions) + eic_left_bound)
                print "     eic_left_bound = " + str(eic_left_bound)
                print "     b = " + str(b)
            # set everything to zero and then set the peaks back to their original intensities
            remove_bad_peak_info_best_we_can_sparse(to_mod_int_matrix, estimate_of_ms_fwhm, eic_left_bound,
                                                    eic_right_bound, a, visualize_post_initial_peaks,
                                                    efficient_next_max)
            for alpha in range(len(all_right_bounds)):
                # don't add anything that has already been removed
                if to_mod_int_matrix[a, eic_left_bound + all_peak_positions[alpha]] < HP.get_epsilon():
                    continue

                put_good_peak_back_mz_low = a - 2 * estimate_of_ms_fwhm
                put_good_peak_back_mz_high = a + 2 * estimate_of_ms_fwhm
                if put_good_peak_back_mz_low < 0:
                    put_good_peak_back_mz_low = 0
                if put_good_peak_back_mz_high > to_mod_int_matrix.shape[0]:
                    put_good_peak_back_mz_high = to_mod_int_matrix.shape[0]

                to_mod_int_matrix[put_good_peak_back_mz_low:put_good_peak_back_mz_high,
                eic_left_bound + all_left_bounds[alpha]:eic_left_bound + all_right_bounds[alpha]] \
                    = pl.copy(int_matrix[put_good_peak_back_mz_low:put_good_peak_back_mz_high,
                              eic_left_bound + all_left_bounds[alpha]:eic_left_bound + all_right_bounds[alpha]])

                efficient_find_next_max.put_back_rows_cols(put_good_peak_back_mz_low, put_good_peak_back_mz_high,
                                                           eic_left_bound + all_left_bounds[alpha],
                                                           eic_left_bound + all_right_bounds[alpha])

        else:
            # We choose the peak whose apex scan number coincides with b and return its characteristics
            # The bounds that are returned from this are the indices in the int_matrix
            rightB, leftB, cur_coef_over_area, cur_wav_peak_sim = get_peakwidth_and_coef_over_area(int_matrix,
                                                                                                   all_peak_positions,
                                                                                                   all_right_bounds,
                                                                                                   all_left_bounds,
                                                                                                   eic_left_bound,
                                                                                                   a,
                                                                                                   b,
                                                                                                   visualize_post_initial_peaks,
                                                                                                   all_coef_over_area_vals,
                                                                                                   all_wavelet_similarity_vals)

            # Look at similarity for different mz slices -> tells you if its a good peak and what
            # the boundaries for removal in m/z are.
            # returns if it is a good peak, and the boundaries
            good_peak_bool, mz_up_bound, mz_low_bound, sim_vals = peak_similarity_test(int_matrix, rightB, leftB,
                                                                                       estimate_of_ms_fwhm, a,
                                                                                       visualize_post_initial_peaks)

            fit_parameters, cur_asym_similarity = cts.asymmetric_gaussian_fit(b - leftB,
                                                                              pl.squeeze(pl.asarray(int_matrix[a,
                                                                                                    leftB:rightB].todense())),
                                                                              visualize_post_initial_peaks)
            # remove the peak and appropriate mz
            to_mod_int_matrix[mz_low_bound:mz_up_bound + 1, leftB:rightB + 1] = 0.0
            efficient_next_max.done_with_rows_cols(mz_low_bound, mz_up_bound + 1, leftB, rightB + 1)

            if visualize_post_initial_peaks:
                plot_color_map_peak_and_box_bounds(int_matrix, rightB, leftB, mz_up_bound, mz_low_bound)

            # if the peak is too close to the edge don't consider it.
            if b < ((rightB - leftB) / 2):
                if VERBOSE:
                    print "Wait -> not a good peak. Peak is less than PW/2 away from boundry."
                good_peak_bool = False
            if cur_asym_similarity < HP.get_peak_asymmetric_fit_similarity_threshold():
                good_peak_bool = False
                if VERBOSE:
                    print "PEAK NOT FOUND REASON:"
                    print "     cur_asym_similarity: " + str(cur_asym_similarity)

            if good_peak_bool:

                peak_mz = unique_mz_list[a] / 10000.0
                peak_rt = rt[b]

                if VERBOSE:
                    print ("peak_mz: " + str(peak_mz))
                    print ("peak_rt: " + str(peak_rt))
                    print ("peak left bound:  " + str(rt[leftB]))
                    print ("peak right bound: " + str(rt[rightB]))
                    print ("peak mz high bound:  " + str(unique_mz_list[mz_up_bound] / 10000.0))
                    print ("peak mz low bound: " + str(unique_mz_list[mz_low_bound] / 10000.0))
                    print "cur_asym_similarity: " + str(cur_asym_similarity)

                result_obj.addRT(peak_rt)
                result_obj.addMZ(peak_mz)
                result_obj.addHeight(int_matrix[a, b])
                result_obj.addMaxMZ(unique_mz_list[mz_up_bound] / 10000.0)
                result_obj.addMinMZ(unique_mz_list[mz_low_bound] / 10000.0)
                result_obj.addMinRT(rt[leftB])
                result_obj.addMaxRT(rt[rightB])
                result_obj.addMZIndex(mz_index_low + a)
                result_obj.addRTIndex(rt_index_low + b)
                result_obj.addMZMinIndex(mz_index_low + mz_low_bound)
                result_obj.addMZMaxIndex(mz_index_low + mz_up_bound)
                result_obj.addRTMinIndex(rt_index_low + leftB)
                result_obj.addRTMaxIndex(rt_index_low + rightB)
                result_obj.addCoefOverArea(cur_coef_over_area)
                result_obj.addAsymGausFitVal(cur_asym_similarity)
                result_obj.addSimilarityAdjMZVals(sim_vals)

        for alpha in investigate_specific_mz:
            if ((unique_mz_list[a] / 10000.0) < (alpha + 0.005)) and ((unique_mz_list[a] / 10000.0) > (alpha - 0.005)):
                VERBOSE = False
                visualize_post_initial_peaks = False
                peakDetector.setVisualize(False)
                peakDetector.setVerbose(False)

        # indices of maximum. Use the modified matrix so that vales can be set to zero to find the next
        # maxima.
        a, b = efficient_next_max.find_max()
        if (a == -1) and (b == -1):
            break
        cur_most_max = int_matrix[a, b]

    return result_obj


def find_isotope_peak(cur_peak_object,
                      int_matrix,
                      isotope_peak_detector,
                      eic_num_scans_leftright,
                      visualize_post_initial_peaks,
                      estimate_of_ms_fwhm,
                      rt,
                      unique_mz_list
                      ):
    """
    Given a peak (presumably monoisotopic but not necessarily) find the isotope(s) of this peak.
    Searches in a narrow m/z strips near where the isotopes are expected to be. If one isotope is found
    it searches for the next. If an isotope is not found it stops.
    Uses relaxed parameters for search.
    Uses a linear scoring mechanism if multiple candidates are found. Compares intensity, m/z distance, rt distance,
    similarity to input peak to choose the best possible candidate.

    :param cur_peak_object: (Peak object) Input peak.
    :param int_matrix: Matrix of all intensities.
    :param isotope_peak_detector: (peak detection object) Configured for detection of isotopes.
    :param eic_num_scans_leftright: (int) How many scans used to make EIC in each direction.
    :param visualize_post_initial_peaks: (boolean) See process?
    :param estimate_of_ms_fwhm: (int) Estimate of mass spec scan full width at half max in number of points.
    :param rt: (list) List of retention times.
    :param unique_mz_list: (list int (floats rounded to 4th decimal place then times 10000 and converted to int))
        list of all unique mz values as integers
    :return: Returns cur_peak_object but with isotope list full of peaks that are isotopes detected.
    """

    global VERBOSE

    # Threshold (in number of scans) of how close the peaks need to be between the isotopes to be
    # consider good
    rt_scan_difference_threshold = 4
    # Threshold between the expected MZ index of the isotope and the index that is found difference
    mz_index_difference_threshold = 2 * estimate_of_ms_fwhm

    plusIso = 1.0000

    count = 0
    found_isotope = True

    if ONLY_VISUALIZE_ISOTOPES_FOR_NORMAL_DETECTED_PEAK:
        visualize_post_initial_peaks = True
        isotope_peak_detector.setVerbose(True)
        isotope_peak_detector.setVisualize(True)
        VERBOSE = True

    while found_isotope:
        found_isotope = False

        count += 1

        addToMZ = count * plusIso

        # Get the data in the narrow rectangle where the isotope should be
        # I think you have to use the it_matrix and not the to_mod_int_matrix because you are interested
        # in peaks with relaxed parameters so their intensities might have been set to zero in the
        # to_int_mod_matrix.

        # get the index of the best guess for the isotopes mz location
        bestMZGuess = cur_peak_object.mz + addToMZ
        scaledIntBestMZGuess = int(round(bestMZGuess * 10000))
        best_mz_idex_guess = unique_mz_list.index(min(unique_mz_list, key=lambda x: abs(x - scaledIntBestMZGuess)))
        actual_mz_at_best_guess_index = unique_mz_list[best_mz_idex_guess]

        if abs(actual_mz_at_best_guess_index - scaledIntBestMZGuess) > 100:
            break

        # best_mz_idex_guess=unique_mz_list.index(scaledIntBestMZGuess)
        if VERBOSE:
            print "New Isotope -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0"
            print "bestMZGuess " + str(bestMZGuess)
            print "scaledIntBestMZGuess " + str(scaledIntBestMZGuess)
            print "best_mz_idex_guess " + str(best_mz_idex_guess)
            print "whats at best_mz_idex_guess " + str(unique_mz_list[best_mz_idex_guess])

        mz_index_low = best_mz_idex_guess - 4 * estimate_of_ms_fwhm
        mz_index_high = best_mz_idex_guess + 4 * estimate_of_ms_fwhm

        if mz_index_low < 0: mz_index_low = 0
        if mz_index_high > len(unique_mz_list): mz_index_high = len(unique_mz_list)

        rt_index = cur_peak_object.rtIndex
        if VERBOSE:
            print "rt_index " + str(rt_index)

        patch_left = rt_index - eic_num_scans_leftright
        patch_right = rt_index + eic_num_scans_leftright
        if patch_left < 0:
            patch_left = 0
        if patch_right >= int_matrix.shape[1]:
            patch_right = int_matrix.shape[1] - 1

        if VERBOSE:
            print "patch_left " + str(patch_left)
            print "patch_right " + str(patch_right)

        re_fill_to_mod_int_matrix = int_matrix[mz_index_low:mz_index_high, patch_left:patch_right].copy()

        # make efficient max point finder object. To do this we need to make a list of the
        # data points in the small portion of the matrix we are investigating.
        cur_data_point_list = []
        # The next 3 lines are to get the indices of the non-zero
        # points in the sparse matrix.
        nonzero_indecies = re_fill_to_mod_int_matrix.nonzero()
        nonzero_row_indx = nonzero_indecies[0]
        nonzero_col_indx = nonzero_indecies[1]
        for i in range(len(nonzero_row_indx)):
            cur_mz_index = nonzero_row_indx[i]
            cur_peak_mz = unique_mz_list[cur_mz_index] / 10000.0
            cur_rt_index = nonzero_col_indx[i]
            cur_data_point = data_point.DataPoint(cur_rt_index,
                                                  cur_mz_index,
                                                  cur_peak_mz,
                                                  re_fill_to_mod_int_matrix[cur_mz_index, cur_rt_index])
            cur_data_point_list.append(cur_data_point)

        cur_efficient_next_max = efficient_find_next_max.EfficientNextMax(cur_data_point_list)

        resultO = get_peak_list(isotope_peak_detector,
                                int_matrix[mz_index_low:mz_index_high, patch_left:patch_right],
                                re_fill_to_mod_int_matrix,
                                eic_num_scans_leftright,
                                visualize_post_initial_peaks,
                                estimate_of_ms_fwhm,
                                0,
                                rt[patch_left:patch_right],
                                unique_mz_list[mz_index_low:mz_index_high],
                                cur_efficient_next_max,
                                mz_index_low=mz_index_low,
                                rt_index_low=patch_left
                                )

        all_peak_objects = []
        for i in range(len(resultO.final_peak_list_rt_vals)):
            cur_peak = peak.Peak()
            cur_peak.setRTIndex(resultO.rt_index[i])
            cur_peak.setMZIndex(resultO.mz_index[i])
            cur_peak.setRT(resultO.final_peak_list_rt_vals[i])
            cur_peak.setMZ(resultO.final_peak_list_mz_vals[i])
            cur_peak.setHeight(resultO.final_peak_list_peak_height_vals[i])
            cur_peak.setMZmin(resultO.final_peak_list_mz_min_vals[i])
            cur_peak.setMZmax(resultO.final_peak_list_mz_max_vals[i])
            cur_peak.setRTMin(resultO.final_peak_list_peak_rt_start_vals[i])
            cur_peak.setRTMax(resultO.final_peak_list_peak_rt_end_vals[i])

            cur_peak.setMZMinIndex(resultO.mz_min_index[i])
            cur_peak.setMZMaxIndex(resultO.mz_max_index[i])
            cur_peak.setRTMinIndex(resultO.rt_min_index[i])
            cur_peak.setRTMaxIndex(resultO.rt_max_index[i])

            cur_peak.setCoefOverArea(resultO.coef_over_areas[i])
            cur_peak.setAsymGausFit(resultO.asym_gaus_fit_vals[i])

            all_peak_objects.append(cur_peak)

        if len(all_peak_objects) == 0:
            # print "len(all_peak_objects) "+str(len(all_peak_objects))
            found_isotope = False

        else:
            # figure out which is the peak is the right one if there are multiple.
            # If there is only one then make sure it satisfies some conditions before recursivly calling
            # this function again.
            # "score" considers both rt and mz distances from expected values
            mz_score_list = []
            rt_core_list = []
            intensity_list = []
            similarity_between_list = []
            similarity_between = 0.0
            found_isotope = False
            for i in range(len(all_peak_objects)):
                this_peak_obj = all_peak_objects[i]
                # Check similarity between monoisotopic peak and the current one
                # The peaks will likely have different boundaries. Lets resize the isotope boundary to
                # match what we assume is the monoisotopic peak boundary -> at least for the comparison.
                mono_iso_intensity = pl.squeeze(pl.asarray(int_matrix[cur_peak_object.mzIndex,
                                                           cur_peak_object.rtMinIndex:cur_peak_object.rtMaxIndex].todense()))
                # Again, WE ARE USING THE (assumed) MONISOTOPIC PEAKS RT BOUNDRIES SO WE CAN DO A GOOD
                # SIMILARITY TEST.
                iso_peak_intensity = pl.squeeze(pl.asarray(
                    int_matrix[this_peak_obj.mzIndex, cur_peak_object.rtMinIndex:cur_peak_object.rtMaxIndex].todense()))
                similarity_between = cts.get_similarity(mono_iso_intensity, iso_peak_intensity)
                similarity_between_list.append(similarity_between)

                # Check MZ value of peak compared to mono. iso + H
                mz_diff_num_index = abs(best_mz_idex_guess - this_peak_obj.mzIndex)

                # Check retention time
                rt_diff_in_num_scans = abs(cur_peak_object.rtIndex - this_peak_obj.rtIndex)

                intensity_list.append(this_peak_obj.height)

                if mz_diff_num_index > mz_index_difference_threshold:
                    mz_score_list.append(0.0)
                    rt_core_list.append(0.0)
                    intensity_list.append(0.0)
                    similarity_between_list.append(0.0)
                    continue
                if rt_diff_in_num_scans > rt_scan_difference_threshold:
                    mz_score_list.append(0.0)
                    rt_core_list.append(0.0)
                    intensity_list.append(0.0)
                    similarity_between_list.append(0.0)
                    continue
                if cur_peak_object.height < this_peak_obj.height:
                    mz_score_list.append(0.0)
                    rt_core_list.append(0.0)
                    intensity_list.append(0.0)
                    similarity_between_list.append(0.0)
                    continue

                # do a linear scoring for rt and mz. 0 distance 1 -> furthest close to threshold values
                # is close to zero
                mz_slope = 1.0 / float(mz_index_difference_threshold + 1)
                rt_slope = 1.0 / float(rt_scan_difference_threshold + 1)
                mz_score = -mz_slope * mz_diff_num_index + 1.0
                rt_score = -rt_slope * rt_diff_in_num_scans + 1.0
                mz_score_list.append(mz_score)
                rt_core_list.append(rt_score)

                if visualize_post_initial_peaks:
                    pl.plot(mono_iso_intensity, c='r')
                    pl.plot(iso_peak_intensity, c='b')
                    pl.title("mono isotope (assumed) intensity (red)\n isotope int. (blue)")
                    pl.show()

            max_iso_intensity = max(intensity_list)
            best_score = 0.01
            for i in range(len(all_peak_objects)):
                this_peak_obj = all_peak_objects[i]
                mz_score = mz_score_list[i]
                rt_score = rt_core_list[i]
                intensity_score = intensity_list[i] / max_iso_intensity
                similarity_score = similarity_between_list[i]

                # Simple average of all different score used to compare peaks in patch
                # being investigated.
                total_score = (mz_score + rt_score + intensity_score + similarity_score) / 4.0
                # Keep track of the best score as you loop over the different peaks.
                if total_score > best_score:
                    best_score = total_score
                    found_isotope = True
                    best_peak = all_peak_objects[i]

                if VERBOSE:
                    print "---------------- isotope information for VERBOSE --------------------"
                    print "cur_peak_object.mzIndex: " + str(cur_peak_object.mzIndex)
                    print "cur_peak_object.rtMinIndex: " + str(cur_peak_object.rtMinIndex)
                    print "cur_peak_object.rtMaxIndex: " + str(cur_peak_object.rtMaxIndex)
                    print "this_peak_obj.mzIndex: " + str(this_peak_obj.mzIndex)
                    print "similarity_between: " + str(similarity_between)
                    print "best guess for what mz should be " + str(best_mz_idex_guess / 10000.0)
                    print "this_peak_obj.mz " + str(this_peak_obj.mz)
                    print "mz_diff_num_index " + str(mz_diff_num_index)
                    print "estimate_of_ms_fwhm " + str(estimate_of_ms_fwhm)
                    print "rt_diff_in_num_scans   " + str(rt_diff_in_num_scans)
                    print "cur_peak_object.rt " + str(cur_peak_object.rt)
                    print "this_peak_obj.rt " + str(this_peak_obj.rt)
                    print "mz_score " + str(mz_score)
                    print "rt_score " + str(rt_score)
                    print "intensity_score " + str(intensity_score)
                    print "similarity_score " + str(similarity_score)
                    print "total_score " + str(total_score)
                    print "---------------------------------------------------------------------"

            if found_isotope:
                cur_peak_object.addIsotope(best_peak)

    # just in case
    return cur_peak_object


def load_data_points_in_lists(dfr, absolute_intensity_thresh, mz_upper_cutoff):
    """
    Fills lists of mz values, RT values, and intensities as well as list of all data points

    :param dfr: Data file reader object
    :param absolute_intensity_thresh: Intensity below which everything is thrown out.
    :return: Lists in this order mz_by_scan, inten_by_scan, rt, list_all_data_points
    """

    mz_by_scan = []
    inten_by_scan = []
    rt = []
    scan_numbers = []
    list_all_data_points = []

    count = 0
    mz, inten = dfr.get_next_scan_mzvals_intensities()
    mz = mz[mz <= mz_upper_cutoff]
    inten = inten[mz <= mz_upper_cutoff]

    if not USE_SMALL_TEST_WINDOW:
        while mz is not None:
            # while (count<100):
            # line below is to skip chemical noise in data Change back to above line HERE
            # if (count>200)and(count<1200):
            sys.stdout.write("\r" + str(count))
            sys.stdout.flush()
            mz_by_scan.append(mz)
            inten_by_scan.append(inten)
            rt.append(dfr.get_rt_from_scan_num(count))
            scan_numbers.append(dfr.get_act_scan_num(count))
            for i in range(len(mz)):
                # Do not look for methabolites with m/z > this mz_upper_cutoff value.
                if (mz[i] > 0.0) and (mz[i] <= mz_upper_cutoff):
                    # if (mz[i]>0.0)and(mz[i]<1000.0):
                    if inten[i] < absolute_intensity_thresh:
                        continue
                    cur_dp = data_point.DataPoint(count, i, mz[i], inten[i])
                    list_all_data_points.append(cur_dp)

            mz, inten = dfr.get_next_scan_mzvals_intensities()
            mz = mz[mz <= mz_upper_cutoff]
            inten = inten[mz <= mz_upper_cutoff]

            count += 1
            # dfr.closeWriter()
    else:
        scan_index_count = 0
        while mz is not None:
            # while (count<100):
            # line below is to skip chemical noise in data Change back to aboveline HERE
            # if (count>200)and(count<1200):

            tf = (mz > MZ_MIN and mz <= MZ_MAX)
            mz = mz[tf]
            inten = inten[tf]

            cur_rt = dfr.get_rt_from_scan_num(count)
            if (cur_rt < (RT_MAX * 60.0)) and (cur_rt > (RT_MIN * 60.0)):
                #                sys.stdout.write("\r"+str(count))
                #                sys.stdout.flush()
                mz_by_scan.append(mz)
                inten_by_scan.append(inten)
                rt.append(cur_rt)
                scan_numbers.append(dfr.get_act_scan_num(count))
                for i in range(len(mz)):
                    # For testing HERE
                    if (mz[i] > MZ_MIN) and (mz[i] < MZ_MAX):
                        # if (mz[i]>0.0)and(mz[i]<1000.0):
                        if inten[i] < absolute_intensity_thresh:
                            continue
                        cur_dp = data_point.DataPoint(scan_index_count, i, mz[i], inten[i])
                        list_all_data_points.append(cur_dp)
                scan_index_count += 1

            mz, inten = dfr.get_next_scan_mzvals_intensities()

            count += 1
    return mz_by_scan, inten_by_scan, rt, scan_numbers, list_all_data_points


def get_unique_mz_values(mz_by_scan, rt, mz_upper_cutoff):
    """
    If we have a nearly perfect grid we can put all the intensities into a matrix. The only tricky
    part is finding all of the unique m/z values... This is tricky because no scan contains all m/z
    values -> lots of gaps. Look at all scans to find all mz values.

    :param mz_by_scan: list of mz values
    :param rt: list of retention time values
    :param mz_upper_cutoff: Example: probably don't want to look at m/z values above 1,000 if looking
    for metabolites. To disregard m/z values over 1,000 mz_upper_cutoff will be 1,000.
    :return: list of unique mz values
    """
    # Doing it with the dictionary and the try/except statement is WAAAYYYYY faster.
    unique_mzs = {}
    print "Building unique_mzs \n"
    print "..."
    for i in range(len(mz_by_scan)):
        # print ("cur scan: " + str(i))
        for j in range(len(mz_by_scan[i])):
            cur_mz = mz_by_scan[i][j]
            cur_mz = int(cur_mz * 10000)
            try:
                unique_mzs[cur_mz]
            except KeyError:
                unique_mzs[cur_mz] = True
    print "Done building unique_mzs \n"

    print ("len(unique_mzs): " + str(len(unique_mzs)))

    unique_mz_list = []

    if not USE_SMALL_TEST_WINDOW:
        for i in unique_mzs:

            # Do not look for methabolites with m/z > 1000.
            if i < (1000 * 10000):
                unique_mz_list.append(i)
    else:
        for i in unique_mzs:
            if (i < int(MZ_MAX * 10000)) and (i > int(MZ_MIN * 10000)):
                unique_mz_list.append(i)

    print "Full mz range of specified region:"
    print "     min(unique_mz_list): " + str(min(unique_mz_list) / 10000.0)
    print "     max(unique_mz_list): " + str(max(unique_mz_list) / 10000.0)
    print "     len(unique_mz_list): " + str(len(unique_mz_list))
    print "Full RT range of specified region:"
    print "     rt[0]: " + str(rt[0])
    print "     rt[-1]: " + str(rt[-1])
    print "     len(rt): " + str(len(rt))

    return unique_mz_list


def get_initial_peaks_and_peak_info(how_many_initial_peak_for_determining_parameters,
                                    eic_num_scans_leftright,
                                    int_matrix,
                                    to_mod_int_matrix,
                                    peakDetector,
                                    rt,
                                    estimate_of_ms_fwhm,
                                    visualize_initial_peaks,
                                    unique_mz_list,
                                    efficient_next_max
                                    ):
    """
    This looks at the entire data set and finds the most intense peaks within it.
    Starting with the highest intensity data point it takes a slice in RT and performs a wavelet
    transform to look for peaks. If it finds a peak at the maxima under investigation it looks at
    adjacent slices in the m/z domain to determine if the peak is good by comparing the similarity
    of those curves with the curve at the maxima. If not good peaks (at all) are found it sets all that
    data to zero in to_mod_int_matrix before performing a search for the next local maxima. If good peaks are
    found that are not at the most intense location, all the data EXCEPT the data corresponding to those good
    peaks is set to zero. A good peak is added to the result list before that data is set to zero. int_matrix
    is never modified so all EICs contain accurate description of data.

    This continues until it has found a given number of peaks. It stores certian properties (see return below)
    of the peaks that can be used to adaptively determine good parameters for the data set.

    :param how_many_initial_peak_for_determining_parameters: The number of peaks to find and to collect
        information from to determine peak picking parameters.
    :type how_many_initial_peak_for_determining_parameters: int
    :param eic_num_scans_leftright: How many scans in RT domain to build EIC around maxima.
    :type eic_num_scans_leftright: int
    :param int_matrix: Matrix of all intensity values (row->m/z, col->rt).
    :type int_matrix: dok_matrix
    :param to_mod_int_matrix: Matrix of all intensities that is slowly modified and set to zero as
        data is analyzed.
    :type to_mod_int_matrix: dok_matrix
    :param peakDetector: Peak detection object configured for finding peaks in EICs.
    :param rt: list of retention times.
    :type rt: list float
    :param estimate_of_ms_fwhm: Estimated value of the full width at half max of a mass spec in number of data points. scan peak.
    :type estimate_of_ms_fwhm: int
    :param visualize_initial_peaks: see the process?
    :type visualize_initial_peaks: bool
    :param unique_mz_list: (list int (floats rounded to 4th decimal place then times 10000 and converted to int))
        list of all unique mz values as integers
    :return: The detected peaks in Result() object form,
            fit parameters for all peaks (2d array),
            array of wavelet similarity values,
            array of coefficient over area values,
            array of the peak widths.
    """

    # store the peak width values for the first most intense peaks
    peak_width_arr = pl.array([])
    # store the coefficient over area values so that the coef over area threshold can be determined
    # adaptively
    coef_over_area_arr = pl.array([])
    # store similarities between wavelet and peak
    wavelet_peak_similarity_arr = pl.array([])
    # peak fit parameters for averaging
    peak_fit_params = pl.array([])

    inital_peaks_results = result.Result()
    num_good_peak_found = 0
    while num_good_peak_found < how_many_initial_peak_for_determining_parameters:
        sys.stdout.write("\r" + "initial peaks: " + str(num_good_peak_found) + " of " + str(
            how_many_initial_peak_for_determining_parameters))
        sys.stdout.flush()

        # indices of maximum. Use the modified matrix so that vales can be set to zero to find the next
        # maxima.
        a, b = efficient_next_max.find_max()

        # get the "EIC" associated with this maxima. Use the original intensity matrix so the eic does
        # not contain a discontinuity from the removal of data
        eic_left_bound = b - eic_num_scans_leftright
        eic_right_bound = b + eic_num_scans_leftright
        if eic_left_bound < 0:
            eic_left_bound = 0
        if eic_right_bound >= int_matrix.shape[1]:
            eic_right_bound = int_matrix.shape[1] - 1
        eic = int_matrix[a, eic_left_bound:eic_right_bound]
        eic = pl.squeeze(pl.asarray(eic.todense()))

        # Now do peaking. It is possible that several peaks will be found so you will have to check
        # each for te with boundaries that contain the index of the maximum found in the
        # "pl.unraveex..." line.

        peakDetector.setSignal(eic)
        peakDetector.setX(rt[eic_left_bound:eic_right_bound])
        all_left_bounds, all_right_bounds, all_peak_positions = peakDetector.findPeaks()
        all_coef_over_area_vals = peakDetector.getallCoefOverAreaVals()
        all_wavelet_similarity_vals = peakDetector.getallWaveletSimilarityVals()

        if VERBOSE: print ("len(all_right_bounds): " + str(len(all_right_bounds)))
        if VERBOSE:
            print "-------"
            for alpha in all_coef_over_area_vals:
                print "coef over area: " + str(alpha)
            print "-------"

        # If there are no good peaks then you can not use similarity between adjacent mz EICs to
        # determine boundary in m/z (for point removal). We will have to use the cts.estimate_fwhm_ms
        if (len(all_right_bounds) == 0) or (b not in (pl.array(all_peak_positions) + eic_left_bound)):
            if VERBOSE: print (
            "remove_bad_peak_info_best_we_can\n either because b not in (all_peak_positions+eic_left_bound)\n or no peak detected using original cwt peak detector.")
            remove_bad_peak_info_best_we_can_sparse(to_mod_int_matrix, estimate_of_ms_fwhm, eic_left_bound,
                                                    eic_right_bound, a, visualize_initial_peaks, efficient_next_max)

        else:
            # The bounds that are returned from this are the indices in the int_matrix
            rightB, leftB, cur_coef_over_area, cur_wav_peak_sim = get_peakwidth_and_coef_over_area(
                pl.asarray(int_matrix.todense()),
                all_peak_positions,
                all_right_bounds,
                all_left_bounds,
                eic_left_bound,
                a,
                b,
                visualize_initial_peaks,
                all_coef_over_area_vals,
                all_wavelet_similarity_vals)

            # Look at similarity for different mz slices -> tells you if its a good peak and what
            # the boundaries for removal in m/z are.
            # returns if it is a good peak, and the boundaries
            good_peak_bool, mz_up_bound, mz_low_bound, sim_vals = peak_similarity_test(int_matrix, rightB, leftB,
                                                                                       estimate_of_ms_fwhm, a,
                                                                                       visualize_initial_peaks)

            fit_parameters, cur_asym_similarity = cts.asymmetric_gaussian_fit(b - leftB,
                                                                              pl.asarray(int_matrix[a,
                                                                                         leftB:rightB].todense())[0],
                                                                              visualize_initial_peaks)

            # remove the peak and appropriate mz
            to_mod_int_matrix[mz_low_bound:mz_up_bound, leftB:rightB] = 0.0
            efficient_next_max.done_with_rows_cols(mz_low_bound, mz_up_bound, leftB, rightB)

            if visualize_initial_peaks:
                plot_color_map_peak_and_box_bounds(pl.asarray(int_matrix.todense()), rightB, leftB, mz_up_bound,
                                                   mz_low_bound)

            cur_peak_width = rightB - leftB
            if good_peak_bool:
                peak_width_arr = pl.append(peak_width_arr, cur_peak_width)
                coef_over_area_arr = pl.append(coef_over_area_arr, cur_coef_over_area)
                wavelet_peak_similarity_arr = pl.append(wavelet_peak_similarity_arr, cur_wav_peak_sim)

                num_good_peak_found += 1

                peak_mz = unique_mz_list[a] / 10000.0
                peak_rt = rt[b]

                if VERBOSE:
                    print ("peak_mz: " + str(peak_mz))
                    print ("peak_rt: " + str(peak_rt))
                    print ("peak left bound:  " + str(rt[leftB]))
                    print ("peak right bound: " + str(rt[rightB]))
                    print ("peak mz low bound:  " + str(unique_mz_list[mz_low_bound] / 10000.0))
                    print ("peak mz high bound: " + str(unique_mz_list[mz_up_bound] / 10000.0))

                inital_peaks_results.addRT(peak_rt)
                inital_peaks_results.addMZ(peak_mz)
                inital_peaks_results.addHeight(int_matrix[a, b])
                inital_peaks_results.addMaxMZ(unique_mz_list[mz_up_bound] / 10000.0)
                inital_peaks_results.addMinMZ(unique_mz_list[mz_low_bound] / 10000.0)
                inital_peaks_results.addMinRT(rt[leftB])
                inital_peaks_results.addMaxRT(rt[rightB])
                inital_peaks_results.addMZIndex(a)
                inital_peaks_results.addRTIndex(b)
                inital_peaks_results.addMZMinIndex(mz_low_bound)
                inital_peaks_results.addMZMaxIndex(mz_up_bound)
                inital_peaks_results.addRTMinIndex(leftB)
                inital_peaks_results.addRTMaxIndex(rightB)

                inital_peaks_results.addCoefOverArea(cur_coef_over_area)
                inital_peaks_results.addAsymGausFitVal(cur_asym_similarity)
                inital_peaks_results.addSimilarityAdjMZVals(sim_vals)

                # use the peak information (height, peakwidth etc) to do a fit and start
                # generating model parameters
                fit_parameters, curSimilarity = cts.asymmetric_gaussian_fit(b - leftB, pl.squeeze(
                    pl.asarray(int_matrix[a, leftB:rightB].todense())), visualize_initial_peaks)
                peak_fit_params = pl.append(peak_fit_params, pl.array(fit_parameters))

    return inital_peaks_results, peak_fit_params, wavelet_peak_similarity_arr, coef_over_area_arr, peak_width_arr


def convert_result_object_to_list_of_peaks(list_of_results):
    all_peak_objects = []

    if type(list_of_results) == list:
        for j in range(len(list_of_results)):
            cur_result = list_of_results[j]
            for i in range(len(cur_result.final_peak_list_rt_vals)):
                cur_peak = peak.Peak()
                cur_peak.setRTIndex(cur_result.rt_index[i])
                cur_peak.setMZIndex(cur_result.mz_index[i])
                cur_peak.setRT(cur_result.final_peak_list_rt_vals[i])
                cur_peak.setMZ(cur_result.final_peak_list_mz_vals[i])
                cur_peak.setHeight(cur_result.final_peak_list_peak_height_vals[i])
                cur_peak.setMZmin(cur_result.final_peak_list_mz_min_vals[i])
                cur_peak.setMZmax(cur_result.final_peak_list_mz_max_vals[i])
                cur_peak.setRTMin(cur_result.final_peak_list_peak_rt_start_vals[i])
                cur_peak.setRTMax(cur_result.final_peak_list_peak_rt_end_vals[i])

                cur_peak.setMZMinIndex(cur_result.mz_min_index[i])
                cur_peak.setMZMaxIndex(cur_result.mz_max_index[i])
                cur_peak.setRTMinIndex(cur_result.rt_min_index[i])
                cur_peak.setRTMaxIndex(cur_result.rt_max_index[i])

                # print "len(resultO.coef_over_areas): " + str(len(resultO.coef_over_areas))
                cur_peak.setCoefOverArea(cur_result.coef_over_areas[i])
                cur_peak.setAsymGausFit(cur_result.asym_gaus_fit_vals[i])
                cur_peak.setSimilarityAdjMZ(cur_result.similarity_adj_mz_vals[i])

                all_peak_objects.append(cur_peak)

    else:
        for i in range(len(list_of_results.final_peak_list_rt_vals)):
            cur_peak = peak.Peak()
            cur_peak.setRTIndex(list_of_results.rt_index[i])
            cur_peak.setMZIndex(list_of_results.mz_index[i])
            cur_peak.setRT(list_of_results.final_peak_list_rt_vals[i])
            cur_peak.setMZ(list_of_results.final_peak_list_mz_vals[i])
            cur_peak.setHeight(list_of_results.final_peak_list_peak_height_vals[i])
            cur_peak.setMZmin(list_of_results.final_peak_list_mz_min_vals[i])
            cur_peak.setMZmax(list_of_results.final_peak_list_mz_max_vals[i])
            cur_peak.setRTMin(list_of_results.final_peak_list_peak_rt_start_vals[i])
            cur_peak.setRTMax(list_of_results.final_peak_list_peak_rt_end_vals[i])

            cur_peak.setMZMinIndex(list_of_results.mz_min_index[i])
            cur_peak.setMZMaxIndex(list_of_results.mz_max_index[i])
            cur_peak.setRTMinIndex(list_of_results.rt_min_index[i])
            cur_peak.setRTMaxIndex(list_of_results.rt_max_index[i])

            # print "len(resultO.coef_over_areas): " + str(len(resultO.coef_over_areas))
            cur_peak.setCoefOverArea(list_of_results.coef_over_areas[i])
            cur_peak.setAsymGausFit(list_of_results.asym_gaus_fit_vals[i])
            cur_peak.setSimilarityAdjMZ(list_of_results.similarity_adj_mz_vals[i])

            all_peak_objects.append(cur_peak)

    return all_peak_objects


def main():
    parser = argparse.ArgumentParser()
    # data file
    parser.add_argument('-f', action='store', dest='f', type=str, required=True)
    # verbose?
    parser.add_argument('-v', action='store_true', dest='v', required=False)
    # this is the intensity below which everything is noise -> setting this to non-zero number
    # drastically speeds things up.
    parser.add_argument('--absoluteintensitythresh', action='store', dest='absoluteintensitythresh', type=float,
                        required=False, default=100.0)
    # minimum peak intensity threshold -> setting this can also speed things up quite a bit.
    parser.add_argument('--peakintensitythresh', action='store', dest='peakintensitythresh', type=float, required=False,
                        default=1000.0)
    # How many initial peaks used for calculating adaptive parameters?
    parser.add_argument('--numinitpeaks', action='store', dest='numinitpeaks', type=int, required=False, default=20)
    # Example: probably don't want to look at m/z values above 1,000 if looking
    # for metabolites. To disregard m/z values over 1,000 set this to 1,000.
    parser.add_argument('--mzupcutoff', action='store', dest='mzupcutoff', type=float, required=False, default=1000)
    # Output results directory paths. Just path
    parser.add_argument('-o', action='store', dest='outputpath', type=str, required=False, default="Results/")
    # Output results directory name
    parser.add_argument('-n', action='store', dest='outputname', type=str, required=False, default="LastRun/")

    inargs = parser.parse_args()

    df_str = inargs.f

    ##########################################################
    ######## Get the output directory set up first ###########
    ########     in case there is some problem     ###########
    ##########################################################
    output_results_path = inargs.outputpath
    if output_results_path[-1] != '/':
        output_results_path += '/'
    output_results_dir_name = inargs.outputname
    if output_results_dir_name[-1] != '/':
        output_results_dir_name += '/'

    result_dir_str = output_results_path + output_results_dir_name

    if os.path.exists(result_dir_str):
        shutil.rmtree(result_dir_str)
    os.mkdir(result_dir_str)

    ##########################################################
    #################### Things to set #######################
    ##########################################################
    global VERBOSE
    VERBOSE = inargs.v

    mz_upper_cutoff = inargs.mzupcutoff

    # All of the parameter will come from these peaks
    how_many_initial_peak_for_determining_parameters = inargs.numinitpeaks

    absolute_intensity_thresh = inargs.absoluteintensitythresh

    # you could also see what the 95th percentile intensity is to get

    ##########################################################
    ### First most intense peak parameters ###################
    ##########################################################
    # initial wavelet parameters
    initial_lowest_wavelet_scale = 1
    initial_highest_wavelet_scale = 10
    # minimum intensity threshold for first most intense peaks that are used in determining
    # parameters
    min_initial_intensity_thresh = 5000
    coef_over_area_initial_thresh = 100
    wavelet_peak_similarity_thresh = "off"  # ignore, done prior to asymmetric gaussian fit
    signal_noise_initial_thresh = "off"
    visualize_initial_peaks = False
    # number of scans used to calculate FWHM of MS peaks. Scans are selected randomly. The highest
    # intensity peak from each random scan is used for a single estimate of FWHM. ex. if this number is
    # 20 then the highest peaks fron 20 random scans are used to calculate FWHM.
    number_of_scans_for_fwhm_calc = 20

    ##########################################################
    ### Current parameters for peak picking after ############
    ### most intense peaks                        ############
    ##########################################################
    min_intensity_thresh = inargs.peakintensitythresh
    sig_noise_thresh = "off"
    visualize_post_initial_peaks = False

    # percent of average peak width (found from initial high intensity peaks)
    # used for determining allowance of peak widths lower bound. The upper bound is less strict and
    # depends on the variation so just use the maximum peak width of the good peaks value.
    # Maybe try using a smaller value here. This appears to be pretty big.
    low_bound_peak_width_percent = 0.75

    ###########################################################
    #################### Done with setting things #############
    ###########################################################






    ###########################################################
    ### Use initialized parameters. No setting in here ########
    ###########################################################
    # number of scans to the left and right of the maximum to be used in the EIC creation
    eic_num_scans_leftright = initial_highest_wavelet_scale * 5

    initial_wavelet_scale_increment = (initial_highest_wavelet_scale - initial_lowest_wavelet_scale) / 9.0

    # make the peak detector object that will find peaks in the EIC slices.
    peakDetector = peakdetector.PeakDetector()
    peakDetector.setWaveletSmallScale(initial_lowest_wavelet_scale)
    peakDetector.setWaveletLargeScale(initial_highest_wavelet_scale)
    peakDetector.setWaveletScaleIncrement(initial_wavelet_scale_increment)
    ###########################################################





    ###########################################################
    ######## Get the all of the data into useful format #######
    ###########################################################
    dfr = easyio.DataFileReader(df_str, True)

    mz_by_scan, inten_by_scan, rt, scan_numbers, list_all_data_points = load_data_points_in_lists(dfr,
                                                                                                  absolute_intensity_thresh,
                                                                                                  mz_upper_cutoff)

    print "Sorting the datapoints list"
    list_all_data_points.sort(key=lambda x: x.intensity, reverse=True)

    # Change to minutes
    rt = pl.array(rt) / 60.0

    # if we have a nearly perfect grid we can put all the intensities into a matrix. The only tricky
    # part is finding al of the unique m/z value... This is tricky because no scan contains all m/z
    # values -> lots of gaps. Look at all scans to find all mz values.
    unique_mz_list = get_unique_mz_values(mz_by_scan, rt, mz_upper_cutoff)

    unique_mz_list.sort()
    # map the values to the index they belong
    mz_to_index_map = {}
    for i in range(len(unique_mz_list)):
        mz_to_index_map[unique_mz_list[i]] = i

    # int_matrix = pl.zeros([len(unique_mz_list),len(mz_by_scan)])
    int_matrix = dok_matrix((len(unique_mz_list), len(mz_by_scan)), dtype=pl.float32)

    # Now that we have everything stored in a simple two dimensional matrix we can get the maximum
    # points really easy. We will want to set the local points near the maximum to zero so we can
    # use this same method to get the next maximum but we will also want to keep the original data
    # for use of boundary determination -> don't want discontinuities in the "EIC" peak picking.
    # to_mod_int_matrix -> to modify intensity matrix
    # to_mod_int_matrix = pl.copy(int_matrix)
    to_mod_int_matrix = dok_matrix((len(unique_mz_list), len(mz_by_scan)), dtype=pl.float32)

    print "Filling sparse matrices with all data\n"
    print "..."

    count_zero_intensity = 0
    countNonZeroIntensity = 0

    c = len(list_all_data_points)
    count = 0
    for i in list_all_data_points:
        if count % (c / 10) == 0:
            print "%.1f percent" % (float(count) / float(c) * 100.0)
        cur_mz = i.mz
        # multiplying by 10000 because with 4 digits of precision we can just represent everything as an
        # integer and divide by 10000 to get the actual number.
        cur_mz = int(cur_mz * 10000)
        cur_scan_index = i.scan_index
        cur_intensity = i.intensity

        cur_mz_index = mz_to_index_map[cur_mz]
        i.mz_index = cur_mz_index

        # set the intensities in the correct place in the matrix.
        int_matrix[cur_mz_index, cur_scan_index] = cur_intensity
        to_mod_int_matrix[cur_mz_index, cur_scan_index] = cur_intensity
        # if not USE_SMALL_TEST_WINDOW:
        #    int_matrix[cur_mz_index,cur_scan_index] = cur_intensity
        #    to_mod_int_matrix[cur_mz_index,cur_scan_index] = cur_intensity
        # else:
        #    int_matrix[cur_mz_index,cur_scan_index-int(MZ_MIN)+1] = cur_intensity
        #    to_mod_int_matrix[cur_mz_index,cur_scan_index-int(MZ_MIN)+1] = cur_intensity

        if cur_intensity < HP.get_epsilon():
            count_zero_intensity += 1
        else:
            countNonZeroIntensity += 1
        count += 1
    print "Done filling sparse matrices with all data\n"
    print "count_zero_intensity: " + str(count_zero_intensity)
    print "countNonZeroIntensity: " + str(countNonZeroIntensity)

    ###########################################################
    ######### Make the object which efficiently gets ##########
    #########          the next maximum              ##########
    ###########################################################
    efficient_next_max = efficient_find_next_max.EfficientNextMax(list_all_data_points)

    ###########################################################
    #### find an estimate of the FWHM of MS peaks #############
    #### and use it to determine some key params  #############
    ####           in the hard_parameters         #############
    ###########################################################
    # we will need this to remove points when we don't find a good peak from an EIC but there is a
    # maximum
    estimate_of_ms_fwhm = cts.estimate_fwhm_ms(int_matrix, number_of_scans_for_fwhm_calc)
    print ("estimate_of_ms_fwhm: " + str(estimate_of_ms_fwhm))
    HP.determine_num_similar_slices_from_estimated_ms_full_width_half_max(estimate_of_ms_fwhm)

    ###########################################################
    ######## Run first peak detection to get most #############
    ######## intense peaks for determining params #############
    ###########################################################
    # Keep setting the peak detector object up
    peakDetector.setVisualize(visualize_initial_peaks)
    peakDetector.setMinIntensityThreshold(min_initial_intensity_thresh)
    peakDetector.setCoefOverAreaThreshold(coef_over_area_initial_thresh)
    peakDetector.setSignalToNoiseThreshold(signal_noise_initial_thresh)
    peakDetector.setWaveletPeakHighestSimilarity(wavelet_peak_similarity_thresh)

    if not USE_HARD_CODED_DETECTION_PARAMETERS:
        # Store the peak width values for the first most intense peaks
        # ---> peak_width_arr = pl.array([])
        # Store the coefficient over area values so that the coef over area threshold can be determined adaptively
        # ---> coef_over_area_arr = pl.array([])
        # Store similarities between wavelet and peak
        # ---> wavelet_peak_similarity_arr = pl.array([])
        # Peak fit parameters for averaging
        # ---> peak_fit_params = pl.array([])
        # All detected peaks
        # ---> inital_peaks_results = result.Result()

        inital_peaks_results, peak_fit_params, wavelet_peak_similarity_arr, coef_over_area_arr, peak_width_arr \
            = get_initial_peaks_and_peak_info(how_many_initial_peak_for_determining_parameters,
                                              eic_num_scans_leftright,
                                              int_matrix,
                                              to_mod_int_matrix,
                                              peakDetector,
                                              rt,
                                              estimate_of_ms_fwhm,
                                              visualize_initial_peaks,
                                              unique_mz_list,
                                              efficient_next_max
                                              )

        ###########################################################
        #### Determine new parameters and do the rest of the ######
        ####                peak detection                   ######
        ###########################################################
        peak_fit_params = peak_fit_params.reshape(-1, 4)
        avgPeakFitParams = [peak_fit_params[:, 0].mean(), peak_fit_params[:, 1].mean(), peak_fit_params[:, 2].mean(),
                            peak_fit_params[:, 3].mean()]

        # From the peaks found above we can determine the best parameters
        avg_peak_width = peak_width_arr.mean()
        std_peak_width = peak_width_arr.std()

        avg_coef_over_area = coef_over_area_arr.mean()
        std_coef_over_area = coef_over_area_arr.std()

        coef_over_area_thresh = avg_coef_over_area / 1.5
        isotope_coef_over_area_thresh = avg_coef_over_area / 1.5

        # similarity with wavelet measures
        avg_sim = wavelet_peak_similarity_arr.mean()  # ignore
        std_sim = wavelet_peak_similarity_arr.std()  # ignore
        print "avg_sim: " + str(avg_sim)
        print "std_sim: " + str(std_sim)

        # bounds of allowed peak widths in units of number of scans
        peak_duration_range \
            = [int(round(avg_peak_width)) - int(round(avg_peak_width * low_bound_peak_width_percent)),
               int(max(peak_width_arr)) + int(round(avg_peak_width * low_bound_peak_width_percent))]

        print ("peak_width_arr: " + str(peak_width_arr))
        print ("avg_peak_width: " + str(avg_peak_width))
        print ("std_peak_width: " + str(std_peak_width))
        print ("coef_over_area_arr: " + str(coef_over_area_arr))
        print ("avg_coef_over_area: " + str(avg_coef_over_area))
        print ("std_coef_over_area: " + str(std_coef_over_area))
        print ("coef_over_area_thresh: " + str(coef_over_area_thresh))
        print ("peak_duration_range: " + str(peak_duration_range))

        # highest_wavelet_scale = avg_peak_width+std_peak_width
        # highest_wavelet_scale = avg_peak_width
        # highest_wavelet_scale = avg_peak_width-std_peak_width
        highest_wavelet_scale = avg_peak_width / 2.0
        print ("highest_wavelet_scale: " + str(highest_wavelet_scale))

        # end of "if USE_HARD_CODED_DETECTION_PARAMETERS:"
    # This section is mostly for testing. Could get rid of eventually or set up to be
    # easier to use.
    else:
        ###########################################################
        ####### YP01 Hard set parameters (mostly for testing) #####
        ###########################################################
        # peak_duration_range = [5,46]
        # isotope_coef_over_area_thresh = 351.34/1.5 # Average/1.5
        # highest_wavelet_scale = 9.975 # avg_peakwidth/2
        # if USE_ISOTOPE_PARAMETERS_FOR_ALL:
        #    coef_over_area_thresh = isotope_coef_over_area_thresh
        # elif not USE_ISOTOPE_PARAMETERS_FOR_ALL:
        #    coef_over_area_thresh = 351.34/1.5 # Average/1.5

        ###########################################################
        ####### DCSM Hard set parameters (mostly for testing) #####
        ###########################################################
        # peak_duration_range = [5,34]
        # isotope_coef_over_area_thresh = 160.48/1.5 # Average/1.5
        # highest_wavelet_scale = 8.25   #avg_peakwidth/2
        # if USE_ISOTOPE_PARAMETERS_FOR_ALL:
        #    coef_over_area_thresh = isotope_coef_over_area_thresh
        # elif not USE_ISOTOPE_PARAMETERS_FOR_ALL:
        #    coef_over_area_thresh = 160.48/1.5 # Average/1.5

        ###########################################################
        ####### Test data hard set parameters                 #####
        ###########################################################
        # peak_duration_range = [2, 100]
        # isotope_coef_over_area_thresh = 1.0  # Average/1.5
        # highest_wavelet_scale = 10.0  # avg_peakwidth/2
        # if USE_ISOTOPE_PARAMETERS_FOR_ALL:
        #     coef_over_area_thresh = isotope_coef_over_area_thresh
        # elif not USE_ISOTOPE_PARAMETERS_FOR_ALL:
        #     coef_over_area_thresh = 1.0  # Average/1.5

        ###########################################################
        ####### Single-cell data hard set parameters          #####
        ###########################################################
        peak_duration_range = [2, 100]
        isotope_coef_over_area_thresh = 89.0  # Average/1.5
        highest_wavelet_scale = 10.0  # avg_peakwidth/2
        if USE_ISOTOPE_PARAMETERS_FOR_ALL:
            coef_over_area_thresh = isotope_coef_over_area_thresh
        elif not USE_ISOTOPE_PARAMETERS_FOR_ALL:
            coef_over_area_thresh = 89.0  # Average/1.5






    # I currently believe this should always be 1.0. Re think this if we start dealing with peaks containing
    # many many more data points than we are currently dealing with.
    lowest_wavelet_scale = 1.0
    wavelet_scale_increment = 1

    # number of scans to the left and right of the maximum to be used in the EIC creation
    eic_num_scans_leftright = int(round(highest_wavelet_scale)) * 5

    ###########################################################
    ###### With the adaptively determined parameters ##########
    ######         detect the rest of the peaks      ##########
    ###########################################################

    # Make the peak detector that will find the peaks using the parameters
    # determined from above.
    the_peak_detector = peakdetector.PeakDetector()

    the_peak_detector.setVisualize(visualize_post_initial_peaks)
    the_peak_detector.setMinIntensityThreshold(min_intensity_thresh)
    the_peak_detector.setCoefOverAreaThreshold(coef_over_area_thresh)
    the_peak_detector.setSignalToNoiseThreshold(sig_noise_thresh)

    the_peak_detector.setWaveletSmallScale(lowest_wavelet_scale)
    the_peak_detector.setWaveletLargeScale(highest_wavelet_scale)
    # the_peak_detector.setWaveletPeakHighestSimilarity("off")
    the_peak_detector.setWaveletPeakHighestSimilarity(0.1)
    the_peak_detector.setWaveletScaleIncrement(wavelet_scale_increment)
    the_peak_detector.setPeakDurationRange(peak_duration_range)

    print "Peak detection"
    print "..."

    resultO = get_peak_list(the_peak_detector,
                            int_matrix,
                            to_mod_int_matrix,
                            eic_num_scans_leftright,
                            visualize_post_initial_peaks,
                            estimate_of_ms_fwhm,
                            min_intensity_thresh,
                            rt,
                            unique_mz_list,
                            efficient_next_max)

    RESULT_LIST.append(resultO)

    print "Done peak detection"

    ###########################################################
    ############# Find larger MZ isotopes of the peaks ########
    ###########################################################
    if not USE_HARD_CODED_DETECTION_PARAMETERS:
        RESULT_LIST.append(inital_peaks_results)

    all_peak_objects = convert_result_object_to_list_of_peaks(RESULT_LIST)

    # THESE ARE ALSO INCLUDED IN all_peak_objects
    if not USE_HARD_CODED_DETECTION_PARAMETERS:
        inital_peaks_objects = convert_result_object_to_list_of_peaks(inital_peaks_results)

    print "len(all_peak_objects) " + str(len(all_peak_objects))
    # print "type(all_peak_objects[0]): " + str(type(all_peak_objects[0]))

    # Make the isotope peak detector. This will be made with the parameters found
    # from the initial peaks. The differen
    isotope_peak_detector = peakdetector.PeakDetector()
    isotope_peak_detector.setVisualize(visualize_post_initial_peaks)
    isotope_peak_detector.setMinIntensityThreshold(0)
    isotope_peak_detector.setCoefOverAreaThreshold(isotope_coef_over_area_thresh)
    isotope_peak_detector.setSignalToNoiseThreshold(sig_noise_thresh)

    isotope_peak_detector.setWaveletSmallScale(lowest_wavelet_scale)
    isotope_peak_detector.setWaveletLargeScale(highest_wavelet_scale)
    isotope_peak_detector.setWaveletPeakHighestSimilarity("off")
    isotope_peak_detector.setWaveletScaleIncrement(wavelet_scale_increment)
    isotope_peak_detector.setPeakDurationRange(peak_duration_range)

    HP.set_use_isotope_parameters_for_all()

    print "Detecting isotopes of all detected peaks"
    print "..."

    for alphao in range(len(all_peak_objects)):
        cur_peak = all_peak_objects[alphao]
        # if ((cur_peak.mz)<(tmp_mz_of_peak+0.005))and((cur_peak.mz)>(tmp_mz_of_peak-0.005)):

        cur_peak = find_isotope_peak(cur_peak,
                                     int_matrix,
                                     isotope_peak_detector,
                                     eic_num_scans_leftright,
                                     visualize_post_initial_peaks,
                                     estimate_of_ms_fwhm,
                                     rt,
                                     unique_mz_list
                                     )
        ISOTOPE_RESULT_PEAK_LIST.append(cur_peak)

    # if cur_peak.hasIsotope and VERBOSE:
    #        print "printing peak and isotope info"
    #        print "peak MZ: " + str(cur_peak.mz)
    #        for betao in range(len(cur_peak.isotopeList)):
    #            print "iso "+str(betao)+" MZ: "+str(cur_peak.isotopeList[betao].mz)

    print "Done detecting isotopes of all detected peaks"

    HP.set_dont_use_isotope_parameters_for_all()

    print "len(ISOTOPE_RESULT_PEAK_LIST): " + str(len(ISOTOPE_RESULT_PEAK_LIST))
    print "len(RESULT_LIST): " + str(len(RESULT_LIST))
    print "result_dir_str: " + str(result_dir_str)
    print "len(rt): " + str(len(rt))
    print "len(scan_numbers): " + str(len(scan_numbers))

    ###########################################################
    ############### write results #############################
    ###########################################################

    write_csv.write(result_dir_str, ISOTOPE_RESULT_PEAK_LIST)
    # This writes the intensities and retention times for each peak so that the jython script can
    # turn the peak list into one that can be read by MZmine 2.
    write_all_peak_data.write(result_dir_str, ISOTOPE_RESULT_PEAK_LIST, pl.asarray(int_matrix.todense()), rt,
                              scan_numbers)

    # write some miscellaneous but useful information
    misc_out_f = open(result_dir_str + "/run_info.txt", "w")
    misc_out_f.write("Data file: " + str(df_str) + "\n")
    misc_out_f.write("Number of initial peaks for determining parameters: " + str(
        how_many_initial_peak_for_determining_parameters) + "\n")
    misc_out_f.write("Absolute intensity threshold: " + str(absolute_intensity_thresh) + "\n")
    misc_out_f.write("Initial lowest wavelet scale: " + str(initial_lowest_wavelet_scale) + "\n")
    misc_out_f.write("Initial highest wavelet scale: " + str(initial_highest_wavelet_scale) + "\n")
    misc_out_f.write("Initial intensity threshold: " + str(min_initial_intensity_thresh) + "\n")
    misc_out_f.write("Initial coefficient over area threshold: " + str(coef_over_area_initial_thresh) + "\n")
    misc_out_f.write("Initial wavelet peak similarity thresh: " + str(wavelet_peak_similarity_thresh) + "\n")
    misc_out_f.write("Initial signal to noise threshold: " + str(signal_noise_initial_thresh) + "\n")
    misc_out_f.write("How many scans used to calculate full width half max of mass spec peaks: " + str(
        number_of_scans_for_fwhm_calc) + "\n")
    misc_out_f.write("######################## Non-initial peak info ##############################\n")
    misc_out_f.write("Signal to noise threshold: " + str(sig_noise_thresh) + "\n")
    misc_out_f.write("Intensity threshold: " + str(min_intensity_thresh) + "\n")
    misc_out_f.write("Peak duration range: " + str(peak_duration_range) + "\n")
    misc_out_f.write("Isotope coefficient over area: " + str(isotope_coef_over_area_thresh) + "\n")
    misc_out_f.write("Coefficient over area: " + str(coef_over_area_thresh) + "\n")
    misc_out_f.write("Highest wavelet scale: " + str(highest_wavelet_scale) + "\n")
    misc_out_f.write("Lowest wavelet scale: " + str(lowest_wavelet_scale) + "\n")
    misc_out_f.write("Wavelet scale increment: " + str(wavelet_scale_increment) + "\n")
    misc_out_f.write("EIC number of scans left and right: " + str(eic_num_scans_leftright) + "\n")
    misc_out_f.close()

    ###########################################################
    ############### Plot Results ##############################
    ###########################################################
    if PLOT_ALL_PEAKS:
        if not USE_HARD_CODED_DETECTION_PARAMETERS:
            peak_plotting.plot_all_the_peaks(result_dir_str, inital_peaks_objects, ISOTOPE_RESULT_PEAK_LIST, int_matrix,
                                             rt)
        else:
            peak_plotting.plot_all_the_peaks(result_dir_str, [], ISOTOPE_RESULT_PEAK_LIST, int_matrix, rt)


if __name__ == "__main__":
    main()




# %run -d main.py -f /Users/xdu4/Documents/Duxiuxia/Analysis/my_projects/adap-3d/raw/DC_010814_StandardsMixTest1_34StandardMix_01.mzXML --absoluteintensitythresh 500 --peakintensitythresh 5000 --numinitpeaks 20 --mzupcutoff 400 -o /Users/xdu4/Documents/Duxiuxia/Analysis/my_projects/adap-3d -n results

# python main.py -f /Users/xdu4/Documents/Duxiuxia/Analysis/my_projects/adap-3d/raw/Ppyralis_spermatophore_pos_20uL.mzXML --absoluteintensitythresh 500 --peakintensitythresh 5000 --numinitpeaks 20 -o /Users/xdu4/Documents/Duxiuxia/Analysis/my_projects/adap-3d -n results

# python main.py -f /Users/xdu4/Documents/Duxiuxia/Analysis/my_projects/adap-3d/raw/test_out.CDF --absoluteintensitythresh 500 --peakintensitythresh 5000 --numinitpeaks 20 -o /Users/xdu4/Documents/Duxiuxia/Analysis/my_projects/adap-3d -n results

# python main.py -f /Users/xdu4/Documents/Duxiuxia/Analysis/my_projects/single_cell/0317EP D11-E7-TR1.mzXML --absoluteintensitythresh 10000 --peakintensitythresh 50000 --numinitpeaks 20 --mzupcutoff 200 -o /Users/xdu4/Documents/Duxiuxia/Analysis/my_projects/single_cell -n results