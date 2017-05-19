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


class RequiredParameters:
    def __init__(self):
        self.epsilon = 1E-8
        self.use_isotope_parameters_for_all = False

        self.peak_params = MeasurementParameters()
        # similarity threshold for determining if two peaks are similar
        #self.peak_params.similarity_threshold = 0.5
        self.peak_params.similarity_threshold = 0.35 # Relaxed
        self.peak_params.asymmetric_fit_similarity_threshold = 0.25
        # the number of similar peaks on either side of the maximum peak required for a peak to be considered
        # a good peak. If this is set to two then a total of 5 peaks must be similar (2 on each side + 1 for
        # maximum peak
        # It turns out that profile data from different places can have some serious differences in the total number
        # of points per peak. Which means these need to be set from the data. Define some functions to do that based
        # off of the average full width at half max estimated from the data.
        self.peak_params.num_sim_either_side_threshold = -1
        self.peak_params.num_sim_total_threshold = -1

        self.iso_params = MeasurementParameters()
        self.iso_params.similarity_threshold = 0.35
        self.iso_params.asymmetric_fit_similarity_threshold = 0.0

        self.iso_params.num_sim_either_side_threshold = -1
        self.iso_params.num_sim_total_threshold = -1

    def get_epsilon(self):
        return self.epsilon

    def set_use_isotope_parameters_for_all(self):
        self.use_isotope_parameters_for_all = True

    def set_dont_use_isotope_parameters_for_all(self):
        self.use_isotope_parameters_for_all = False

    def get_peak_similarity_threshold(self):
        if self.use_isotope_parameters_for_all:
            return self.iso_params.similarity_threshold
        else:
            return self.peak_params.similarity_threshold

    def get_peak_asymmetric_fit_similarity_threshold(self):
        if self.use_isotope_parameters_for_all:
            return self.iso_params.asymmetric_fit_similarity_threshold
        else:
            return self.peak_params.asymmetric_fit_similarity_threshold

    def get_peak_num_sim_either_side_threshold(self):
        if self.use_isotope_parameters_for_all:
            return self.iso_params.num_sim_either_side_threshold
        else:
            return self.peak_params.num_sim_either_side_threshold

    def get_peak_num_sim_total_threshold(self):
        if self.use_isotope_parameters_for_all:
            return self.iso_params.num_sim_total_threshold
        else:
            return self.peak_params.num_sim_total_threshold

    def get_iso_similarity_threshold(self):
        return self.iso_params.similarity_threshold

    def get_iso_asymmetric_fit_similarity_threshold(self):
        return self.iso_params.asymmetric_fit_similarity_threshold

    def get_iso_num_sim_either_side_threshold(self):
        return self.iso_params.num_sim_either_side_threshold

    def get_iso_num_sim_total_threshold(self):
        return self.iso_params.num_sim_total_threshold

    def determine_num_similar_slices_from_estimated_ms_full_width_half_max(self,fwhm):
        # don't round or anything -> want it to be automatically rounded down for now.
        # These numbers have come from the extensive testing I have done on the DCSM
        # (standard mixture file) as well as the YP01 file (blood plasma).
        self.peak_params.num_sim_either_side_threshold = fwhm/2
        self.peak_params.num_sim_total_threshold = 2*self.peak_params.num_sim_either_side_threshold

        self.iso_params.num_sim_either_side_threshold = fwhm/4
        self.iso_params.num_sim_total_threshold = 1 + 2*self.iso_params.num_sim_either_side_threshold


class MeasurementParameters:

    def __init__(self):
        self.similarity_threshold = 0.0
        self.asymmetric_fit_similarity_threshold = 0.0
        self.num_sim_either_side_threshold = 0
        self.num_sim_total_threshold = 0

