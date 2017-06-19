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


class EfficientNextMax():

    def __init__(self,data_point_list_in):
        self.sorted_data_point_list = []
        self.index_to_data_point_dic = {}
        # where are we in the list as far as finding the next max
        self.cur_max_index = 0

        self.set_up_list_and_dict(data_point_list_in)

    def set_up_list_and_dict(self, data_point_list_in):
        self.sorted_data_point_list = sorted(data_point_list_in, key=lambda x: x.intensity, reverse=True)
        for i in self.sorted_data_point_list:
            self.index_to_data_point_dic[(i.mz_index,i.scan_index)] = i

    # just like slicing this will be inclusive of the lower bound and exclusive of the higher.
    # i.e [low,up)
    def done_with_rows_cols(self, row_low, row_high, col_low, col_high):
        for i in range(row_low,row_high):
            for j in range(col_low,col_high):
                if (i, j) in self.index_to_data_point_dic:
                    self.index_to_data_point_dic[(i, j)].been_removed = True

    def put_back_rows_cols(self, row_low, row_high, col_low, col_high):
        for i in range(row_low,row_high):
            for j in range(col_low,col_high):
                if (i, j) in self.index_to_data_point_dic:
                    self.index_to_data_point_dic[(i, j)].been_removed = False

    def find_max(self):
        if len(self.sorted_data_point_list) == 0:
            return -1, -1

        cur_data_point = self.sorted_data_point_list[self.cur_max_index]
        while cur_data_point.been_removed:

            self.cur_max_index+=1
            if self.cur_max_index>=len(self.sorted_data_point_list):
                self.cur_max_index -= 1
                return -1, -1

            cur_data_point = self.sorted_data_point_list[self.cur_max_index]

        return cur_data_point.mz_index, cur_data_point.scan_index
