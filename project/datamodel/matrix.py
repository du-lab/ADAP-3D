from scipy.sparse import dok_matrix
import efficient_find_next_max
import pylab as pl
import data_point
import sys
import numpy as np

class Matrix():


    def __init__(self, dfr, parameters, mz_constant_diff):

        self.parameters = parameters

        self.mz_constant_diff = mz_constant_diff

        self.list_all_data_points = self.load_data_points_in_list(dfr)

        self.get_unique_mz_values()

        self.int_matrix = dok_matrix((len(self.unique_mz_list), len(self.mz_by_scan)), dtype=pl.float32)
#        self.to_mod_int_matrix = dok_matrix((len(self.unique_mz_list), len(self.mz_by_scan)), dtype=pl.float32)

        self.create_matrix()

        self.efficient_next_max = efficient_find_next_max.EfficientNextMax(self.list_all_data_points)

        self.index_to_data_point_dict = self.efficient_next_max.index_to_data_point_dic


    def load_data_points_in_list(self, dfr):

        absolute_intensity_thresh = self.parameters['absolute_intensity_thresh']

        self.mz_by_scan = []
        self.inten_by_scan = []
        self.rt = []
        self.scan_numbers = []
        list_all_data_points = []

        self.count = 0

        mz, inten = dfr.get_next_scan_mzvals_intensities()

        while mz is not None:
            sys.stdout.write("\r" + str(self.count))
            sys.stdout.flush()
            self.mz_by_scan.append(mz)
            self.inten_by_scan.append(inten)
            self.rt.append(dfr.get_rt_from_scan_num(self.count))
            self.scan_numbers.append(dfr.get_act_scan_num(self.count))
            for i in range(len(mz)):
                if inten[i] < absolute_intensity_thresh:
                    continue
                cur_dp = data_point.DataPoint(self.count, i, mz[i], inten[i])
                list_all_data_points.append(cur_dp)

            mz, inten = dfr.get_next_scan_mzvals_intensities()

            self.count += 1

        return list_all_data_points


    def get_unique_mz_values(self):

        self.unique_mzs = {}
        self.unique_mz_list = []

        print "Building unique_mzs \n"
        print "..."

        for i in range(len(self.mz_by_scan)):
            for j in range(len(self.mz_by_scan[i])):
                cur_mz = self.mz_by_scan[i][j]
                cur_mz = int(cur_mz * self.parameters['mz_factor'])
                try:
                    self.unique_mzs[cur_mz]
                except KeyError:
                    self.unique_mzs[cur_mz] = True

        print "Done building unique_mzs \n"

        print ("len(unique_mzs): " + str(len(self.unique_mzs)))

        for i in self.unique_mzs:
            self.unique_mz_list.append(i)

        print "Full mz range of specified region:"
        print "     min(unique_mz_list): " + str(min(self.unique_mz_list) / 10000.0)
        print "     max(unique_mz_list): " + str(max(self.unique_mz_list) / 10000.0)
        print "     len(unique_mz_list): " + str(len(self.unique_mz_list))
        print "Full RT range of specified region:"
        print "     rt[0]: " + str(self.rt[0])
        print "     rt[-1]: " + str(self.rt[-1])
        print "     len(rt): " + str(len(self.rt))


    def create_matrix(self):

        self.list_all_data_points.sort(key=lambda x: x.intensity, reverse=True)
        self.rt = pl.array(self.rt) / 60.0
        self.unique_mz_list.sort()

        self.mz_to_index_map = {}

        count_2 = 0
        c = len(self.list_all_data_points)

        for i in range(len(self.unique_mz_list)):
            self.mz_to_index_map[self.unique_mz_list[i]] = i

        for i in self.list_all_data_points:
            if count_2 % (c / 10) == 0:
                print "%.1f percent" % (float(count_2) / float(c) * 100.0)

            cur_mz = i.mz
            cur_mz = int(cur_mz * self.parameters['mz_factor'])
            cur_scan_index = i.scan_index
            cur_intensity = i.intensity

            cur_mz_index = self.mz_to_index_map[cur_mz]
            i.mz_index = cur_mz_index

            self.int_matrix[cur_mz_index, cur_scan_index] = cur_intensity
#            self.to_mod_int_matrix[cur_mz_index, cur_scan_index] = cur_intensity

            count_2 += 1


    def find_max(self):

        return self.efficient_next_max.find_max()


    def index_to_mz(self, index):

        return self.unique_mz_list[index]


    def iterate_unique_mz_list(self, first_mz_boundary, second_mz_boundary, first_scan_boundary, second_scan_boundary):

        mz_end_found = False
        mz_start_found = False

        for unique_mz in self.unique_mz_list:
            if unique_mz <= second_mz_boundary:
                if unique_mz >= first_mz_boundary:
                    for scan_index in range(first_scan_boundary, second_scan_boundary + 1):
                        if not mz_start_found:
                            possible_mz_start = self.mz_to_index_map[unique_mz]
                            try:
                                if not self.index_to_data_point_dict[possible_mz_start, scan_index].been_removed:
                                    mz_start = possible_mz_start
                                    mz_start_found = True
                            except KeyError:
                                pass
                        possible_mz_end = self.mz_to_index_map[unique_mz]
                        try:
                            if not self.index_to_data_point_dict[possible_mz_end, scan_index].been_removed:
                                mz_end = possible_mz_end
                                mz_end_found = True
                                break
                        except KeyError:
                            continue
                else:
                    continue
            else:
                break

        if mz_start_found is False or mz_end_found is False:
            return None, None, None

        inten_array = self.int_matrix[mz_start:mz_end + 1, first_scan_boundary:second_scan_boundary + 1].toarray().max(axis=0)

        return mz_start, mz_end, inten_array


    def construct_EIC(self, int_mz_value, first_scan_boundary, second_scan_boundary):

        mz_int_tolerance = self.parameters['mz_factor'] * self.parameters['mz_tolerance']

        first_mz_boundary = int_mz_value - mz_int_tolerance
        second_mz_boundary = int_mz_value + mz_int_tolerance

        mz_start, mz_end, inten_array = self.iterate_unique_mz_list(first_mz_boundary, second_mz_boundary, first_scan_boundary, second_scan_boundary)

        if mz_start is None or mz_end is None:
            return [], [], None, None

        rt_array = []

        for scan in range(first_scan_boundary, second_scan_boundary + 1):
            rt = self.rt[scan]
            rt_array.append(rt)

        return np.array(rt_array), inten_array, mz_start, mz_end


    def iterate_to_find_inbetween(self, mz_value, expected_mz_value, mz_range, scan_index):

        first_scan_boundary = scan_index
        second_scan_boundary = scan_index

        found_points = 1
        int_inbetween_mz_value_list = []

        for scalar in mz_range:

            inbetween_mz_value = mz_value + (scalar * self.mz_constant_diff)
            int_inbetween_mz_value = inbetween_mz_value * self.parameters['mz_factor']

            mz_int_tolerance = self.parameters['mz_factor'] * self.parameters['mz_tolerance']

            first_mz_boundary = int_inbetween_mz_value - mz_int_tolerance
            second_mz_boundary = int_inbetween_mz_value + mz_int_tolerance

            min_mz_index, max_mz_index, inten_array = self.iterate_unique_mz_list(first_mz_boundary, second_mz_boundary, first_scan_boundary, second_scan_boundary)

            if min_mz_index is None or max_mz_index is None:
                if inbetween_mz_value == expected_mz_value:
                    return 0, []
                continue

            elif max(inten_array) == 0:
                if inbetween_mz_value == expected_mz_value:
                    return 0, []
                continue

            found_points = found_points + 1
            int_inbetween_mz_value_list.append(int_inbetween_mz_value)

        int_inbetween_mz_value_list = sorted(int_inbetween_mz_value_list)

        return found_points, int_inbetween_mz_value_list

    def find_inbetween_mz_values(self, mz_value, expected_mz_value, scan_index, mz_scale):

        if mz_value > expected_mz_value:
            mz_range = range(0, mz_scale)

            found_point_value, int_inbetween_mz_value_list = self.iterate_to_find_inbetween(expected_mz_value, expected_mz_value, mz_range, scan_index)


        elif mz_value < expected_mz_value:
            mz_range = range(mz_scale, 0, -1)

            found_point_value, int_inbetween_mz_value_list = self.iterate_to_find_inbetween(mz_value, expected_mz_value, mz_range, scan_index)

        else:
            print "This shouldn't happen"
            print mz_value
            print expected_mz_value
            found_point_value = 0
            int_inbetween_mz_value_list = []

        if (found_point_value / mz_scale) < self.parameters['found_values_between_peaks_threshold']:
            return False, []
        else:
            return True, int_inbetween_mz_value_list


    def remove_cur_max(self, mz_start, mz_end, first_scan_boundary, second_scan_boundary):

        efficient_find_next_max.EfficientNextMax.done_with_rows_cols(self.efficient_next_max, mz_start, mz_end + 1, first_scan_boundary, second_scan_boundary + 1)