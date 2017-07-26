from scipy.sparse import dok_matrix
import efficient_find_next_max
import pylab as pl
import data_point
import sys

class Matrix():

    def __init__(self, dfr, parameters):

        self.parameters = parameters

        self.list_all_data_points = self.load_data_points_in_list(dfr)

        self.get_unique_mz_values()

        self.int_matrix = dok_matrix((len(self.unique_mz_list), len(self.mz_by_scan)), dtype=pl.float32)
        self.to_mod_int_matrix = dok_matrix((len(self.unique_mz_list), len(self.mz_by_scan)), dtype=pl.float32)

        self.create_matrix()

        self.efficient_next_max = efficient_find_next_max.EfficientNextMax(self.list_all_data_points)

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
                self.list_all_data_points.append(cur_dp)

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
            self.to_mod_int_matrix[cur_mz_index, cur_scan_index] = cur_intensity

            count_2 += 1

    def find_max(self):
        self.mz_index, self.scan_index = self.efficient_next_max.find_max()
