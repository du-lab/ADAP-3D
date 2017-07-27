from easyIOmassspec import easyio
from peakdetectioncwt import peakdetector
import matrix
import os

def main():

    parameters = {'absolute_intensity_thresh': 100.0,
                  'mz_factor': 10000.0,
                  'peak_range': [15, 5],
                  'scan_boundary': 20,
                  'mz_tolerance': 0.01,
                  'initial_lowest_wavelet_scale': 1,
                  'initial_highest_wavelet_scale': 10,
                  'initial_wavelet_scale_increment': 1,
                  'peak_intensity_threshold': 5000}

    input_dir = '/Users/jasonzhou/Desktop/adap-3d/project/Data/'

    peakDetector = peakdetector.PeakDetector()
    peakDetector.setWaveletSmallScale(parameters['initial_lowest_wavelet_scale'])
    peakDetector.setWaveletLargeScale(parameters['initial_highest_wavelet_scale'])
    peakDetector.setWaveletScaleIncrement(parameters['initial_wavelet_scale_increment'])
    peakDetector.setVisualize(False)

    for file_name in os.listdir(input_dir):
        if file_name == '.DS_Store':
            continue
        df_str = input_dir + file_name
        dfr = easyio.DataFileReader(df_str, False)
        matrix_variables = matrix.Matrix(dfr, parameters)

        left_bounds_checked = []
        right_bounds_checked = []
        scan_index_list = []

        for data_point in matrix_variables.list_all_data_points:

            if data_point.intensity < parameters['peak_intensity_threshold']:

                break

            if not data_point.been_removed:

                mz_index, scan_index = matrix_variables.find_max()

                int_mz_value = matrix_variables.index_to_mz(mz_index)
                mz_value = matrix_variables.index_to_mz(mz_index) / parameters['mz_factor']

                first_scan_boundary = scan_index - parameters['scan_boundary']
                second_scan_boundary = scan_index + parameters['scan_boundary']

                rt_array, inten_array = matrix_variables.construct_EIC(int_mz_value, first_scan_boundary, second_scan_boundary, parameters)

                peakDetector.setSignal(inten_array)
                peakDetector.setX(rt_array)
                all_left_bounds, all_right_bounds, all_peak_positions = peakDetector.findPeaks()

                if len(all_peak_positions) == 0:
                    matrix_variables.remove_cur_max(mz_index, scan_index, first_scan_boundary, second_scan_boundary)

                else:
                    for index,left_bound in enumerate(all_left_bounds):
                        if matrix_variables.rt[scan_index] > rt_array[left_bound] and matrix_variables.rt[scan_index] < rt_array[all_right_bounds[index]]:
                            left_bounds_checked.append(left_bound)
                            right_bounds_checked.append(all_right_bounds[index])
                            scan_index_list.append(scan_index)
                        else:
                            matrix_variables.remove_cur_max(mz_index, scan_index)



            else:
                continue

        for mz_scale in range(parameters['peak_range'][0], parameters['peak_range'][1], -1):
            expected_value_1 = mz_value - (mz_scale * 1.0033)
            expected_value_2 = mz_value + (mz_scale * 1.0033)


        stop = 3


if __name__ == "__main__":
    main()