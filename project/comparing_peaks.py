from scipy.optimize import curve_fit
from easyIOmassspec import easyio
from generalcurvetools import curve_tools
from datamodel import matrix
import matplotlib.pyplot as plt
import numpy as np
import math
import os
import time

path = '/Users/jasonzhou/Desktop/Working Folder/'
mz_constant_diff = 1.00335

def plot_chromatogram(rt_array_list, inten_array_list, name, title, rt_of_max_inten):

    fig_combined = plt.figure()
    ax_combined = fig_combined.add_subplot(111)

    for index, array in enumerate(rt_array_list):
        ax_combined.plot(rt_array_list[index], inten_array_list[index])

    fig_combined.suptitle(title + ' | RT: ' + str(rt_of_max_inten))
    ax_combined.set_xlabel("Retention Time")
    ax_combined.set_ylabel("Intensity")
    plt.savefig(path + 'comparing_peaks_results/' + name + '.pdf')
    plt.close()


def gaussian(x, a, b, c):
    return a * np.exp(-(x-b) ** 2 / (2 * c))


def calc_width_of_gaussian(popt, parameters):

    variance = popt[2]
    mean = popt[1]

    gaussian_intensity_percentage = parameters['gaussian_intensity_percentage']
    ln_of_intensity_percentage = math.log(gaussian_intensity_percentage)

    radical = ln_of_intensity_percentage * -2 * variance

    range_value_1 = -1 * math.sqrt(radical)
    range_value_2 = math.sqrt(radical)

    rt_boundary_1 = range_value_1 + mean
    rt_boundary_2 = range_value_2 + mean
    rt_boundaries = [rt_boundary_1, rt_boundary_2]

    boundary_width = rt_boundary_2 - rt_boundary_1

    return rt_boundaries, boundary_width


def find_second_peak_and_remove_all_values(mz_value,
                                           first_scan_boundary,
                                           second_scan_boundary,
                                           scan_index,
                                           bounded_rt_array,
                                           bounded_inten_array,
                                           matrix_variables,
                                           rt_of_max_inten,
                                           parameters,
                                           count):

    for mz_scale in range(parameters['peak_range'][0], parameters['peak_range'][1], -1):

        name = str(count)
        title = 'Original: ' + str(mz_value) + ' | Expected: '

        expected_value_1 = mz_value - (mz_scale * mz_constant_diff)
        expected_value_2 = mz_value + (mz_scale * mz_constant_diff)

        int_expected_value_1 = expected_value_1 * parameters['mz_factor']
        int_expected_value_2 = expected_value_2 * parameters['mz_factor']

        found_rt_array_1, found_inten_array_1, mz_start_1, mz_end_1 = matrix_variables.construct_EIC(int_expected_value_1,
                                                                               first_scan_boundary,
                                                                               second_scan_boundary)

        found_rt_array_2, found_inten_array_2, mz_start_1, mz_end_1 = matrix_variables.construct_EIC(int_expected_value_2,
                                                                               first_scan_boundary,
                                                                               second_scan_boundary)

        if len(found_rt_array_1) == 0 and len(found_rt_array_2) == 0:
            continue

        elif len(found_rt_array_1) == 0:
            if sum(found_inten_array_2) == 0.0:
                continue
            sim = curve_tools.get_similarity(bounded_inten_array, found_inten_array_2)
            if sim > parameters['peak_similarity_threshold']:
                passes_found_point_threshold = matrix_variables.remove_inbetween_mz_values(mz_value, expected_value_2, first_scan_boundary, second_scan_boundary, scan_index, mz_scale)

                if not passes_found_point_threshold:
                    break

                title = title + str(expected_value_2)

                inten_array_list = [bounded_inten_array, found_inten_array_2]
                rt_array_list = [bounded_rt_array, found_rt_array_2]
                plot_chromatogram(rt_array_list, inten_array_list, name, title, rt_of_max_inten)

                break

        elif len(found_rt_array_2) == 0:
            if sum(found_inten_array_1) == 0.0:
                continue
            sim = curve_tools.get_similarity(bounded_inten_array, found_inten_array_1)
            if sim > parameters['peak_similarity_threshold']:
                passes_found_point_threshold = matrix_variables.remove_inbetween_mz_values(mz_value, expected_value_1, first_scan_boundary, second_scan_boundary, scan_index, mz_scale)

                if not passes_found_point_threshold:
                    break

                title = title + str(expected_value_1)

                inten_array_list = [bounded_inten_array, found_inten_array_1]
                rt_array_list = [bounded_rt_array, found_rt_array_1]
                plot_chromatogram(rt_array_list, inten_array_list, name, title, rt_of_max_inten)

                break

        else:
            if sum(found_inten_array_1) + sum(found_inten_array_2) == 0.0:
                continue

            sim_1 = curve_tools.get_similarity(bounded_inten_array, found_inten_array_1)
            sim_2 = curve_tools.get_similarity(bounded_inten_array, found_inten_array_2)

            if sim_1 > sim_2:
                if sim_1 > parameters['peak_similarity_threshold']:
                    passes_found_point_threshold = matrix_variables.remove_inbetween_mz_values(mz_value, expected_value_1, first_scan_boundary, second_scan_boundary, scan_index, mz_scale)

                    if not passes_found_point_threshold:
                        break

                    title = title + str(expected_value_1)

                    inten_array_list = [bounded_inten_array, found_inten_array_1]
                    rt_array_list = [bounded_rt_array, found_rt_array_1]
                    plot_chromatogram(rt_array_list, inten_array_list, name, title, rt_of_max_inten)

                    break

            elif sim_2 > sim_1:
                if sim_2 > parameters['peak_similarity_threshold']:
                    passes_found_point_threshold = matrix_variables.remove_inbetween_mz_values(mz_value, expected_value_2, first_scan_boundary,  second_scan_boundary, scan_index, mz_scale)

                    if not passes_found_point_threshold:
                        break

                    title = title + str(expected_value_2)

                    inten_array_list = [bounded_inten_array, found_inten_array_2]
                    rt_array_list = [bounded_rt_array, found_rt_array_2]
                    plot_chromatogram(rt_array_list, inten_array_list, name, title, rt_of_max_inten)

                    break

            else:
                if sim_1 > parameters['peak_similarity_threshold'] and sim_2 > parameters['peak_similarity_threshold']:

                    print "HOW DOES THIS HAPPEN? WHAT ARE THE ODDS?"
                    print sim_1
                    print sim_2
                    print mz_value
                    print expected_value_1
                    print expected_value_2

                    passes_found_point_threshold_1 = matrix_variables.remove_inbetween_mz_values(mz_value, expected_value_1, first_scan_boundary, second_scan_boundary, scan_index, mz_scale)
                    passes_found_point_threshold_2 = matrix_variables.remove_inbetween_mz_values(mz_value, expected_value_2, first_scan_boundary, second_scan_boundary, scan_index, mz_scale)

                    if not passes_found_point_threshold_1 or not passes_found_point_threshold_2:
                        break

                    title = title + 'First - ' + str(expected_value_1) + ', Second -' + str(expected_value_2)

                    inten_array_list = [bounded_inten_array, found_inten_array_1, found_inten_array_2]
                    rt_array_list = [bounded_rt_array, found_rt_array_1, found_rt_array_2]

                    plot_chromatogram(rt_array_list, inten_array_list, name, title, rt_of_max_inten)

                    break


def main():

    time_1 = time.time()

    parameters = {'absolute_intensity_thresh': 100.0,
                  'mz_factor': 10000.0,
                  'peak_range': [15, 5],
                  'scan_boundary': 20,
                  'mz_tolerance': 0.001,
                  'peak_intensity_threshold': 5000,
                  'gaussian_error_tolerance': 0.02,
                  'gaussian_intensity_percentage': 0.05,
                  'low_boundary_range': 0.02,
                  'high_boundary_range': 0.1,
                  'peak_similarity_threshold': 0.5,
                  'found_values_between_peaks_threshold': 0.66}

    input_dir = path + 'Data/'

    for file_name in os.listdir(input_dir):
        if file_name == '.DS_Store':
            continue

        df_str = input_dir + file_name
        dfr = easyio.DataFileReader(df_str, False)
        matrix_variables = matrix.Matrix(dfr, parameters, mz_constant_diff)

        count = 0

        for data_point in matrix_variables.list_all_data_points:

            if data_point.intensity < parameters['peak_intensity_threshold']:

                break

            if not data_point.been_removed:

                mz_index, scan_index = matrix_variables.find_max()

                int_mz_value = matrix_variables.index_to_mz(mz_index)
                mz_value = matrix_variables.index_to_mz(mz_index) / parameters['mz_factor']
                max_inten = matrix_variables.int_matrix[mz_index, scan_index]
                rt_of_max_inten = matrix_variables.rt[scan_index]

                first_scan_boundary = scan_index - parameters['scan_boundary']
                second_scan_boundary = scan_index + parameters['scan_boundary']

                first_scan_boundary = max(0, first_scan_boundary)
                second_scan_boundary = min(matrix_variables.int_matrix.shape[1] - 1, second_scan_boundary)

                rt_array, inten_array, mz_start, mz_end = matrix_variables.construct_EIC(int_mz_value, first_scan_boundary, second_scan_boundary)

                if len(rt_array) == 0:
                    matrix_variables.remove_cur_max(mz_index, mz_index + 1, first_scan_boundary, second_scan_boundary)
                    continue

                if sum(inten_array) == 0.0:
                    matrix_variables.remove_cur_max(mz_start, mz_end, first_scan_boundary, second_scan_boundary)
                    continue

                try:
                    popt, pcov = curve_fit(gaussian,
                                           rt_array,
                                           inten_array,
                                           bounds=([max_inten * 0.75, rt_of_max_inten - 0.01, 0],
                                                   [max_inten * 1.25, rt_of_max_inten + 0.01, np.inf]),
                                           absolute_sigma=False)

                except RuntimeError:
                    matrix_variables.remove_cur_max(mz_start, mz_end, first_scan_boundary, second_scan_boundary)
                    continue

                # errors = np.sqrt(np.diag(pcov))
                # c_value = errors[2]
                #
                # if c_value > parameters['c_value_filter']:
                #     matrix_variables.remove_cur_max(mz_start, mz_end, mz_index, first_scan_boundary, second_scan_boundary)
                #     continue

                # fig = plt.figure()
                # axes = fig.add_subplot(111)
                # axes.plot(rt_array, inten_array)
                # fig.suptitle(str(mz_value))
                # plt.savefig('/Users/jasonzhou/Desktop/Testing/' + str(count) + '.pdf')
                # plt.close()

                gaussian_values = gaussian(rt_array, max_inten, rt_of_max_inten, popt[2])
                delta = inten_array - gaussian_values
                error = np.mean(delta ** 2)
                norm_error = error / (max(inten_array) ** 2)

                if norm_error > parameters['gaussian_error_tolerance']:
                    matrix_variables.remove_cur_max(mz_start, mz_end, first_scan_boundary, second_scan_boundary)
                    continue

                rt_boundaries, boundary_width = calc_width_of_gaussian(popt, parameters)

                if boundary_width < parameters['low_boundary_range'] or boundary_width > parameters['high_boundary_range']:
                    matrix_variables.remove_cur_max(mz_start, mz_end, first_scan_boundary, second_scan_boundary)
                    continue

                else:

                    bounded_rt_array = []
                    bounded_inten_array = []

                    for index, rt in enumerate(rt_array):
                        if rt > rt_boundaries[0] and rt < rt_boundaries[1]:
                            bounded_rt_array.append(rt)
                            bounded_inten_array.append(inten_array[index])

                    for rt_index, listed_rt in enumerate(matrix_variables.rt):
                        if listed_rt == min(bounded_rt_array):
                            first_scan_boundary = rt_index
                        elif listed_rt == max(bounded_rt_array):
                            second_scan_boundary = rt_index
                            break
                        else:
                            continue

                    find_second_peak_and_remove_all_values(mz_value,
                                                           first_scan_boundary,
                                                           second_scan_boundary,
                                                           scan_index,
                                                           bounded_rt_array,
                                                           bounded_inten_array,
                                                           matrix_variables,
                                                           rt_of_max_inten,
                                                           parameters,
                                                           count)

                    matrix_variables.remove_cur_max(mz_start, mz_end, first_scan_boundary, second_scan_boundary)

                    print count
                    count = count + 1

                    continue

            else:
                continue

    time_2 = time.time()

    print time_2 - time_1




if __name__ == "__main__":
    main()