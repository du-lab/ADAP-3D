import math
import os
import glob
import shutil
import time
import csv

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

from iso_dist_generator import get_carbon_num_to_mono_mass as mono_mass
from datamodel import matrix
from easyIOmassspec import easyio
from generalcurvetools import curve_tools
from bisect import bisect_left

path = '/Users/jasonzhou/Desktop/Working Folder/'
mz_constant_diff = 1.00335


def export_data(writer, mz_list, rt, first_ratio, second_ratio):

    row_to_write = []

    monoisotopic_mass = str(max(mz_list))
    ratio = str(first_ratio) + ':' + str(second_ratio)

    row_to_write.append(monoisotopic_mass)
    row_to_write.append(str(rt))
    row_to_write.append(mz_list)
    row_to_write.append(ratio)

    writer.writerows([row_to_write])


def plot_chromatogram(found_mz_list, mz_value, second_mz, bounded_rt_array, intensity_array_list, count,
                      scan_index, chromatogram_title, feature_title, matrix_variables, parameters):
    fig_combined = plt.figure()

    ax_combined = fig_combined.add_subplot(211)

    for intensity_array in intensity_array_list:
        ax_combined.plot(bounded_rt_array, intensity_array)

    ax_combined.set_title(chromatogram_title)

    ax_combined.set_xlabel("Retention Time")
    ax_combined.set_ylabel("Intensity")

    feature_ax = fig_combined.add_subplot(212)

    starting_mz_value = min(found_mz_list) - parameters['extension_of_feature']
    ending_mz_value = max(found_mz_list) + parameters['extension_of_feature']

    complete_mz_list = matrix_variables.mz_by_scan[scan_index]
    complete_intensity_list = matrix_variables.inten_by_scan[scan_index]

    for index, mz in enumerate(complete_mz_list):
        if mz > starting_mz_value:
            if mz < ending_mz_value:
                intensity = complete_intensity_list[index]
                if mz == mz_value or mz == second_mz:
                    feature_ax.plot([mz, mz], [0, intensity], 'b')
                elif mz in found_mz_list:
                    feature_ax.plot([mz, mz], [0, intensity], 'r')
                else:
                    feature_ax.plot([mz, mz], [0, intensity], 'k', alpha=0.3)
            else:
                break

    feature_ax.set_title(feature_title)

    feature_ax.set_xlabel("M/z")
    feature_ax.set_ylabel("Intensity")

    fig_combined.tight_layout()

    plt.savefig(path + 'comparing_peaks_results/' + str(count) + '.pdf')
    plt.close()


def get_norm_dot_product(experimental, theoretical):
    dot_product = np.dot(experimental, theoretical)

    exp_norm_factor = np.linalg.norm(experimental)
    theo_norm_factor = np.linalg.norm(theoretical)

    normalizing_factor = exp_norm_factor * theo_norm_factor

    norm_dot_product = dot_product / normalizing_factor

    return norm_dot_product


def get_carbon_number_factor(carbon_number, parameters):
    adjusted_carbon_number = carbon_number - 1
    exponent = -adjusted_carbon_number * parameters['carbon_number_parameter']
    inversed_factor = np.exp(exponent)
    carbon_number_factor = 1 - inversed_factor

    return carbon_number_factor


def find_mz_error(mz_value, mz_list, greater_mz, scan_index, matrix_variables):
    mz_to_error_dict = {}

    for mz_index in mz_list:
        mz = matrix_variables.index_to_data_point_dict[mz_index, scan_index].mz
        if greater_mz:
            scalar = round((mz_value - mz) / mz_constant_diff)
            correct_mz = mz_value - (scalar * mz_constant_diff)
            mz_error = abs(correct_mz - mz)
            mz_to_error_dict[mz_index] = mz_error
        else:
            scalar = round((mz - mz_value) / mz_constant_diff)
            correct_mz = mz_value + (scalar * mz_constant_diff)
            mz_error = abs(correct_mz - mz)
            mz_to_error_dict[mz] = mz_error

    return mz_to_error_dict


def remove_all_values(starting_mz_list, ending_mz_list, first_scan_boundary, second_scan_boundary, matrix_variables):
    for index, starting_index in enumerate(starting_mz_list):
        if starting_index is None:
            continue
        ending_index = ending_mz_list[index]
        matrix_variables.remove_cur_max(starting_index, ending_index, first_scan_boundary, second_scan_boundary + 1)


def gaussian(x, a, b, c):
    return a * np.exp(-(x - b) ** 2 / (2 * c))


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


def estimate_carbon_numbers(parameters):
    plt.figure()

    exp_c13_abundance = parameters['exp_c13_abundance']
    exp_c12_abundance = 1 - exp_c13_abundance

    smallest_range = parameters['peak_range'][1]
    largest_range = parameters['peak_range'][0]

    carbon_width_to_range_dict = {}

    for n in range(smallest_range, largest_range + 1):

        dist1 = []
        dist2 = []

        for k in range(0, n + 1):
            choose = math.factorial(n) / (math.factorial(n - k) * math.factorial(k))

            dist1_a = exp_c13_abundance ** k
            dist1_b = exp_c12_abundance ** (n - k)

            dist2_a = exp_c12_abundance ** k
            dist2_b = exp_c13_abundance ** (n - k)

            dist1.append(choose * dist1_a * dist1_b)
            dist2.append(choose * dist2_a * dist2_b)

        max_prob1 = 0
        max_index1 = None
        max_prob2 = 0
        max_index2 = None

        for index, prob in enumerate(dist1):

            if prob >= max_prob1:
                max_prob1 = prob
                max_index1 = index

        for index, prob in enumerate(dist2):

            if prob > max_prob2:
                max_prob2 = prob
                max_index2 = index

        width_difference = abs(max_index1 - max_index2)
        if width_difference in carbon_width_to_range_dict:
            carbon_width_to_range_dict[width_difference].append(n)
        else:
            carbon_width_to_range_dict[width_difference] = [n]

        plt.plot(n, width_difference, '.')

    plt.savefig(path + 'carbon_curve.pdf')
    plt.close()

    return carbon_width_to_range_dict


def create_theoretical_dist(carbon_range, first_ratio, second_ratio, parameters):
    exp_c13_abundance = parameters['exp_c13_abundance']
    exp_c12_abundance = 1 - exp_c13_abundance

    theoretical_dist = []

    for k in range(0, carbon_range + 1):
        choose = math.factorial(carbon_range) / (math.factorial(carbon_range - k) * math.factorial(k))

        dist1_a = exp_c13_abundance ** k
        dist1_b = exp_c12_abundance ** (carbon_range - k)

        dist2_a = exp_c12_abundance ** k
        dist2_b = exp_c13_abundance ** (carbon_range - k)

        dist1 = choose * dist1_a * dist1_b * first_ratio
        dist2 = choose * dist2_a * dist2_b * second_ratio

        theoretical_dist.append(dist1 + dist2)

    return theoretical_dist


def get_experimental_dist(carbon_range, carbon_width, mz_value, scan_index, data_point_dict, parameters, matrix_variables, neg):

    peak_check_range = (carbon_range - carbon_width) / 2

    number_of_found_points = 0

    mz_list = []
    experimental_intensity_list = []

    start = 0 - peak_check_range
    end = carbon_width + peak_check_range

    for scalar in range(start, end + 1):

        if scalar not in data_point_dict:
            if scalar in range(start, carbon_width + 1):
                data_point = None
            else:
                if neg:
                    expected_mz = mz_value - (scalar * mz_constant_diff)
                else:
                    expected_mz = mz_value + (scalar * mz_constant_diff)

                data_point = matrix_variables.check_existence(expected_mz, scan_index, need_mz_index=False)

        else:
            data_point = data_point_dict[scalar]

        if data_point is None:
            if neg:
                mz_list.append(mz_value - (scalar * mz_constant_diff))
            else:
                mz_list.append(mz_value + (scalar) * mz_constant_diff)
            experimental_intensity_list.append(0.0)

        else:
            if neg:
                mz_list.insert(0, data_point.mz)
                experimental_intensity_list.insert(0, data_point.intensity)
            else:
                mz_list.append(data_point.mz)
                experimental_intensity_list.append(data_point.intensity)

            number_of_found_points = number_of_found_points + 1

    if number_of_found_points / float(carbon_range + 1) <= parameters['found_values_between_peaks_threshold']:
        return None, None

    return mz_list, experimental_intensity_list


def get_data_point_dicts(mz_value, scan_index, matrix_variables, carbon_width_to_range_dict):

    max_range = max(carbon_width_to_range_dict)
    min_range = min(carbon_width_to_range_dict)

    neg_data_point_dict = {}
    pos_data_point_dict = {}

    for mz_scale in range(0, max_range + 1):

        neg_expected_value = mz_value - (mz_scale * mz_constant_diff)
        pos_expected_value = mz_value + (mz_scale * mz_constant_diff)

        neg_found_data_point = matrix_variables.check_existence(neg_expected_value, scan_index, need_mz_index=False)
        pos_found_data_point = matrix_variables.check_existence(pos_expected_value, scan_index, need_mz_index=False)

        if neg_found_data_point is not None:
            neg_data_point_dict[mz_scale] = neg_found_data_point

        if pos_found_data_point is not None:
            pos_data_point_dict[mz_scale] = pos_found_data_point

    if max(neg_data_point_dict) <= min_range:
        neg_data_point_dict = None

    if max(pos_data_point_dict) <= min_range:
        pos_data_point_dict = None

    return neg_data_point_dict, pos_data_point_dict


def get_similarity_and_EIC(mz, bounded_inten_array, first_scan_boundary, second_scan_boundary, matrix_variables, parameters):

    int_mz = mz * parameters['mz_factor']

    rt_array, intensity_array, mz_start, mz_end = matrix_variables.construct_EIC(int_mz, first_scan_boundary, second_scan_boundary)

    similarity = curve_tools.get_similarity(bounded_inten_array, intensity_array)

    return similarity, intensity_array, mz_start, mz_end


def find_second_peak(count, mz_value, scan_index, max_inten, carbon_width_to_range_dict,
                     carbon_number_to_mono_mass, bounded_inten_array, first_scan_boundary,
                     second_scan_boundary, matrix_variables, parameters):

    if os.path.exists(path + 'theoretical_distributions/' + str(count)):
        shutil.rmtree(path + 'theoretical_distributions/' + str(count))
    os.mkdir(path + 'theoretical_distributions/' + str(count))

    min_similarity = parameters['min_similarity']

    found_mz_list = None
    found_intensity_list = None

    found_second_mz = None
    found_second_intensity = None

    found_first_ratio = None
    found_second_ratio = None

    found_inten_array_list = None
    found_starting_index_list = None
    found_ending_index_list = None

    neg_data_point_dict, pos_data_point_dict = get_data_point_dicts(mz_value, scan_index, matrix_variables, carbon_width_to_range_dict)

    if neg_data_point_dict is not None:

        neg_similarity_to_intensity_dict = {}
        neg_intensity_array_list = []
        neg_starting_index_list = []
        neg_ending_index_list = []

        for mz_scale in neg_data_point_dict:

            neg_found_data_point = neg_data_point_dict[mz_scale]

            intensity = neg_found_data_point.intensity
            mz = neg_found_data_point.mz

            similarity, intensity_array, starting_index, ending_index = get_similarity_and_EIC(mz, bounded_inten_array,
                                                                                               first_scan_boundary,
                                                                                               second_scan_boundary,
                                                                                               matrix_variables,
                                                                                               parameters)

            neg_similarity_to_intensity_dict[similarity] = intensity

            good_similarity = True

            for similarity in neg_similarity_to_intensity_dict:
                if similarity < parameters['similarity_of_EIC_threshold']:
                    if neg_similarity_to_intensity_dict[similarity] > intensity * parameters['intensity_check_for_EIC']:
                        good_similarity = False

            if not good_similarity:
                continue

            neg_intensity_array_list.append(intensity_array)
            neg_starting_index_list.append(starting_index)
            neg_ending_index_list.append(ending_index)

            carbon_width = mz_scale

            if carbon_width < parameters['peak_range'][1]:
                continue

            first_ratio = intensity / (max_inten + intensity)
            second_ratio = max_inten / (max_inten + intensity)

            for carbon_range in carbon_width_to_range_dict[carbon_width]:

                mz_list, experimental_intensity_list = get_experimental_dist(carbon_range, carbon_width, mz_value, scan_index, neg_data_point_dict, parameters, matrix_variables, neg=True)

                if mz_list is None:
                    continue

                if intensity not in experimental_intensity_list:
                    print "The second peak intensity was not found within the experimental distribution. " \
                          "Check find_second_peak."

                if carbon_range in carbon_number_to_mono_mass:
                    if min(mz_list) < min(carbon_number_to_mono_mass[carbon_range]) or max(mz_list) > max(carbon_number_to_mono_mass[carbon_range]):
                        continue

                theoretical_dist = create_theoretical_dist(carbon_range, first_ratio, second_ratio, parameters)

                norm_dot_product = get_norm_dot_product(experimental_intensity_list, theoretical_dist)
                carbon_number_factor = get_carbon_number_factor(carbon_range, parameters)

                if norm_dot_product < parameters['dot_product_filter']:
                    continue

                inf_similarity = norm_dot_product * (1 - parameters['carbon_number_influence']) + (carbon_number_factor * parameters['carbon_number_influence'])

                fig = plt.figure()
                ax1 = fig.add_subplot(211)
                for index, inten in enumerate(theoretical_dist):
                    plt.plot([mz_list[index], mz_list[index]], [0, inten])
                ax1.set_title("Influenced Similarity: " + str(inf_similarity) + " | Theoretical")
                ax2 = fig.add_subplot(212)
                for index, inten in enumerate(experimental_intensity_list):
                    plt.plot([mz_list[index], mz_list[index]], [0, inten])
                ax2.set_title("Dot Product: " + str(norm_dot_product) + " | Experimental")
                plt.tight_layout()
                plt.savefig(path + 'theoretical_distributions/' + str(count) + '/' + str(mz_scale) + ' | ' + str(carbon_range))
                plt.close()

                if inf_similarity > min_similarity:

                    min_similarity = inf_similarity

                    found_mz_list = mz_list
                    found_intensity_list = experimental_intensity_list

                    found_second_mz = mz
                    found_second_intensity = intensity

                    found_first_ratio = first_ratio
                    found_second_ratio = second_ratio

                    found_inten_array_list = neg_intensity_array_list
                    found_starting_index_list = neg_starting_index_list
                    found_ending_index_list = neg_ending_index_list

    if pos_data_point_dict is not None:

        pos_similarity_to_intensity_dict = {}
        pos_intensity_array_list = []
        pos_starting_index_list = []
        pos_ending_index_list = []

        for mz_scale in pos_data_point_dict:

            pos_found_data_point = pos_data_point_dict[mz_scale]

            intensity = pos_found_data_point.intensity
            mz = pos_found_data_point.mz

            similarity, intensity_array, starting_index, ending_index = get_similarity_and_EIC(mz, bounded_inten_array,
                                                                                               first_scan_boundary,
                                                                                               second_scan_boundary,
                                                                                               matrix_variables,
                                                                                               parameters)

            pos_similarity_to_intensity_dict[similarity] = intensity

            good_similarity = True

            for similarity in pos_similarity_to_intensity_dict:
                if similarity < parameters['similarity_of_EIC_threshold']:
                    if pos_similarity_to_intensity_dict[similarity] > intensity * parameters['intensity_check_for_EIC']:
                        good_similarity = False

            if not good_similarity:
                continue

            pos_intensity_array_list.append(intensity_array)
            pos_starting_index_list.append(starting_index)
            pos_ending_index_list.append(ending_index)

            carbon_width = mz_scale

            if carbon_width < parameters['peak_range'][1]:
                continue

            first_ratio = max_inten / (max_inten + intensity)
            second_ratio = intensity / (max_inten + intensity)

            for carbon_range in carbon_width_to_range_dict[carbon_width]:

                mz_list, experimental_intensity_list = get_experimental_dist(carbon_range, carbon_width, mz_value,
                                                                             scan_index, pos_data_point_dict,
                                                                             parameters, matrix_variables, neg=False)

                if mz_list is None:
                    continue

                if intensity not in experimental_intensity_list:
                    print "The second peak intensity was not found within the experimental distribution. " \
                          "Check find_second_peak."

                if carbon_range in carbon_number_to_mono_mass:
                    if min(mz_list) < min(carbon_number_to_mono_mass[carbon_range]) or min(mz_list) > max(carbon_number_to_mono_mass[carbon_range]):
                        continue

                theoretical_dist = create_theoretical_dist(carbon_range, first_ratio, second_ratio, parameters)

                norm_dot_product = get_norm_dot_product(experimental_intensity_list, theoretical_dist)
                carbon_number_factor = get_carbon_number_factor(carbon_range, parameters)

                if norm_dot_product < parameters['dot_product_filter']:
                    continue

                inf_similarity = norm_dot_product * (1 - parameters['carbon_number_influence']) + (carbon_number_factor * parameters['carbon_number_influence'])

                fig = plt.figure()
                ax1 = fig.add_subplot(211)
                for index, inten in enumerate(theoretical_dist):
                    plt.plot([mz_list[index], mz_list[index]], [0, inten])
                ax1.set_title("Influenced Similarity: " + str(inf_similarity) + " | Theoretical")
                ax2 = fig.add_subplot(212)
                for index, inten in enumerate(experimental_intensity_list):
                    plt.plot([mz_list[index], mz_list[index]], [0, inten])
                ax2.set_title("Dot Product: " + str(norm_dot_product) + " | Experimental")
                plt.tight_layout()
                plt.savefig(
                    path + 'theoretical_distributions/' + str(count) + '/' + str(mz_scale) + ' | ' + str(carbon_range))
                plt.close()

                if inf_similarity > min_similarity:

                    min_similarity = inf_similarity

                    found_mz_list = mz_list
                    found_intensity_list = experimental_intensity_list

                    found_second_mz = mz
                    found_second_intensity = intensity

                    found_first_ratio = first_ratio
                    found_second_ratio = second_ratio

                    found_inten_array_list = pos_intensity_array_list
                    found_starting_index_list = pos_starting_index_list
                    found_ending_index_list = pos_ending_index_list

    if len(glob.glob(path + 'theoretical_distributions/' + str(count) + '/*')) == 0:
        shutil.rmtree(path + 'theoretical_distributions/' + str(count))

    if found_mz_list is not None:
        print min_similarity

    return found_mz_list, found_intensity_list, found_second_mz, found_second_intensity, found_first_ratio, found_second_ratio, found_inten_array_list, found_starting_index_list, found_ending_index_list





def main():
    time_1 = time.time()

    parameters = {'exp_c13_abundance': 0.05,
                  'absolute_intensity_thresh': 5000.0,
                  'mz_factor': 10000.0,
                  'peak_range': [35, 3],
                  'peak_check_threshold': 0.95,
                  'scan_boundary': 20,
                  'mz_tolerance': 0.001,
                  'peak_intensity_threshold': 10000,
                  'gaussian_error_tolerance': 0.15,
                  'gaussian_intensity_percentage': 0.05,
                  'dot_product_filter': 0.95,
                  'carbon_number_parameter': 0.05,
                  'carbon_number_influence': 0.5,
                  'min_similarity': 0.5,
                  'low_boundary_range': 0.02,
                  'high_boundary_range': 0.5,
                  'intensity_check_for_EIC': 0.2,
                  'similarity_of_EIC_threshold': 0.4,
                  'found_values_between_peaks_threshold': 0.5,
                  'extension_of_feature': 2}

    carbon_width_to_range_dict = estimate_carbon_numbers(parameters)

    carbon_number_to_mono_mass = mono_mass.get_mono_mass_dict()

    input_dir = path + 'Data/'

    if os.path.exists(path + 'comparing_peaks_results'):
        shutil.rmtree(path + 'comparing_peaks_results')

    os.mkdir(path + 'comparing_peaks_results')

    if os.path.exists(path + 'theoretical_distributions'):
        shutil.rmtree(path + 'theoretical_distributions')

    os.mkdir(path + 'theoretical_distributions')

    output_file = path + 'results'

    if os.path.exists(output_file):
        shutil.rmtree(output_file)

    os.mkdir(output_file)

    for file_name in os.listdir(input_dir):
        if file_name == '.DS_Store':
            continue

        unique_output_file = output_file + '/' + file_name

        if os.path.exists(unique_output_file):
            shutil.rmtree(unique_output_file)

        os.mkdir(unique_output_file)

        results_file = unique_output_file + '/results.csv'
        output = open(results_file, 'wb')
        writer = csv.writer(output)

        column_names = [['monoisotopic_mass', 'retention_time', 'm/z list', 'case_to_control']]
        writer.writerows(column_names)

        df_str = input_dir + file_name
        dfr = easyio.DataFileReader(df_str, False)
        matrix_variables = matrix.Matrix(dfr, parameters, mz_constant_diff)

        count = 0

        for data_point in matrix_variables.list_all_data_points:

            if data_point.intensity < parameters['peak_intensity_threshold']:
                break

            if not data_point.been_removed:

                mz_index = data_point.mz_index
                scan_index = data_point.scan_index

                mz_value = data_point.mz
                int_mz_value = mz_value * parameters['mz_factor']

                max_inten = data_point.intensity
                rt_of_max_inten = matrix_variables.rt[scan_index]

                first_scan_boundary = scan_index - parameters['scan_boundary']
                second_scan_boundary = scan_index + parameters['scan_boundary']

                first_scan_boundary = max(0, first_scan_boundary)
                second_scan_boundary = min(matrix_variables.int_matrix.shape[1] - 1, second_scan_boundary)

                rt_array, inten_array, original_mz_start, original_mz_end = matrix_variables.construct_EIC(int_mz_value,
                                                                                                           first_scan_boundary,
                                                                                                           second_scan_boundary)

                if len(rt_array) == 0 or sum(inten_array) == 0:
                    matrix_variables.remove_cur_max(mz_index, mz_index, scan_index, scan_index)
                    continue

                try:
                    popt, pcov = curve_fit(gaussian,
                                           rt_array,
                                           inten_array,
                                           bounds=([max_inten * 0.75, rt_of_max_inten - 0.01, 0],
                                                   [max_inten * 1.25, rt_of_max_inten + 0.01, np.inf]),
                                           absolute_sigma=False)

                except RuntimeError:
                    matrix_variables.remove_cur_max(mz_index, mz_index, scan_index, scan_index)
                    continue

                gaussian_values = gaussian(rt_array, max_inten, rt_of_max_inten, popt[2])
                delta = inten_array - gaussian_values
                error = np.mean(delta ** 2)
                norm_error = error / (max(inten_array) ** 2)

                if norm_error > parameters['gaussian_error_tolerance']:
                    matrix_variables.remove_cur_max(mz_index, mz_index, scan_index, scan_index)
                    continue

                rt_boundaries, boundary_width = calc_width_of_gaussian(popt, parameters)

                if boundary_width < parameters['low_boundary_range'] or boundary_width > parameters['high_boundary_range']:
                    matrix_variables.remove_cur_max(mz_index, mz_index, scan_index, scan_index)
                    continue

                first_index = bisect_left(rt_array, rt_boundaries[0])
                last_index = bisect_left(rt_array, rt_boundaries[1])

                bounded_rt_array = []
                bounded_inten_array = []

                for index in range(first_index, last_index):
                    bounded_rt_array.append(rt_array[index])
                    bounded_inten_array.append(inten_array[index])

                first_scan_boundary = bisect_left(matrix_variables.rt, min(bounded_rt_array))
                second_scan_boundary = bisect_left(matrix_variables.rt, max(bounded_rt_array))

                found_mz_list, found_intensity_list, second_mz, second_intensity, found_first_ratio, found_second_ratio, found_inten_array_list, found_starting_index_list, found_ending_index_list = find_second_peak(count, mz_value, scan_index, max_inten, carbon_width_to_range_dict, carbon_number_to_mono_mass, bounded_inten_array, first_scan_boundary, second_scan_boundary, matrix_variables, parameters)

                if found_mz_list is None:
                    matrix_variables.remove_cur_max(mz_index, mz_index, scan_index, scan_index)
                    continue

                chromatogram_title = 'Original: ' + str(mz_value) + ' | Expected: ' + str(second_mz)
                feature_title = 'RT: ' + str(rt_of_max_inten) + ' | Scan Index: ' + str(scan_index)

                plot_chromatogram(found_mz_list, mz_value, second_mz, bounded_rt_array, found_inten_array_list, count,
                                  scan_index, chromatogram_title, feature_title, matrix_variables, parameters)

                export_data(writer, found_mz_list, rt_of_max_inten, found_first_ratio, found_second_ratio)

                remove_all_values(found_starting_index_list, found_ending_index_list, first_scan_boundary,
                                  second_scan_boundary, matrix_variables)

                print count
                count = count + 1

    time_2 = time.time()

    print time_2 - time_1


if __name__ == "__main__":
    main()
