import math
import os
import shutil
import time

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

from iso_dist_generator import get_carbon_num_to_mono_mass as mono_mass
from datamodel import matrix
from easyIOmassspec import easyio
from generalcurvetools import curve_tools

path = '/Users/jasonzhou/Desktop/Working Folder/'
mz_constant_diff = 1.00335


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
        ending_index = ending_mz_list[index]
        matrix_variables.remove_cur_max(starting_index, ending_index + 1, first_scan_boundary, second_scan_boundary + 1)


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


def get_experimental_dist(carbon_range, carbon_width, mz_value, second_mz, scan_value, matrix_variables):
    peak_check_range = (carbon_range - carbon_width) / 2

    if mz_value < second_mz:
        start_mz = mz_value - (peak_check_range * mz_constant_diff)
        return matrix_variables.check_for_experimental(start_mz, carbon_range, scan_value, lower_mz=True)

    elif mz_value > second_mz:
        start_mz = mz_value + (peak_check_range * mz_constant_diff)
        return matrix_variables.check_for_experimental(start_mz, carbon_range, scan_value, lower_mz=False)

    else:
        "This really should not be a problem."
        return None, None, None


def find_second_peak(count, mz_value, scan_index, max_inten, carbon_width_to_range_dict,
                     carbon_number_to_mono_mass, matrix_variables, parameters):
    if os.path.exists(path + str(count)):
        shutil.rmtree(path + str(count))
    os.mkdir(path + str(count))

    left_peak_exists = False
    left_check = mz_value - mz_constant_diff
    left_point = matrix_variables.check_existence(left_check, scan_index, need_mz_index=False)

    right_peak_exists = False
    right_check = mz_value + mz_constant_diff
    right_point = matrix_variables.check_existence(right_check, scan_index, need_mz_index=False)

    if left_point is not None:
        left_inten = left_point.intensity
        if left_inten > max_inten * parameters['peak_check_threshold']:
            left_peak_exists = True

    if right_point is not None:
        right_inten = right_point.intensity
        if right_inten > max_inten * parameters['peak_check_threshold']:
            right_peak_exists = True

    min_similarity = parameters['min_similarity']

    found_mz_list = None
    found_intensity_list = None
    found_mz_index_list = None

    found_second_mz = None
    found_second_intensity = None

    for mz_scale in range(max(carbon_width_to_range_dict) * -1, max(carbon_width_to_range_dict) + 1):

        if mz_scale in range((parameters['peak_range'][1] * -1) + 1, parameters['peak_range'][1]):
            continue

        expected_value = mz_value + (mz_scale * mz_constant_diff)

        found_data_point = matrix_variables.check_existence(expected_value, scan_index, need_mz_index=False)

        if found_data_point is None:
            continue

        carbon_width = abs(mz_scale)

        second_intensity = found_data_point.intensity
        second_mz = found_data_point.mz

        if mz_scale < 0:
            if left_peak_exists:
                try:
                    if len(carbon_width_to_range_dict[carbon_width - 1]) > 1:
                        carbon_width = carbon_width - 1
                except KeyError:
                    pass
            first_ratio = second_intensity / (max_inten + second_intensity)
            second_ratio = max_inten / (max_inten + second_intensity)

        elif mz_scale > 0:
            if right_peak_exists:
                try:
                    if len(carbon_width_to_range_dict[carbon_width - 1]) > 1:
                        carbon_width = carbon_width - 1
                except KeyError:
                    pass
            first_ratio = max_inten / (max_inten + second_intensity)
            second_ratio = second_intensity / (max_inten + second_intensity)

        else:
            print "This really should not be a problem."
            continue

        for carbon_range in carbon_width_to_range_dict[carbon_width]:

            mz_list, experimental_intensity_list, mz_index_list = get_experimental_dist(carbon_range, carbon_width,
                                                                                        mz_value, second_mz,
                                                                                        scan_index, matrix_variables)

            if mz_list is None:
                continue

            if second_intensity not in experimental_intensity_list:
                print "The second peak intensity was not found within the experimental distribution." \
                      "Check find_second_peak."

            if carbon_range in carbon_number_to_mono_mass:
                if min(mz_list) < min(carbon_number_to_mono_mass[carbon_range]) or min(mz_list > max(carbon_number_to_mono_mass[carbon_range])):
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
            plt.savefig(path + str(count) + '/' + str(mz_scale) + ' | ' + str(carbon_range))
            plt.close()

            if inf_similarity > min_similarity:
                min_similarity = inf_similarity

                found_mz_list = mz_list
                found_intensity_list = experimental_intensity_list
                found_mz_index_list = mz_index_list

                found_second_mz = second_mz
                found_second_intensity = second_intensity

                # fig = plt.figure()
                # fig.add_subplot(211)
                # plt.plot(mz_list, experimental_intensity_list, 'b')
                # fig.add_subplot(212).set_title('Theoretical')
                # plt.plot(mz_list, theoretical_dist, 'r')
                # plt.suptitle(str(similarity))
                # plt.tight_layout()
                # plt.show()
                # plt.close()

    if found_mz_list is not None:
        print min_similarity

    return found_mz_list, found_intensity_list, found_mz_index_list, found_second_mz, found_second_intensity


def check_EIC(mz_list, intensity_list, bounded_inten_array, second_intensity, first_scan_boundary, second_scan_boundary,
              matrix_variables, parameters):

    is_similar = True

    inten_array_list = []
    starting_index_list = []
    ending_index_list = []

    similarity_list = []

    for index, mz in enumerate(mz_list):

        if intensity_list[index] > second_intensity * parameters['intensity_check_for_EIC']:

            int_mz = mz * parameters['mz_factor']

            rt_array, inten_array, mz_start, mz_end = matrix_variables.construct_EIC(int_mz, first_scan_boundary,
                                                                                     second_scan_boundary)

            similarity = curve_tools.get_similarity(bounded_inten_array, inten_array)

            if similarity < parameters['similarity_of_EIC_threshold']:
                is_similar = False
                break

            inten_array_list.append(inten_array)
            starting_index_list.append(mz_start)
            ending_index_list.append(mz_end)

            similarity_list.append(similarity)

    if is_similar:
        return inten_array_list, starting_index_list, ending_index_list, similarity_list
    else:
        return None, None, None, None


def main():
    time_1 = time.time()

    parameters = {'exp_c13_abundance': 0.05,
                  'absolute_intensity_thresh': 5000.0,
                  'mz_factor': 10000.0,
                  'peak_range': [35, 3],
                  'peak_check_threshold': 0.95,
                  'scan_boundary': 20,
                  'mz_tolerance': 0.0005,
                  'peak_intensity_threshold': 5000,
                  'gaussian_error_tolerance': 0.1,
                  'gaussian_intensity_percentage': 0.05,
                  'dot_product_filter': 0.85,
                  'carbon_number_parameter': 0.05,
                  'carbon_number_influence': 0.5,
                  'min_similarity': 0.5,
                  'low_boundary_range': 0.02,
                  'high_boundary_range': 0.1,
                  'intensity_check_for_EIC': 0.2,
                  'similarity_of_EIC_threshold': 0.25,
                  'found_values_between_peaks_threshold': 0.67,
                  'extension_of_feature': 2}

    carbon_width_to_range_dict = estimate_carbon_numbers(parameters)

    carbon_number_to_mono_mass = mono_mass.get_mono_mass_dict()

    input_dir = path + 'Data/'

    if os.path.exists(path + 'comparing_peaks_results'):
        shutil.rmtree(path + 'comparing_peaks_results')

    os.mkdir(path + 'comparing_peaks_results')

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
                    matrix_variables.remove_cur_max(mz_index, mz_index + 1, first_scan_boundary, second_scan_boundary)
                    continue

                try:
                    popt, pcov = curve_fit(gaussian,
                                           rt_array,
                                           inten_array,
                                           bounds=([max_inten * 0.75, rt_of_max_inten - 0.01, 0],
                                                   [max_inten * 1.25, rt_of_max_inten + 0.01, np.inf]),
                                           absolute_sigma=False)

                except RuntimeError:
                    matrix_variables.remove_cur_max(mz_index, mz_index + 1, first_scan_boundary, second_scan_boundary)
                    continue

                gaussian_values = gaussian(rt_array, max_inten, rt_of_max_inten, popt[2])
                delta = inten_array - gaussian_values
                error = np.mean(delta ** 2)
                norm_error = error / (max(inten_array) ** 2)

                if norm_error > parameters['gaussian_error_tolerance']:
                    matrix_variables.remove_cur_max(mz_index, mz_index + 1, first_scan_boundary, second_scan_boundary)
                    continue

                rt_boundaries, boundary_width = calc_width_of_gaussian(popt, parameters)

                if boundary_width < parameters['low_boundary_range'] or boundary_width > parameters['high_boundary_range']:
                    matrix_variables.remove_cur_max(mz_index, mz_index + 1, first_scan_boundary, second_scan_boundary)
                    continue

                bounded_rt_array = []
                bounded_inten_array = []

                for index, rt in enumerate(rt_array):
                    if rt_boundaries[0] < rt < rt_boundaries[1]:
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

                found_mz_list, found_intensity_list, found_mz_index_list, second_mz, second_intensity = find_second_peak(
                    count, mz_value, scan_index, max_inten, carbon_width_to_range_dict, carbon_number_to_mono_mass,
                    matrix_variables, parameters)

                if found_mz_list is None:
                    matrix_variables.remove_cur_max(mz_index, mz_index + 1, first_scan_boundary, second_scan_boundary)
                    continue

                intensity_array_list, starting_mz_list, ending_mz_list, similarity_list = check_EIC(found_mz_list,
                                                                                                    found_intensity_list,
                                                                                                    bounded_inten_array,
                                                                                                    second_intensity,
                                                                                                    first_scan_boundary,
                                                                                                    second_scan_boundary,
                                                                                                    matrix_variables,
                                                                                                    parameters)

                if intensity_array_list is None:
                    matrix_variables.remove_cur_max(mz_index, mz_index + 1, first_scan_boundary, second_scan_boundary)
                    continue

                print similarity_list

                chromatogram_title = 'Original: ' + str(mz_value) + ' | Expected: ' + str(second_mz)
                feature_title = 'RT: ' + str(rt_of_max_inten) + ' | Scan Index: ' + str(scan_index)

                plot_chromatogram(found_mz_list, mz_value, second_mz, bounded_rt_array, intensity_array_list, count,
                                  scan_index, chromatogram_title, feature_title, matrix_variables, parameters)

                remove_all_values(starting_mz_list, ending_mz_list, first_scan_boundary, second_scan_boundary,
                                  matrix_variables)

                print count
                count = count + 1

    time_2 = time.time()

    print time_2 - time_1


if __name__ == "__main__":
    main()
