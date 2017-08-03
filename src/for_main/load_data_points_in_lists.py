import os, sys

cwd = os.getcwd()
sys.path.append(cwd + "/src/")

from easyIOmassspec import easyio


def load_data_points_in_lists(file_name, absolute_intensity_thresh, mz_upper_cutoff):
    """
    Fills lists of mz values, RT values, and intensities as well as list of all data points

    :param dfr: Data file reader object
    :param absolute_intensity_thresh: Intensity below which everything is thrown out.
    :return: Lists in this order mz_by_scan, inten_by_scan, rt, list_all_data_points
    """

    file_reader = easyio.DataFileReader(file_name, True)

    mz_by_scan = []
    inten_by_scan = []
    rt = []
    scan_numbers = []
    list_all_data_points = []

    count = 0
    mz, inten = file_reader.get_next_scan_mzvals_intensities()
    if not USE_SMALL_TEST_WINDOW:
        while mz is not None:
            # while (count<100):
            # line below is to skip chemical noise in data Change back to above line HERE
            # if (count>200)and(count<1200):
            sys.stdout.write("\r" + str(count))
            sys.stdout.flush()
            mz_by_scan.append(mz)
            inten_by_scan.append(inten)
            rt.append(file_reader.get_rt_from_scan_num(count))
            scan_numbers.append(file_reader.get_act_scan_num(count))
            for i in range(len(mz)):
                # Do not look for methabolites with m/z > this mz_upper_cutoff value.
                if (mz[i] > 0.0) and (mz[i] < mz_upper_cutoff):
                    # if (mz[i]>0.0)and(mz[i]<1000.0):
                    if inten[i] < absolute_intensity_thresh:
                        continue
                    cur_dp = data_point.DataPoint(count, i, mz[i], inten[i])
                    list_all_data_points.append(cur_dp)

            mz, inten = file_reader.get_next_scan_mzvals_intensities()

            count += 1
            # dfr.closeWriter()
    else:
        scan_index_count = 0
        while mz is not None:
            # while (count<100):
            # line below is to skip chemical noise in data Change back to aboveline HERE
            # if (count>200)and(count<1200):
            cur_rt = file_reader.get_rt_from_scan_num(count)
            if (cur_rt < (RT_MAX * 60.0)) and (cur_rt > (RT_MIN * 60.0)):
                #                sys.stdout.write("\r"+str(count))
                #                sys.stdout.flush()
                mz_by_scan.append(mz)
                inten_by_scan.append(inten)
                rt.append(cur_rt)
                scan_numbers.append(file_reader.get_act_scan_num(count))
                for i in range(len(mz)):
                    # For testing HERE
                    if (mz[i] > MZ_MIN) and (mz[i] < MZ_MAX):
                        # if (mz[i]>0.0)and(mz[i]<1000.0):
                        if inten[i] < absolute_intensity_thresh:
                            continue
                        cur_dp = data_point.DataPoint(scan_index_count, i, mz[i], inten[i])
                        list_all_data_points.append(cur_dp)
                scan_index_count += 1

            mz, inten = file_reader.get_next_scan_mzvals_intensities()

            count += 1
    return mz_by_scan, inten_by_scan, rt, scan_numbers, list_all_data_points

