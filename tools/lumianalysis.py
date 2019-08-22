from tools.luminometer import Luminometer as L
from tools.detectorsratio import DetectorsRatio as Ratios
from tools.detectorsratio import MultipleDetectorsRatio as MRatios
from tools.linearityanalysis import LinearityAnalysis
import settings as setts
from tools.by_nls_test import ByNlsTest as BNLS
from tools import full_run_utilities as frutils
import numpy as np


class LAnalysis:
    def __init__(self, dets_file_labels: list, input_dir: str, vs_all_analysis=False,
                 run_linearity_analysis=False, mixed_data=False, run_stddev_test=False, c_years = False,
                 exclusion = False, all_and_excluded_analysis=False) -> None:

        if exclusion:
            print('Executing exclusion mode ...\n')
            n_files = len(dets_file_labels)
            det_read = []
            years_and_dir = input_dir.split('/')
            year = years_and_dir[1]
            intersection = {}
            bad_runs_in_detectors = {}

            for k in range (0, n_files, 1):
                intersection[k] = []

            counts = 0
            datatype = type(1)

            for i in range (0, n_files - 1, 1):
                for j in range (i + 1, n_files, 1):
                    input_path = 'plots/' + year + '/' + dets_file_labels[i] + '-' + dets_file_labels[j] +\
                                 '/txt/Bad_runs.txt'
                    try:
                        det_read.append(np.loadtxt(input_path, dtype=datatype))
                    except IOError as errIO:
                        print('\n Missing ' + dets_file_labels[i] + '-' + dets_file_labels[j] + ' bad runs \n')
                        det_missing = [dets_file_labels[i], dets_file_labels[j]]
                        LAnalysis(det_missing, input_dir=input_dir, mixed_data=mixed_data,
                                  run_stddev_test=run_stddev_test, c_years=c_years, exclusion=False)
                        det_read.append(np.loadtxt(input_path, dtype=datatype))

                    intersection[i].append(counts)
                    intersection[j].append(counts)
                    counts = counts + 1

            for i in range (0, n_files, 1):
                bad_runs_in_detectors[i] = np.intersect1d(det_read[intersection[i][0]], det_read[intersection[i][1]])
                for j in range (2, len(intersection[i]), 1):
                    bad_runs_in_detectors[i] = np.intersect1d(bad_runs_in_detectors[i], det_read[intersection[i][j]])

            filePath = 'plots/'
            txtfileName = 'Bad_runs_in_detectors'
            fileout = open(filePath + txtfileName + ".txt", "w+")
            for i in range(0, n_files, 1):
                fileout.write(dets_file_labels[i] + ': ' + str(bad_runs_in_detectors[i]) + '\n')
            fileout.close()

            quit()

        n_files = 0
        if all_and_excluded_analysis:
            all_data_detcs = []
        else:
            all_data_detcs = None

        if c_years:
            years_and_dir = input_dir.split(',')
            n_files = len(dets_file_labels)
            files_path = []
            detcs = []
            for i in range(1, len(years_and_dir), 1):
                files_path.append(years_and_dir[0] + years_and_dir[i] + '/')

            for det in dets_file_labels:
                detcs.append(frutils.merge_lumies(det, files_path, mixed_data=mixed_data))
        else:
            n_files = len(dets_file_labels)
            detcs = []
            for det in dets_file_labels:
                try:
                    detcs.append(L(det, input_dir + det + ".csv", mixed_data=mixed_data, full_data_analysis=all_and_excluded_analysis))

                except IOError as errIO:
                    print(errIO)
                    print('Please check if default input folder is correctly created: ' + setts.csv_input_base_dir)
                    print('Also check that your .csv file is in the correct year folder: ' + input_dir)
                    raise

        if n_files == 2:
            print("******** 2 detector comparison choose!! ******")

            if c_years:
                years_and_dir = input_dir.split(',')
                years = years_and_dir[1]
                for i in range(2, len(years_and_dir), 1):
                    years = years + ',' + years_and_dir[i]
                ratios12 = Ratios(detcs[0], detcs[1], year=years, c_years=c_years)
            else:
                ratios12 = Ratios(detcs[0], detcs[1])

            if run_stddev_test:
                #test = BNLS(detcs[0], detcs[1], array_by_step = 5)
                test = BNLS(detcs[0], detcs[1])
                # Standard deviation test
                test.plot_stddev_vs_nls()
                test.plot_nls_histograms()
                test.save_plots()
                quit()

            # Fill stability plots
            ratios12.plot_ratio_vs_date()
            ratios12.plot_ratio_vs_lumi2()
            ratios12.plot_nls_ratio_vs_lumi2()
            ratios12.plot_ratio_hist()
            ratios12.plot_nls_ratio_hist()
            ratios12.plot_ratio_hist_weighted()
            ratios12.plot_nls_ratio_hist_weighted()

            # Extra
            ratios12.plot_nls_ratio_vs_run()
            ratios12.plot_nls_ratio_vs_fill()
            ratios12.plot_ratio_vs_run()
            ratios12.plot_ratio_vs_fill()
            ratios12.plot_bad_fills()
            ratios12.plot_bad_runs()

            # all/excluded analysis plots
            ratios12.plot_all_and_excluded_by_detc()

            ratios12.save_plots()

            if run_linearity_analysis:
                # Linearity
                lin_analysis = LinearityAnalysis(ratios12)
                lin_analysis.plot_hist_sbil()
                lin_analysis.plot_ratio_vs_all_data_sbil()
                lin_analysis.save_plots()

        elif n_files == 3:
            print("******** 3 detector comparison choose!! ******")

            if c_years:
                years_and_dir = input_dir.split(',')
                years = years_and_dir[1]
                for i in range(2, len(years_and_dir), 1):
                    years = years + ',' + years_and_dir[i]
                ratios123 = MRatios(detcs[0], detcs[1], detcs[2], year=years, c_years=c_years)
            else:
                ratios123 = MRatios(detcs[0], detcs[1], detcs[2])

            # Fill plots
            ratios123.plot_nls_ratios_vs_date()
            ratios123.plot_ratios_vs_date()

            ratios123.plot_ratios_vs_lumi3()
            ratios123.plot_nls_ratios_vs_lumi3()

            # Mine
            # ratios123.plot_nls_ratios_vs_run()
            # ratios123.plot_nls_ratios_vs_fill()
            ratios123.plot_ratios_vs_run()
            # ratios123.plot_ratios_vs_fill()

            # Save plots
            ratios123.save_plots()

# if n_files == 2:
#
#     try:
#         detc1 = L(det1_label, base_input_path + det1_label + ".csv")
#         detc2 = L(det2_label, base_input_path + det2_label + ".csv")
#     except IOError as errIO:
#         print(errIO)
#         print('Please check if default input folder is correctly created: ' + setts.csv_input_base_dir)
#         print('Also check that your .csv file is in the correct year folder: ' + base_input_path)
#         raise
#
#     ratios12 = Ratios(detc1, detc2)
#
#     # Fill plots
#     ratios12.plot_nls_ratio_hist()
#
#     # Save plots
#     ratios12.save_plots()
#
# elif n_files == 3:
#     print("******** 3 detector comparison choose!! ******")
#     det1_label = args[0]
#     det2_label = args[1]
#     det3_label = args[2]
#
#     try:
#         detc1 = L(det1_label, base_input_path + det1_label + ".csv")
#         detc2 = L(det2_label, base_input_path + det2_label + ".csv")
#         detc3 = L(det3_label, base_input_path + det3_label + ".csv")
#     except IOError as errIO:
#         print(errIO)
#         print('Please check if default input folder is correctly created: ' + setts.csv_input_base_dir)
#         print('Also check that your .csv file is in the correct year folder: ' + base_input_path)
#         raise
#
#
# else:
#     raise IOError('No detectors names entered')
