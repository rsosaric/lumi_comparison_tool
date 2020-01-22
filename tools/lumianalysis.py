from tools.luminometer import Luminometer as L
from tools.detectorsratio import DetectorsRatio as Ratios
from tools.detectorsratio import MultipleDetectorsRatio as MRatios
import settings as setts
from tools.by_nls_test import ByNlsTest as BNLS
from tools import full_run_utilities as frutils


class LAnalysis:
    def __init__(self, dets_file_labels: list, input_dir: str, lumi_type: str,
                 mixed_data=False, run_stddev_test=False, c_years=False,
                 exclusion=False, all_and_excluded_analysis=False) -> None:

        n_files = 0

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
                    detcs.append(L(det, input_dir + det + ".csv", mixed_data=mixed_data,
                                   full_data_analysis=all_and_excluded_analysis))

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
                ratios12 = Ratios(detcs[0], detcs[1], year=years, c_years=c_years, lumi_type=lumi_type,
                                  compute_by_run_by_fill=setts.get_by_fill_by_run_plots)
            else:
                ratios12 = Ratios(detcs[0], detcs[1], lumi_type=lumi_type,
                                  compute_by_run_by_fill=setts.get_by_fill_by_run_plots)

            if run_stddev_test:
                # nls vs. Standard deviation test
                test = BNLS(detcs[0], detcs[1])

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
            if setts.get_extra_plots:
                ratios12.plot_nls_ratio_vs_run()
                ratios12.plot_nls_ratio_vs_fill()
                ratios12.plot_ratio_vs_run()
                ratios12.plot_ratio_vs_fill()
                ratios12.plot_bad_fills()
                ratios12.plot_bad_runs()

            # By fill, by run
            if setts.get_by_fill_by_run_plots:
                ratios12.plot_by_fill_ratio_vs_fill()
                ratios12.plot_by_run_ratio_vs_run()
                # ratios12.plot_by_run_ratio_vs_run_with_errors()

            # all/excluded analysis plots
            if all_and_excluded_analysis:
                ratios12.plot_all_and_excluded_by_detc()

            ratios12.save_plots()

        elif n_files == 3:
            print("******** 3 detector comparison choose!! ******")

            if c_years:
                years_and_dir = input_dir.split(',')
                years = years_and_dir[1]
                for i in range(2, len(years_and_dir), 1):
                    years = years + ',' + years_and_dir[i]
                ratios123 = MRatios(detcs[0], detcs[1], detcs[2], year=years, c_years=c_years, lumi_type=lumi_type)
            else:
                ratios123 = MRatios(detcs[0], detcs[1], detcs[2], lumi_type=lumi_type)

            # Fill plots
            ratios123.plot_nls_ratios_vs_date()
            ratios123.plot_ratios_vs_date()

            ratios123.plot_ratios_vs_lumi3()
            ratios123.plot_nls_ratios_vs_lumi3()

            # ratios123.plot_nls_ratios_vs_run()
            # ratios123.plot_nls_ratios_vs_fill()
            ratios123.plot_ratios_vs_run()
            # ratios123.plot_ratios_vs_fill()

            if all_and_excluded_analysis:
                ratios123.plot_all_and_excluded_vs_lumi2()

            # Save plots
            ratios123.save_plots()