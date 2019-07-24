from tools.luminometer import Luminometer as L
from tools.detectorsratio import DetectorsRatio as Ratios
from tools.detectorsratio import MultipleDetectorsRatio as MRatios
from tools.linearityanalysis import LinearityAnalysis
import settings as setts

class LAnalysis:
    def __init__(self, dets_file_labels: list, input_dir: str, vs_all_analysis = False,
                 run_linearity_analysis=False, mixed_data=False) -> None:

        n_files = len(dets_file_labels)
        detcs = []
        for det in dets_file_labels:
            try:
                detcs.append(L(det, input_dir + det + ".csv", mixed_data=mixed_data))
            except IOError as errIO:
                print(errIO)
                print('Please check if default input folder is correctly created: ' + setts.csv_input_base_dir)
                print('Also check that your .csv file is in the correct year folder: ' + input_dir)
                raise
        if n_files == 2:
            print("******** 2 detector comparison choose!! ******")

            ratios12 = Ratios(detcs[0], detcs[1])

            # Fill stability plots
            ratios12.plot_ratio_vs_date()
            ratios12.plot_ratio_vs_lumi2()
            ratios12.plot_nls_ratio_vs_lumi2()
            ratios12.plot_ratio_hist()
            ratios12.plot_nls_ratio_hist()
            ratios12.plot_ratio_hist_weighted()

            # Mine
            ratios12.plot_nls_ratio_vs_run()
            ratios12.plot_nls_ratio_vs_fill()
            ratios12.plot_ratio_vs_run()
            ratios12.plot_ratio_vs_fill()

            ratios12.save_plots()

            if run_linearity_analysis:
                # Linearity
                lin_analysis = LinearityAnalysis(ratios12)
                lin_analysis.plot_hist_sbil()
                lin_analysis.plot_ratio_vs_all_data_sbil()
                lin_analysis.save_plots()

        elif n_files == 3:
            print("******** 3 detector comparison choose!! ******")

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
