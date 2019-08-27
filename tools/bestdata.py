import numpy as np
import pandas as pd
import settings as setts
from tools.luminometer import Luminometer as L
from tools.detectorsratio import DetectorsRatio as Ratios
import tools.plotting_tools as plotting

class BestDataAnalysis():
    def __init__(self, dets_file_labels: list, input_dir: str, c_years: bool = False) -> None:
        print('Executing physics data analysis ...\n')
        mixed_data = True

        # Class variables
        self.__detector_ratio_label = "Ratios"
        self.__detector_pair_percent_dict = None


        n_files = len(dets_file_labels)
        if n_files != 2:
            raise AssertionError("BestDataAnalysis only implemented for 2 'detectors'")

        detcs = []
        for det in dets_file_labels:
            try:
                detcs.append(L(det, input_dir + det + ".csv", mixed_data=mixed_data))

            except IOError as errIO:
                print(errIO)
                print('Please check if default input folder is correctly created: ' + setts.csv_input_base_dir)
                print('Also check that your .csv file is in the correct year folder: ' + input_dir)
                raise

        ratios12 = Ratios(detcs[0], detcs[1])
        self.__ratios = ratios12
        self.__label_ratio_normalized = self.__ratios.label_ratio + "_normalized"
        self.fill_detector_ratio_label_column()
        self.normalized_detector_ratios_by_pair()

        #TODO: detector usage percent in physics and physics_compre

        print ("\n Detector Pairs usage: \n")
        #print (pd.DataFrame(self.detector_pair_percent_dict, columns=["detector pair", "usage(%)"]))
        print (self.detector_pair_percent_dict)

        self.__output_dir = self.__ratios.output_dir
        self.__sns_plots = {}

        ratios12.plot_ratio_hist()
        ratios12.plot_nls_ratio_hist_weighted()
        ratios12.plot_nls_ratio_vs_lumi2()
        ratios12.plot_nls_ratio_vs_run()
        ratios12.plot_nls_ratio_vs_fill()

        ratios12.plot_bad_fills()

        self.plot_hist_detectors()
        self.plot_vs_lumi2_by_pair()
        self.plot_normalized_vs_lumi2_by_pair()
        self.plot_normalized_hist_detectors()

        ratios12.save_plots()
        self.save_plots()

    def fill_detector_ratio_label_column(self):
        data = self.__ratios.common_data_filtered
        data[self.__detector_ratio_label] = data[self.__ratios.det1.detector_name_label] + "/" + data[
            self.__ratios.det2.detector_name_label]

        # Filling detector pairs percents
        temp_array = np.unique(data[self.__detector_ratio_label], return_counts=True)
        keys = []
        percents = []
        sum = 0.0
        for key in temp_array[0]:
            keys.append(key)
        for count in temp_array[1]:
            percents.append(count)
            sum += count

        percents_np = np.true_divide(np.array(percents), sum / 100.)

        self.__detector_pair_percent_dict = {}
        for i in range(0, len(keys)):
            self.__detector_pair_percent_dict[keys[i]] = percents_np[i]

    def normalized_detector_ratios_by_pair(self):
        data_to_use = self.__ratios.common_data_filtered

        temp_array = []

        for index_data in range(0, len(data_to_use)):
            ratio_mean_factor = setts.normalization_factor[int(self.__ratios.year)][data_to_use[self.__detector_ratio_label][index_data]]
            temp_array.append(data_to_use[self.__ratios.label_ratio][index_data]/ratio_mean_factor)

        data_to_use[self.__label_ratio_normalized] = np.array(temp_array)

    def plot_hist_detectors(self):
        hist_detectors = plotting.snsplot_hist_all_and_excluded(self.__ratios.common_data_filtered,
                                                               x_data_label=self.__ratios.label_ratio,
                                                               conditional_label=self.__detector_ratio_label,
                                                               # conditional_label_extra=self.det2.excluded_label,
                                                               xlabel="Detector ratios",
                                                               ylabel="Counts",
                                                               xmin=setts.ratio_min,
                                                               xmax=setts.ratio_max,
                                                               bins=setts.nbins,
                                                               energy_year_label=self.__ratios.year_energy_label)
        self.__sns_plots['hist_detectors'] = hist_detectors

    def plot_normalized_hist_detectors(self):
        normalized_hist_detectors = plotting.snsplot_hist_all_and_excluded(self.__ratios.common_data_filtered,
                                                               x_data_label=self.__label_ratio_normalized,
                                                               conditional_label=self.__detector_ratio_label,
                                                               # conditional_label_extra=self.det2.excluded_label,
                                                               xlabel="Detector ratios",
                                                               ylabel="Counts",
                                                               xmin=setts.ratio_min,
                                                               xmax=setts.ratio_max,
                                                               bins=setts.nbins,
                                                               energy_year_label=self.__ratios.year_energy_label)
        self.__sns_plots['normalized_hist_detectors'] = normalized_hist_detectors

    def plot_vs_lumi2_by_pair(self):
        vs_lumi2_by_pair = plotting.snsplot_detector_all_and_excluded(self.__ratios.common_data_filtered,
                                                                       x_data_label=self.__ratios.accumulated_rec_lumi2_label,
                                                                       y_data_label=self.__ratios.label_ratio,
                                                                       conditional_label=self.__detector_ratio_label,
                                                                       # conditional_label_extra=self.det2.excluded_label,
                                                                       xlabel="Integrated luminosity [$" +
                                                                              self.__ratios.lumi_unit + "^{-1}$]",
                                                                       ylabel=self.__ratios.label_ratio + " ratios in "
                                                                              + str(
                                                                           self.__ratios.nls) + ' LS',
                                                                       ymin=setts.ratio_min,
                                                                       ymax=setts.ratio_max,
                                                                       energy_year_label=self.__ratios.year_energy_label,
                                                                       leg_col=4)
        self.__sns_plots['vs_lumi2_by_pair'] = vs_lumi2_by_pair

    def plot_normalized_vs_lumi2_by_pair(self):
        normalized_vs_lumi2_by_pair = plotting.snsplot_detector_all_and_excluded(self.__ratios.common_data_filtered,
                                                                       x_data_label=self.__ratios.accumulated_rec_lumi2_label,
                                                                       y_data_label=self.__label_ratio_normalized,
                                                                       conditional_label=self.__detector_ratio_label,
                                                                       # conditional_label_extra=self.det2.excluded_label,
                                                                       xlabel="Integrated luminosity [$" +
                                                                              self.__ratios.lumi_unit + "^{-1}$]",
                                                                       ylabel=self.__ratios.label_ratio + " ratios in "
                                                                              + str(
                                                                           self.__ratios.nls) + ' LS',
                                                                       ymin=setts.ratio_min,
                                                                       ymax=setts.ratio_max,
                                                                       energy_year_label=self.__ratios.year_energy_label,
                                                                       leg_col=4)
        self.__sns_plots['normalized_vs_lumi2_by_pair'] = normalized_vs_lumi2_by_pair

    def save_plots(self):
        print('\n\n Saving plots:')
        #plotting.save_plots(self.__plt_plots, self.__output_dir)
        plotting.save_plots(self.__sns_plots, self.__output_dir)

    @property
    def detector_pair_percent_dict(self):
        if self.__detector_pair_percent_dict is None:
            raise AssertionError("Variable not filled")
        return self.__detector_pair_percent_dict