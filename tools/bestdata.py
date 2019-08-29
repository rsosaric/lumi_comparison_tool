import numpy as np
import pandas as pd
import settings as setts
from tools.luminometer import Luminometer as L
from tools.detectorsratio import DetectorsRatio as Ratios
import tools.plotting_tools as plotting
from tools import lumi_tools as ltools
from tools import full_run_utilities as frutils

class BestDataAnalysis():
    def __init__(self, dets_file_labels: list, input_dir: str, c_years: bool = False) -> None:
        print('Executing physics data analysis ...\n')
        mixed_data = True

        # dict [int year][str label_ratio] = mean_vale
        self.__ratio_pair_mean_info = {}

        # Class variables
        self.__detector_ratio_label = "Ratios"
        self.__detector_pair_percent_dict = None

        n_files = len(dets_file_labels)
        detcs = []

        if c_years:
            years_and_dir = input_dir.split(',')
            self.__input_dir = years_and_dir[0]
            files_path = []
            for i in range(1, len(years_and_dir), 1):
                files_path.append(years_and_dir[0] + years_and_dir[i] + '/')

            for det in dets_file_labels:
                detcs.append(frutils.merge_lumies(det, files_path, mixed_data=mixed_data))
            years = years_and_dir[1]
            for i in range(2, len(years_and_dir), 1):
                years = years + ',' + years_and_dir[i]
            ratios12 = Ratios(detcs[0], detcs[1], year=years, c_years=c_years)
        else:
            self.__input_dir = input_dir
            if n_files != 2:
                raise AssertionError("BestDataAnalysis only implemented for 2 'detectors'")
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
        self.__year = ratios12.year
        self.__label_ratio_normalized = self.__ratios.label_ratio + "_normalized"
        self.__label_ratio_nls_normalized = self.__ratios.label_ratio + "_nls_normalized"
        self.fill_detector_ratio_label_column()
        self.normalized_detector_ratios_by_pair()
        self.__common_data_filtered_no_nan = self.__ratios.common_data_filtered.dropna()
        self.__common_data_filtered_no_nan_index = self.__common_data_filtered_no_nan.reset_index(drop=True)

        # -> Initializing stats vars
        # single ls ratios
        # lumi weighted
        self.__nls_ratios_lw_mean = None
        self.__nls_ratios_lw_stdv = None
        self.__nls_ratios_lw_stdv_dof_corr = None

        self.fill_stats()

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
        self.plot_nls_ratio_hist_weighted_normalized()

        ratios12.save_plots()
        self.save_plots()

        #print (self.__ratio_pair_mean_info)

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
        print ("Adding normalized data ...")
        data_to_use = self.__ratios.common_data_filtered

        temp_array = []
        temp_nls_array = []

        for index_data in range(0, len(data_to_use)):
            index_year = int(data_to_use["date"][index_data].year)
            index_label_ratio = data_to_use[self.__detector_ratio_label][index_data]
            ratio_mean_factor = self.get_detector_pair_mean(detector_pair_label=index_label_ratio, year=index_year)
            temp_array.append(data_to_use[self.__ratios.label_ratio][index_data]/ratio_mean_factor)

            if np.isnan(data_to_use[self.__ratios.by_nls_label_ratio][index_data]):
                temp_nls_array.append(data_to_use[self.__ratios.by_nls_label_ratio][index_data])
            else:
                temp_nls_array.append(data_to_use[self.__ratios.by_nls_label_ratio][index_data]/ratio_mean_factor)

        data_to_use[self.__label_ratio_normalized] = np.array(temp_array)
        data_to_use[self.__label_ratio_nls_normalized] = np.array(temp_nls_array)


    def __get_detector_pair_mean_value_from_stat_file(self, detector_pair_label: str, year: int):
        pair_mean_value = None
        # splitting detectors:
        dets_labels = detector_pair_label.split("/")
        det1_label = ltools.convert_detector_name(dets_labels[0]).lower()
        det2_label = ltools.convert_detector_name(dets_labels[1]).lower()

        base_year_path = setts.csv_input_base_dir + str(year)
        base_year_path_plots = setts.default_output_dir + str(year)

        stats_file_name = base_year_path_plots + '/' + det1_label + '-' + det2_label + '/stats.csv'

        if setts.clean_run:
            print ("Run clean option have been set. Producing the file for "+ str(stats_file_name) + "...")
            det1 = L(det1_label, base_year_path + '/' + det1_label + ".csv")
            det2 = L(det2_label, base_year_path + '/' + det2_label + ".csv")

            ratio12_temp = Ratios(det1, det2)
            pair_mean_value = ratio12_temp.nls_ratios_lw_mean
        else:
            try:
                stats_file_pd = pd.read_csv(stats_file_name)
                print ("\n reading mean from stat file for: [" + str(year) +"]" + detector_pair_label + "\n")
                pair_mean_value = stats_file_pd["nls_ratios_lw_mean"][0]
            except:
                print (" \n File " + str(stats_file_name) + " not found. Producing the file ... \n")
                det1 = L(det1_label, base_year_path + '/' + det1_label + ".csv")
                det2 = L(det2_label, base_year_path + '/' + det2_label + ".csv")

                ratio12_temp = Ratios(det1, det2)
                pair_mean_value = ratio12_temp.nls_ratios_lw_mean

        return pair_mean_value

    def get_detector_pair_mean(self, detector_pair_label: str, year: int):
        if ltools.keys_exists(self.__ratio_pair_mean_info, year, detector_pair_label):
            pair_mean_value = self.__ratio_pair_mean_info[year][detector_pair_label]
        else:
            pair_mean_value = self.__get_detector_pair_mean_value_from_stat_file(detector_pair_label=detector_pair_label, year=year)
            if not ltools.keys_exists(self.__ratio_pair_mean_info, year):
                self.__ratio_pair_mean_info[year] = {}
            self.__ratio_pair_mean_info[year][detector_pair_label] = pair_mean_value

        return pair_mean_value



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
                                                                       ylabel=self.__ratios.label_ratio + " ratios",
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
                                                                       ylabel=self.__ratios.label_ratio + " ratios",
                                                                       ymin=setts.ratio_min,
                                                                       ymax=setts.ratio_max,
                                                                       energy_year_label=self.__ratios.year_energy_label,
                                                                       leg_col=4)
        self.__sns_plots['normalized_vs_lumi2_by_pair'] = normalized_vs_lumi2_by_pair

    def plot_nls_ratio_hist_weighted_normalized(self):
        ratio_hist_lumi2_w = plotting.hist_from_pandas_frame(data_frame=self.__common_data_filtered_no_nan,
                                                             col_label=self.__label_ratio_nls_normalized,
                                                             nbins=setts.nbins,
                                                             xlabel=self.__ratios.label_ratio + " ratios in " + str(
                                                                 self.__ratios.nls) + ' LS',
                                                             ylabel="Integrated luminosity [$" +
                                                                    self.__ratios.lumi_unit + "^{-1}$]",
                                                             # title='Detectors Ratios Histogram (lumi weighted)',
                                                             xmin=setts.ratio_min, xmax=setts.ratio_max,
                                                             #mean=self.__nls_ratios_lw_mean,
                                                             stdv=self.__nls_ratios_lw_stdv_dof_corr,
                                                             energy_year_label=self.__ratios.year_energy_label,
                                                             weight_label=self.__ratios.by_nls_lumi_label)
        self.__sns_plots['nls_ratio_hist_lw_normalized'] = ratio_hist_lumi2_w[0][0].get_figure()

    def save_plots(self):
        print('\n\n Saving plots:')
        #plotting.save_plots(self.__plt_plots, self.__output_dir)
        plotting.save_plots(self.__sns_plots, self.__output_dir)

    @property
    def detector_pair_percent_dict(self):
        if self.__detector_pair_percent_dict is None:
            raise AssertionError("Variable not filled")
        return self.__detector_pair_percent_dict


    def fill_stats(self):
        data = self.__common_data_filtered_no_nan_index

        nls_lw_stats = ltools.get_w_stats(data[self.__label_ratio_nls_normalized],
                                           data[self.__ratios.by_nls_lumi_label],
                                           min_val=setts.ratio_min, max_val=setts.ratio_max)

        self.__nls_ratios_lw_mean = nls_lw_stats.mean
        self.__nls_ratios_lw_stdv = nls_lw_stats.std_mean
        self.__nls_ratios_lw_stdv_dof_corr = nls_lw_stats.std