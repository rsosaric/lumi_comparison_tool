from tools.luminometer import Luminometer as L
import pandas as pd
import numpy as np
import settings as setts
from tools import lumi_tools as ltools
import tools.plotting_tools as plotting
from tools.detectorsratio import DetectorsRatio


class MultipleDetectorsRatio:
    def __init__(self, detcs_input, lumi_type: str,
                 year: str = None, energy: str = None, c_years=False,
                 nls: int = None, fill_norm_ratios: bool = False) -> None:

        number_of_dets = len(detcs_input)
        last_det_id = number_of_dets - 1
        ref_detector: L = detcs_input[last_det_id]
        self.__lumi_unit = ref_detector.lumi_unit

        print("Taking " + ref_detector.name + " as reference detector \n")

        list_of_detcs_energies = []
        list_of_detcs_years = []
        list_of_detcs_lumi_unit = []
        self.__comparison_name = None

        for i_det in detcs_input:
            list_of_detcs_energies.append(i_det.energy)
            list_of_detcs_years.append(i_det.year)
            list_of_detcs_lumi_unit.append(i_det.lumi_unit)

            if self.__comparison_name:
                self.__comparison_name += ("-" + i_det.name)
            else:
                self.__comparison_name = i_det.name

        if year is None:
            if len(np.unique(list_of_detcs_years)) == 1:
                self.year = ref_detector.year
            else:
                raise ValueError("Comparing detectors from different years makes no sense :( -->> " +
                                 str(list_of_detcs_years))
        else:
            self.year = year

        if energy is None:
            if len(np.unique(list_of_detcs_energies)) == 1:
                self.energy = ref_detector.energy
            else:
                raise ValueError("Comparing detectors from different energies makes no sense :( -->> " +
                                 str(list_of_detcs_energies))
        else:
            self.energy = energy

        assert len(np.unique(list_of_detcs_lumi_unit)) == 1
        self.__lumi_unit = ref_detector.lumi_unit

        if c_years:
            self.__output_dir = setts.default_output_dir + year.replace(',', '-') + "_" + str(self.energy) + 'TeV/' \
                                + self.__comparison_name + '/'
            self.__year_energy_label = str(self.energy) + 'TeV(' + year + ')'
        else:
            self.__output_dir = setts.default_output_dir + str(self.year) + "_" + str(self.energy) + 'TeV/' \
                                + self.__comparison_name + '/'
            self.__year_energy_label = str(self.energy) + 'TeV(' + str(self.year) + ')'

        # Checking all_data_analysis_included for all detectors
        self.__all_data_analysis_included = True
        for det in detcs_input:
            if not det.all_data_analysis_included:
                self.__all_data_analysis_included = False

        self.__ratios = []
        self.__ratios_all = []

        # sorting detectors by data size (except for ref --> end)
        dets_size = []
        detcs_dict = {}
        for i_det in detcs_input:
            detcs_dict[i_det.name] = i_det
            if i_det.name != ref_detector.name:
                dets_size.append((i_det.name, len(i_det.data)))
        dets_size.sort(key=lambda x: x[1])
        self.__dets = []
        for i_dets_size in dets_size:
            self.__dets.append(detcs_dict[i_dets_size[0]])
        self.__dets.append(ref_detector)


        # TODO: implement full_data comparison
        for i in range(0, number_of_dets):
            for j in range(i + 1, number_of_dets):
                print('Setting ' + str(self.__dets[i].name) + '/' + str(self.__dets[j].name) + ' ...')
                self.__ratios.append(DetectorsRatio(self.__dets[i], self.__dets[j],
                                                    lumi_type=lumi_type, fill_norm_ratios=fill_norm_ratios))

        # for i in range(number_of_dets - 1, -1, -1):
        #     for j in range(i - 1, -1, -1):
        #         print(i,j)
        #         print('Setting ' + str(self.__dets[j].name) + '/' + str(self.__dets[i].name) + ' ...')
        #         self.__ratios.append(DetectorsRatio(self.__dets[j], self.__dets[i],
        #                                             lumi_type=lumi_type, nls=nls, fill_norm_ratios=fill_norm_ratios))

        # self.__ref_ratio = self.__ratios[len(self.__ratios) - 1]
        number_of_ratios = len(self.__ratios)
        self.__ref_ratio = self.__ratios[-1]
        print("\n Taking " + self.__ref_ratio.label_ratio + " as reference ratio \n")

        self.__nls = self.__ref_ratio.nls

        self.__nls_ratio_col_names = []
        self.__nls_norm_ratio_col_names = []
        self.__ratio_plotting_labels = []

        self.__combinations_to_ref_ratios_labels = []
        self.__combinations_to_ref_nls_ratios_labels = []
        self.__combinations_to_ref_norm_nls_ratios_labels = []

        for i_ratio in self.__ratios:
            self.__nls_ratio_col_names.append(i_ratio.by_nls_label_ratio)
            self.__nls_norm_ratio_col_names.append(i_ratio.by_nls_label_ratio_norm)
            self.__ratio_plotting_labels.append(i_ratio.label_ratio)

            # selecting only ratios with reference detector in denominator
            if ref_detector.name == i_ratio.det2.name:
                self.__combinations_to_ref_ratios_labels.append(i_ratio.label_ratio)
                self.__combinations_to_ref_nls_ratios_labels.append(i_ratio.by_nls_label_ratio)
                self.__combinations_to_ref_norm_nls_ratios_labels.append(i_ratio.by_nls_label_ratio_norm)

        self.__lumi3_col_name = self.__ref_ratio.accumulated_lumi2_label
        self.__ratio_col_names = self.__ratio_plotting_labels

        keys_for_merging = ['fill', 'run', 'ls', 'time', 'date']

        # if len(self.__dets) == 3:
        #     merge_tmp = pd.merge(self.__ratios[0].common_data_filtered,
        #                          self.__ratios[1].common_data_filtered,
        #                          on=keys_for_merging,
        #                          how='outer').merge(self.__ratios[2].common_data_filtered,
        #                                             on=keys_for_merging, how='outer')
        #     if self.__all_data_analysis_included:
        #         merge_all_tmp = pd.merge(self.__ratios[0].data_all, self.__ratios[1].data_all,
        #                                  on=keys_for_merging,
        #                                  how='outer').merge(self.__ratios[2].data_all,
        #                                                     on=keys_for_merging, how='outer')
        #         merge_all_tmp = merge_all_tmp.reset_index(drop=True)
        #         ltools.check_and_clean_after_merging(merge_all_tmp)
        #         self.__all_data = merge_all_tmp
        #
        #         self.__exclusion_info_labels = []
        #         for i_ratio in self.__ratios:
        #             self.__exclusion_info_labels.append(i_ratio.by_nls_binary_exclusion_info_label)
        # else:
        #     raise ValueError('Number of detectors not implemented yet :(')

        merge_tmp = pd.merge(self.__ratios[0].common_data_filtered,
                             self.__ratios[1].common_data_filtered,
                             on=keys_for_merging,
                             how='outer')

        for i in range(2, number_of_ratios):
            merge_tmp = pd.merge(merge_tmp,
                                 self.__ratios[i].common_data_filtered,
                                 on=keys_for_merging,
                                 how='outer')

        if self.__all_data_analysis_included:
            merge_all_tmp = pd.merge(self.__ratios[0].data_all, self.__ratios[1].data_all,
                                     on=keys_for_merging,
                                     how='outer')
            for i in range(2, number_of_ratios):
                merge_all_tmp = pd.merge(merge_all_tmp, self.__ratios[i].data_all,
                                         on=keys_for_merging,
                                         how='outer')
            merge_all_tmp = merge_all_tmp.reset_index(drop=True)
            ltools.check_and_clean_after_merging(merge_all_tmp)
            self.__all_data = merge_all_tmp

            self.__exclusion_info_labels = []
            for i_ratio in self.__ratios:
                self.__exclusion_info_labels.append(i_ratio.by_nls_binary_exclusion_info_label)

        merge_tmp = merge_tmp.reset_index(drop=True)
        ltools.check_and_clean_after_merging(merge_tmp)

        self.__combined_data = merge_tmp
        # print(list(merge_tmp))

        self.__plt_plots = {}
        self.__sns_plots = {}

    def plot_ratios_vs_date(self):
        print('creating ratios_vs_date plot ...')
        ratios_vs_date = plotting.scatter_plot_from_pandas_frame(data_frame=self.combined_data,
                                                                 y_data_label=self.__ratio_col_names,
                                                                 x_data_label='date',
                                                                 ylabel="ratios",
                                                                 ymin=setts.ratio_min, ymax=setts.ratio_max,
                                                                 energy_year_label=self.__year_energy_label,
                                                                 # legend_labels=self.__nls_ratio_plotting_labels
                                                                 )
        self.__plt_plots['ratios_vs_date'] = ratios_vs_date.get_figure()

    def plot_ratios_vs_lumi3(self):
        print("creating ratios_vs_lumi3 plot ...")
        ratios_vs_lumi3 = plotting.scatter_plot_from_pandas_frame(data_frame=self.combined_data,
                                                                  y_data_label=self.__ratio_col_names,
                                                                  x_data_label=self.__lumi3_col_name,
                                                                  ylabel="ratios",
                                                                  xlabel="Integrated luminosity [$" +
                                                                         self.__lumi_unit + "^{-1}$]",
                                                                  ymin=setts.ratio_min, ymax=setts.ratio_max,
                                                                  energy_year_label=self.__year_energy_label,
                                                                  # legend_labels=self.__nls_ratio_plotting_labels
                                                                  )
        self.__plt_plots['ratios_vs_lumi3'] = ratios_vs_lumi3.get_figure()

    def plot_nls_ratios_vs_date(self):
        print('creating by_nls_ratios_vs_date plot ...')
        ratios_vs_date = plotting.scatter_plot_from_pandas_frame(data_frame=self.combined_data,
                                                                 y_data_label=self.__nls_ratio_col_names,
                                                                 x_data_label='date',
                                                                 ylabel="ratios in " + str(self.__nls) + ' LS',
                                                                 ymin=setts.ratio_min, ymax=setts.ratio_max,
                                                                 energy_year_label=self.__year_energy_label,
                                                                 legend_labels=self.__ratio_plotting_labels)
        self.__plt_plots['by_nls_ratios_vs_date'] = ratios_vs_date.get_figure()

    def plot_nls_ratios_vs_date_only2(self):
        print('creating by_nls_ratios_vs_date_only2 plot ...')
        ratios_vs_date = plotting.scatter_plot_from_pandas_frame(data_frame=self.combined_data,
                                                                 y_data_label=
                                                                 self.__combinations_to_ref_nls_ratios_labels,
                                                                 x_data_label='date',
                                                                 ylabel="ratios in " + str(self.__nls) + ' LS',
                                                                 ymin=setts.ratio_min, ymax=setts.ratio_max,
                                                                 energy_year_label=self.__year_energy_label,
                                                                 legend_labels=
                                                                 self.__combinations_to_ref_ratios_labels)
        self.__plt_plots['by_nls_ratios_vs_date_only2'] = ratios_vs_date.get_figure()

    def plot_nls_ratios_vs_lumi3(self):
        print('creating nls_ratios_vs_lumi3 plot ...')
        nls_ratios_vs_lumi3 = plotting.scatter_plot_from_pandas_frame(data_frame=self.combined_data,
                                                                      y_data_label=self.__nls_ratio_col_names,
                                                                      x_data_label=self.__lumi3_col_name,
                                                                      ylabel="ratios in " + str(self.__nls) + ' LS',
                                                                      xlabel="Integrated luminosity [$" +
                                                                             self.__lumi_unit + "^{-1}$]",
                                                                      ymin=setts.ratio_min, ymax=setts.ratio_max,
                                                                      energy_year_label=self.__year_energy_label,
                                                                      legend_labels=self.__ratio_plotting_labels)
        self.__plt_plots['nls_ratios_vs_lumi3'] = nls_ratios_vs_lumi3.get_figure()

    def plot_nls_ratios_vs_lumi3_only2(self):
        print('creating nls_ratios_vs_lumi3_only2 plot ...')
        nls_ratios_vs_lumi3 = plotting.scatter_plot_from_pandas_frame(data_frame=self.combined_data,
                                                                      y_data_label=
                                                                      self.__combinations_to_ref_nls_ratios_labels,
                                                                      x_data_label=self.__lumi3_col_name,
                                                                      ylabel="ratios in " + str(self.__nls) + ' LS',
                                                                      xlabel="Integrated luminosity [$" +
                                                                             self.__lumi_unit + "^{-1}$]",
                                                                      ymin=setts.ratio_min, ymax=setts.ratio_max,
                                                                      energy_year_label=self.__year_energy_label,
                                                                      legend_labels=
                                                                      self.__combinations_to_ref_ratios_labels)
        self.__plt_plots['nls_ratios_vs_lumi3_only2'] = nls_ratios_vs_lumi3.get_figure()

    def plot_nls_ratios_vs_lumi3_only2_norm(self):
        print('creating nls_ratios_vs_lumi3_only2_norm plot ...')
        nls_ratios_vs_lumi3 = plotting.scatter_plot_from_pandas_frame(data_frame=self.combined_data,
                                                                      y_data_label=
                                                                      self.__combinations_to_ref_norm_nls_ratios_labels,
                                                                      x_data_label=self.__lumi3_col_name,
                                                                      ylabel="ratios in " + str(self.__nls) + ' LS',
                                                                      xlabel="Integrated luminosity [$" +
                                                                             self.__lumi_unit + "^{-1}$]",
                                                                      ymin=setts.ratio_min, ymax=setts.ratio_max,
                                                                      energy_year_label=self.__year_energy_label,
                                                                      legend_labels=
                                                                      self.__combinations_to_ref_ratios_labels)
        self.__plt_plots['normalized_nls_ratios_vs_lumi3_only2'] = nls_ratios_vs_lumi3.get_figure()

    def plot_ratios_vs_run(self):
        print('creating ratios_vs_run plot ...')
        ratios_vs_run = plotting.scatter_plot_from_pandas_frame(data_frame=self.combined_data,
                                                                y_data_label=self.__ratio_col_names,
                                                                x_data_label='run',
                                                                xlabel="run",
                                                                ylabel="ratios",
                                                                ymin=setts.ratio_min, ymax=setts.ratio_max,
                                                                energy_year_label=self.__year_energy_label,
                                                                # legend_labels=self.__nls_ratio_plotting_labels
                                                                )
        self.__plt_plots['ratios_vs_run'] = ratios_vs_run.get_figure()

    # all/excluded data plots
    def plot_all_and_excluded_vs_lumi2(self):
        print('creating all_and_excluded_vs_lumi2 plot ...')
        all_and_excluded_vs_lumi2 = plotting.snsplot_detector_all_and_excluded(self.__all_data,
                                                                               x_data_label=self.__lumi3_col_name,
                                                                               y_data_label=self.__nls_ratio_col_names,
                                                                               conditional_label=
                                                                               self.__exclusion_info_labels,
                                                                               # conditional_label_extra=self.det2.excluded_label,
                                                                               xlabel="Integrated luminosity [$" +
                                                                                      self.__lumi_unit + "^{-1}$]",
                                                                               ylabel="ratios in " + str(
                                                                                   self.__nls) + ' LS',
                                                                               ymin=setts.ratio_min,
                                                                               ymax=setts.ratio_max,
                                                                               energy_year_label=
                                                                               self.__year_energy_label)
        self.__sns_plots["all_and_excluded_vs_lumi2"] = all_and_excluded_vs_lumi2

    # TODO: create combined hists plot

    @property
    def combined_data(self):
        return self.__combined_data

    def save_plots(self):
        print('\n\n Saving plots:')
        plotting.save_plots(self.__plt_plots, self.__output_dir)
        plotting.save_plots(self.__sns_plots, self.__output_dir)
