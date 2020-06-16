import pandas as pd
import settings_bpm as setts
import tools.plotting_tools as plotting
from scipy import interpolate
import numpy as np
import tools.lumi_tools as ltools


class BPM:
    __doros_name = "DOROS"
    __arcBPM_name = "arcBPM"
    __nominal_name = "Nominal"
    __allowed_detectors = (__nominal_name, __doros_name)
    __col_b1h = "B1_H"
    __col_b1v = "B1_V"
    __col_b2h = "B2_H"
    __col_b2v = "B2_V"
    __col_b1hR = "B1_H_R"
    __col_b1vR = "B1_V_R"
    __col_b2hR = "B2_H_R"
    __col_b2vR = "B2_V_R"
    __col_b1hL = "B1_H_L"
    __col_b1vL = "B1_V_L"
    __col_b2hL = "B2_H_L"
    __col_b2vL = "B2_V_L"

    __col_time = "timestamp"
    __col_time_min = "t [min]"

    __distance_unit = r'$\mu\,m$'

    __scale_position = {
        __doros_name: 1000.,
        __nominal_name: 1000.,
        __arcBPM_name: 1.
    }
    __scale_time = 0.016666667
    __cols_b1_names = (__col_b1h, __col_b1v)
    __cols_b2_names = (__col_b2h, __col_b2v)
    __cols_b1_names_LR = (__col_b1hL, __col_b1vL, __col_b1hR, __col_b1vR)
    __cols_b2_names_LR = (__col_b2hL, __col_b2vL, __col_b2hR, __col_b2vR)

    __cols_to_be_position_scaled = (__col_b1h, __col_b1v, __col_b2h, __col_b2v,
                                    __col_b1hR, __col_b1vR, __col_b2hR, __col_b2vR,
                                    __col_b1hL, __col_b1vL, __col_b2hL, __col_b2vL)

    __ref_col_names = (__col_b1h, __col_b1v, __col_b2h, __col_b2v)
    __ref_col_names_LR = (__col_b1hR, __col_b1vR, __col_b2hR, __col_b2vR,
                          __col_b1hL, __col_b1vL, __col_b2hL, __col_b2vL)
    __ref_col_diff_suffix_label = '_diff'
    __ref_col_lscale_suffix_label = '_lscale'
    __ref_col_deflection_suffix_label = '_deflection'
    __ref_col_orbit_drift_suffix_label = '_orbit_drift'

    __col_names_to_read = {
        __doros_name: {"sec": __col_time,
                       "LHC.BPM.1L5.B1_DOROS:POS_H": __col_b1hL,
                       "LHC.BPM.1L5.B1_DOROS:POS_V": __col_b1vL,
                       "LHC.BPM.1L5.B2_DOROS:POS_H": __col_b2hL,
                       "LHC.BPM.1L5.B2_DOROS:POS_V": __col_b2vL,
                       "LHC.BPM.1R5.B1_DOROS:POS_H": __col_b1hR,
                       "LHC.BPM.1R5.B1_DOROS:POS_V": __col_b1vR,
                       "LHC.BPM.1R5.B2_DOROS:POS_H": __col_b2hR,
                       "LHC.BPM.1R5.B2_DOROS:POS_V": __col_b2vR},
        __nominal_name: {"sec": __col_time,
                         "set_nominal_B1xingPlane": __col_b1h,
                         "set_nominal_B1sepPlane": __col_b1v,
                         "set_nominal_B2xingPlane": __col_b2h,
                         "set_nominal_B2sepPlane": __col_b2v}
    }

    __col_correction_ini_time = 'timestamp_start'
    __col_correction_end_time = 'timestamp_end'
    __col_correction_x = 'correction_x'
    __col_correction_y = 'correction_y'

    def __init__(self, name: str, fill: int, data_file_name: str = None,
                 get_only_data: bool = False, nominal_data=None) -> None:
        self.name = name

        if self.name not in BPM.__allowed_detectors:
            raise AssertionError("Detector " + self.name + " no implemented")
        if self.name not in list(BPM.__col_names_to_read):
            raise AssertionError("Detector " + self.name + " data reading setting not available")

        self.__output_dir = "plots_bpm/" + str(fill) + "/" + self.name + "/"
        self.__plt_plots = {}
        self.__sns_plots = {}
        self.__nominal_data = None
        self.__path_to_data = []
        self.__offsets = {}
        self.__min_timestamp = None
        self.__max_timestamp = None
        self.__do_deep_studies = False
        self.__offsets_time_split = []

        # column names
        self.__ref_col_lscale_names = []
        self.__ref_col_lscale_names_LR = []
        self.__ref_col_orbit_drift_names = []
        self.__ref_col_orbit_drift_names_LR = []
        self.__ref_col_diff_names = []
        self.__ref_col_diff_names_LR = []
        self.__ref_col_diff_names_b1_LR = []
        self.__ref_col_diff_names_b2_LR = []
        self.__ref_col_lscale_diff_names = []
        self.__ref_col_lscale_diff_names_LR = []
        self.__ref_col_lscale_diff_names_b1_LR = []
        self.__ref_col_lscale_diff_names_b2_LR = []
        self.__ref_col_deflection_diff_names = []
        self.__ref_col_deflection_diff_names_LR = []
        self.__ref_col_orbit_drift_diff_names = []
        self.__ref_col_orbit_drift_diff_names_LR = []
        self.__ref_col_lscale_and_deflection_diff_names = []
        self.__ref_col_lscale_and_deflection_diff_names_LR = []
        self.__ref_col_deflection_and_lscale_diff_names = []
        self.__ref_col_deflection_and_lscale_diff_names_LR = []
        self.__ref_col_lscale_names_b1_LR = []
        self.__ref_col_lscale_names_b2_LR = []
        self.__ref_col_deflection_diff_names_b1_LR = []
        self.__ref_col_deflection_diff_names_b2_LR = []
        self.__ref_col_lscale_and_deflection_diff_names_b1_LR = []
        self.__ref_col_lscale_and_deflection_diff_names_b2_LR = []
        self.__ref_col_deflection_and_lscale_diff_names_b1_LR = []
        self.__ref_col_deflection_and_lscale_diff_names_b2_LR = []
        self.__ref_col_deflection_and_orbit_and_lscale_diff_names = []
        self.__ref_col_deflection_and_orbit_and_lscale_diff_names_LR = []

        self.fill_col_names()

        if self.name == BPM.__nominal_name:
            self.__is_nominal = True
            self.__interpolation_to_nominal_time = None
        else:
            self.__is_nominal = False
            self.__nominal_det = nominal_data
            self.__nominal_data = nominal_data.data
            self.__interpolation_to_nominal_time = pd.DataFrame()

        self.__rename_col_dict = BPM.__col_names_to_read[self.name]
        self.__cols_to_read_from_file = list(self.__rename_col_dict)
        self.__cols_after_renaming = [self.__rename_col_dict[x] for x in self.__cols_to_read_from_file]

        if fill:
            self.__fill = fill
            self.__path_to_data.append(self.getDataPathFromFillNumber(fill))
        elif data_file_name:
            self.__path_to_data.append(data_file_name)

        self.__settings = setts.config_dict[self.name][fill]
        self.__settings_list = list(self.__settings)

        if setts.conf_label_deep_studies in self.__settings_list:
            if self.__settings[setts.conf_label_deep_studies]:
                self.__do_deep_studies = True
                print("Deep studies will be performed.")

        self.__y_range = None
        self.__y_diff_range = None

        if setts.conf_label_y_range in list(setts.config_dict[self.name][fill]):
            self.__y_range = setts.config_dict[self.name][fill][setts.conf_label_y_range]
        else:
            self.__y_range = [setts.l_min_plot, setts.l_max_plot]

        if setts.conf_label_y_diff_range in list(setts.config_dict[self.name][fill]):
            self.__y_diff_range = setts.config_dict[self.name][fill][setts.conf_label_y_diff_range]
        else:
            self.__y_diff_range = [setts.l_diff_zoom_min_plot, setts.l_diff_zoom_max_plot]

        if setts.conf_label_beam_deflection in self.__settings_list:
            self.__apply_deflection = True
        else:
            self.__apply_deflection = False

        if setts.conf_label_length_scale in self.__settings_list:
            self.__apply_lscale = True
        else:
            self.__apply_lscale = False

        if setts.conf_label_orbit_drift in self.__settings_list:
            self.__apply_orbit_drift = True
        else:
            self.__apply_orbit_drift = False

        print("Analysing " + self.name + " using data: " + str(self.__path_to_data))
        in_data = pd.DataFrame()
        n_file = 0
        for file_path in self.__path_to_data:
            if n_file == 0:
                in_data = pd.read_csv(file_path, usecols=self.__cols_to_read_from_file)
            else:
                in_data.append(pd.read_csv(file_path, usecols=self.__cols_to_read_from_file), ignore_index=True)
            n_file += 1

        in_data.rename(columns=self.__rename_col_dict, inplace=True)

        if setts.conf_label_min_time in list(setts.config_dict[self.name][fill]) \
                and setts.conf_label_max_time in list(setts.config_dict[self.name][fill]):
            self.__min_timestamp = setts.config_dict[self.name][fill][setts.conf_label_min_time]
            self.__max_timestamp = setts.config_dict[self.name][fill][setts.conf_label_max_time]

            in_data = in_data[
                (in_data[BPM.__col_time] >= self.__min_timestamp) &
                (in_data[BPM.__col_time] <= self.__max_timestamp)]

        in_data.sort_values(BPM.__col_time, inplace=True)
        in_data.dropna(inplace=True)
        in_data.reset_index(drop=True, inplace=True)
        self.__zero_time = in_data[BPM.__col_time][0]
        if setts.conf_label_zero_time in setts.config_dict[self.name][fill]:
            self.__zero_time = setts.config_dict[self.name][fill][setts.conf_label_zero_time]

        self.__in_data_def_format = self.covert_cols_to_default_format(in_data)

        self.__in_data_def_format = self.apply_offsets(self.__in_data_def_format)
        # print(self.__in_data_def_format)

        self.__base_cols = list(BPM.__ref_col_names)
        if not self.__is_nominal:
            self.get_interpolation_to_nominal_time()
            if self.__do_deep_studies:
                self.__base_cols.extend(list(BPM.__ref_col_names_LR))
            self.get_cols_diff(in_data=self.__interpolation_to_nominal_time, cols_names=self.__base_cols)
            self.__ref_col_diff_names_final_result = self.__ref_col_diff_names
            self.add_time_min_col(self.__interpolation_to_nominal_time)
        else:
            self.__interpolation_to_nominal_time = self.__in_data_def_format

        self.apply_corrections_to_diff(basic_cols=self.__base_cols,
                                       apply_lscale=self.__apply_lscale,
                                       apply_deflection=self.__apply_deflection,
                                       apply_orbit_drift=self.__apply_orbit_drift)

        # Plotting
        if not get_only_data:
            print("Final results columns: " + str(self.__ref_col_diff_names_final_result))
            self.plot_detector_data()
            if setts.conf_label_special_time_intervals in list(setts.config_dict[self.name][self.__fill]):
                for xrange in setts.config_dict[self.name][self.__fill][setts.conf_label_special_time_intervals]:
                    self.plot_detector_data(xrange=xrange)
            # if setts.conf_label_length_scale in self.__settings_list:
            #     self.plot_detector_data(extra_name_suffix="_scale_corr", cols_to_plot=self.__ref_col_lscale_names,
            #                             data_to_use=self.__interpolation_to_nominal_time)

            if self.name != BPM.__nominal_name:
                self.plot_detector_all_data()
                self.plot_detector_data_diff_with_nominal()
                if self.__do_deep_studies:
                    self.plot_detector_data_diff_with_nominal(cols_to_plot=self.__ref_col_diff_names_LR,
                                                              data_to_use=self.__interpolation_to_nominal_time,
                                                              extra_name_suffix="_LR")
                    self.plot_detector_data_diff_with_nominal(cols_to_plot=self.__ref_col_diff_names_b1_LR,
                                                              data_to_use=self.__interpolation_to_nominal_time,
                                                              extra_name_suffix="_b1_LR")
                    self.plot_detector_data_diff_with_nominal(cols_to_plot=self.__ref_col_diff_names_b2_LR,
                                                              data_to_use=self.__interpolation_to_nominal_time,
                                                              extra_name_suffix="_b2_LR")
                if setts.conf_label_special_time_intervals in list(setts.config_dict[self.name][self.__fill]):
                    for xrange in setts.config_dict[self.name][self.__fill][setts.conf_label_special_time_intervals]:
                        self.plot_detector_data_diff_with_nominal(xrange=xrange)
                        self.plot_detector_data_diff_with_nominal(extra_name_suffix="_final",
                                                                  cols_to_plot=self.__ref_col_diff_names_final_result,
                                                                  data_to_use=self.__interpolation_to_nominal_time,
                                                                  xrange=xrange,
                                                                  legend_labels=self.__ref_col_names)
                if self.__apply_lscale:
                    self.plot_detector_data_diff_with_nominal(extra_name_suffix="_scale_corr",
                                                              cols_to_plot=self.__ref_col_lscale_diff_names,
                                                              data_to_use=self.__interpolation_to_nominal_time)
                if self.__apply_deflection:
                    # print(self.__ref_col_deflection_diff_names)
                    # print(list(self.__interpolation_to_nominal_time))
                    self.plot_detector_data_diff_with_nominal(extra_name_suffix="_deflection_corr",
                                                              cols_to_plot=self.__ref_col_deflection_diff_names,
                                                              data_to_use=self.__interpolation_to_nominal_time)
                if self.__apply_orbit_drift:
                    self.plot_detector_data_diff_with_nominal(extra_name_suffix="_orbit_drift",
                                                              cols_to_plot=self.__ref_col_orbit_drift_diff_names,
                                                              data_to_use=self.__interpolation_to_nominal_time)

                # if self.__apply_deflection and self.__apply_lscale:
                #     self.plot_detector_data_diff_with_nominal(extra_name_suffix="_deflection_and_lscale_corr",
                #                                               cols_to_plot=self.__ref_col_deflection_and_lscale_diff_names,
                #                                               data_to_use=self.__interpolation_to_nominal_time)
                #     if self.__do_deep_studies:
                #         self.plot_detector_data_diff_with_nominal(extra_name_suffix="_deflection_and_lscale_corr_LR_b1",
                #                                                   cols_to_plot=self.__ref_col_deflection_and_lscale_diff_names_b1_LR,
                #                                                   data_to_use=self.__interpolation_to_nominal_time)
                #         self.plot_detector_data_diff_with_nominal(extra_name_suffix="_deflection_and_lscale_corr_LR_b2",
                #                                                   cols_to_plot=self.__ref_col_deflection_and_lscale_diff_names_b2_LR,
                #                                                   data_to_use=self.__interpolation_to_nominal_time)
                #         if setts.conf_label_special_time_intervals in self.__settings_list:
                #             short_name_deflection_and_lscale_diff_names_b1_LR = exclude_suffix_from_names(
                #                         self.__ref_col_deflection_and_lscale_diff_names_b1_LR,
                #                         suffix=BPM.__ref_col_deflection_suffix_label + BPM.__ref_col_lscale_suffix_label
                #                                + BPM.__ref_col_diff_suffix_label)
                #             short_name_deflection_and_lscale_diff_names_b2_LR = exclude_suffix_from_names(
                #                 self.__ref_col_deflection_and_lscale_diff_names_b2_LR,
                #                 suffix=BPM.__ref_col_deflection_suffix_label + BPM.__ref_col_lscale_suffix_label
                #                        + BPM.__ref_col_diff_suffix_label)
                #             for xrange in setts.config_dict[self.name][self.__fill][
                #                 setts.conf_label_special_time_intervals]:
                #                 self.plot_detector_data_diff_with_nominal(xrange=xrange,
                #                                                           cols_to_plot=self.__ref_col_diff_names_b1_LR,
                #                                                           data_to_use=self.__interpolation_to_nominal_time,
                #                                                           extra_name_suffix="_b1_LR")
                #                 self.plot_detector_data_diff_with_nominal(
                #                     extra_name_suffix="_deflection_and_lscale_corr_b1_LR",
                #                     cols_to_plot=self.__ref_col_deflection_and_lscale_diff_names_b1_LR,
                #                     data_to_use=self.__interpolation_to_nominal_time,
                #                     xrange=xrange,
                #                     legend_labels=short_name_deflection_and_lscale_diff_names_b1_LR)
                #
                #                 self.plot_detector_data_diff_with_nominal(
                #                     xrange=xrange,
                #                     cols_to_plot=self.__ref_col_diff_names_b2_LR,
                #                     data_to_use=self.__interpolation_to_nominal_time,
                #                     extra_name_suffix="_b2_LR")
                #                 self.plot_detector_data_diff_with_nominal(
                #                     extra_name_suffix="_deflection_and_lscale_corr_b2_LR",
                #                     cols_to_plot=self.__ref_col_deflection_and_lscale_diff_names_b2_LR,
                #                     data_to_use=self.__interpolation_to_nominal_time,
                #                     xrange=xrange,
                #                     legend_labels=short_name_deflection_and_lscale_diff_names_b2_LR)

            self.save_plots()

        #   Saving data to file
        cols_to_save = [BPM.__col_time, BPM.__col_time_min] + list(BPM.__ref_col_names) + self.__ref_col_diff_names
        if self.__is_nominal:
            data_to_save = self.__in_data_def_format
        else:
            data_to_save = self.__interpolation_to_nominal_time
        ltools.save_columns_from_pandas_to_file(data_to_save, cols_to_save, self.__output_dir +
                                                "plotted_data.csv")
        # ltools.save_columns_from_pandas_to_file(data_to_save, self.__output_dir +
        #                                         "plotted_data_all.csv")

    @property
    def name(self):
        return self.__name

    @name.setter
    def name(self, name):
        self.__name = name

    @property
    def data(self):
        return self.__in_data_def_format

    def fill_col_names(self):
        self.__ref_col_lscale_names = [s + BPM.__ref_col_lscale_suffix_label for s in BPM.__ref_col_names]
        self.__ref_col_lscale_names_LR = [s + BPM.__ref_col_lscale_suffix_label for s in BPM.__ref_col_names_LR]
        self.__ref_col_lscale_names_b1_LR = [s + BPM.__ref_col_lscale_suffix_label for s in BPM.__cols_b1_names_LR]
        self.__ref_col_lscale_names_b2_LR = [s + BPM.__ref_col_lscale_suffix_label for s in BPM.__cols_b2_names_LR]
        self.__ref_col_orbit_drift_names = [s + BPM.__ref_col_orbit_drift_suffix_label for s in BPM.__ref_col_names]
        self.__ref_col_orbit_drift_names_LR = [s + BPM.__ref_col_orbit_drift_suffix_label for s in BPM.__ref_col_names_LR]

        if self.name == BPM.__nominal_name:
            self.__ref_col_diff_names = []
            self.__ref_col_diff_names_LR = []
        else:
            self.__ref_col_diff_names = [s + BPM.__ref_col_diff_suffix_label for s in BPM.__ref_col_names]
            self.__ref_col_diff_names_LR = [s + BPM.__ref_col_diff_suffix_label for s in BPM.__ref_col_names_LR]
            self.__ref_col_diff_names_b1_LR = [s + BPM.__ref_col_diff_suffix_label for s in BPM.__cols_b1_names_LR]
            self.__ref_col_diff_names_b2_LR = [s + BPM.__ref_col_diff_suffix_label for s in BPM.__cols_b2_names_LR]

            self.__ref_col_lscale_diff_names = \
                [s + BPM.__ref_col_diff_suffix_label for s in self.__ref_col_lscale_names]
            self.__ref_col_lscale_diff_names_LR = \
                [s + BPM.__ref_col_diff_suffix_label for s in self.__ref_col_lscale_names_LR]

            self.__ref_col_lscale_diff_names_b1_LR = \
                [s + BPM.__ref_col_diff_suffix_label for s in self.__ref_col_lscale_names_b1_LR]
            self.__ref_col_lscale_diff_names_b2_LR = \
                [s + BPM.__ref_col_diff_suffix_label for s in self.__ref_col_lscale_names_b2_LR]

            self.__ref_col_deflection_diff_names = \
                [s + BPM.__ref_col_deflection_suffix_label for s in self.__ref_col_diff_names]
            self.__ref_col_deflection_diff_names_LR = \
                [s + BPM.__ref_col_deflection_suffix_label for s in self.__ref_col_diff_names_LR]
            self.__ref_col_deflection_diff_names_b1_LR = \
                [s + BPM.__ref_col_deflection_suffix_label for s in self.__ref_col_diff_names_b1_LR]
            self.__ref_col_deflection_diff_names_b2_LR = \
                [s + BPM.__ref_col_deflection_suffix_label for s in self.__ref_col_diff_names_b2_LR]

            self.__ref_col_orbit_drift_diff_names = \
                [s + BPM.__ref_col_orbit_drift_suffix_label for s in self.__ref_col_diff_names]
            self.__ref_col_orbit_drift_diff_names_LR = \
                [s + BPM.__ref_col_orbit_drift_suffix_label for s in self.__ref_col_diff_names_LR]

            self.__ref_col_deflection_and_lscale_diff_names = \
                [s + BPM.__ref_col_deflection_suffix_label + BPM.__ref_col_lscale_suffix_label +
                 BPM.__ref_col_diff_suffix_label for s in BPM.__ref_col_names]
            self.__ref_col_deflection_and_lscale_diff_names_LR = \
                [s + BPM.__ref_col_deflection_suffix_label + BPM.__ref_col_lscale_suffix_label +
                 BPM.__ref_col_diff_suffix_label for s in BPM.__ref_col_names_LR]
            self.__ref_col_deflection_and_lscale_diff_names_b1_LR = \
                [s + BPM.__ref_col_deflection_suffix_label + BPM.__ref_col_lscale_suffix_label +
                 BPM.__ref_col_diff_suffix_label for s in BPM.__cols_b1_names_LR]
            self.__ref_col_deflection_and_lscale_diff_names_b2_LR = \
                [s + BPM.__ref_col_deflection_suffix_label + BPM.__ref_col_lscale_suffix_label +
                 BPM.__ref_col_diff_suffix_label for s in BPM.__cols_b2_names_LR]

            self.__ref_col_deflection_and_orbit_and_lscale_diff_names = \
                [s + BPM.__ref_col_deflection_suffix_label + BPM.__ref_col_orbit_drift_suffix_label +
                 BPM.__ref_col_lscale_suffix_label +
                 BPM.__ref_col_diff_suffix_label for s in BPM.__ref_col_names]
            self.__ref_col_deflection_and_orbit_and_lscale_diff_names_LR = \
                [s + BPM.__ref_col_deflection_suffix_label + BPM.__ref_col_orbit_drift_suffix_label +
                 BPM.__ref_col_lscale_suffix_label +
                 BPM.__ref_col_diff_suffix_label for s in BPM.__ref_col_names_LR]

            self.__ref_col_diff_names_final_result = self.__ref_col_diff_names

    def getDataPathFromFillNumber(self, fill: int):
        path_to_fill_data = setts.data_base_folder + self.name + "/Fill" + str(fill) + "/"
        if self.name in list(setts.config_dict):
            data_path = path_to_fill_data + setts.config_dict[self.name][fill][setts.conf_label_data_file_path]
        else:
            raise AssertionError("detector not configured in settings")

        return data_path

    def covert_cols_to_default_format(self, in_data: pd.DataFrame, rescale: bool = True):
        data_in_unified_format = in_data
        fill = self.__fill

        if rescale:
            rescale_factor = 1.
            time0 = self.__zero_time
            for i_col in data_in_unified_format.columns:
                if i_col == BPM.__col_time:
                    data_in_unified_format[BPM.__col_time_min] = (data_in_unified_format[
                                                                      BPM.__col_time] - time0) * BPM.__scale_time
                elif i_col in BPM.__cols_to_be_position_scaled:
                    rescale_factor = BPM.__scale_position[self.name]
                    data_in_unified_format[i_col] = data_in_unified_format[i_col] * rescale_factor
                else:
                    print("     *** WARNING -->> column " + i_col + " not scaled!")

        if setts.conf_label_special_col_sign in self.__settings_list:
            for col_name in list(setts.config_dict[self.name][self.__fill][setts.conf_label_special_col_sign]):
                col_factor = setts.config_dict[self.name][self.__fill][setts.conf_label_special_col_sign][col_name]
                data_in_unified_format[col_name] = data_in_unified_format[col_name] * col_factor
                print(" --> Column " + col_name + " has been multiplied by " + str(col_factor))

        return data_in_unified_format

    def apply_offsets(self, in_data: pd.DataFrame):
        data_in_unified_format = in_data

        if self.name == BPM.__doros_name:
            offsets_computed = False
            offsets_for_LR_cols_available = False
            if setts.conf_label_compute_offsets in self.__settings_list:
                if self.__settings[setts.conf_label_compute_offsets]:
                    self.compute_offsets()
                    offsets_computed = True
            if setts.conf_label_offset_values in self.__settings_list:
                self.get_offsets_from_config()

            # print(self.__offsets_time_split)
            # print(self.__offsets)
            offsets_time = self.__offsets_time_split
            n_offsets_intervals = len(offsets_time)
            cols_with_offsets = list(self.__offsets)

            # check if LR cols offsets are available, and if not filling dumb values with 0.
            if all(item in cols_with_offsets for item in BPM.__ref_col_names_LR):
                offsets_for_LR_cols_available = True
                print("LR cols offsets values loaded")
            else:
                dumb_offset_list = [0.] * n_offsets_intervals
                print(dumb_offset_list)
                for col_name in BPM.__ref_col_names_LR:
                    self.__offsets[col_name] = dumb_offset_list

            # Check if arrays are ordered in time
            if not sorted(offsets_time) == offsets_time:
                raise AssertionError("Offsets values most be ordered in time")

            col_b1h = []
            col_b1v = []
            col_b2h = []
            col_b2v = []

            col_b1h_L = []
            col_b1v_L = []
            col_b2h_L = []
            col_b2v_L = []

            col_b1h_R = []
            col_b1v_R = []
            col_b2h_R = []
            col_b2v_R = []

            # read offsets:
            pos_in_offsets_time_array = 0
            print("Applying offsets for time: " + str(offsets_time[pos_in_offsets_time_array]))
            off_b1h = self.__offsets[BPM.__col_b1h][0]
            off_b1v = self.__offsets[BPM.__col_b1v][0]
            off_b2h = self.__offsets[BPM.__col_b2h][0]
            off_b2v = self.__offsets[BPM.__col_b2v][0]

            off_b1hL = self.__offsets[BPM.__col_b1hL][0]
            off_b1vL = self.__offsets[BPM.__col_b1vL][0]
            off_b2hL = self.__offsets[BPM.__col_b2hL][0]
            off_b2vL = self.__offsets[BPM.__col_b2vL][0]
            off_b1hR = self.__offsets[BPM.__col_b1hR][0]
            off_b1vR = self.__offsets[BPM.__col_b1vR][0]
            off_b2hR = self.__offsets[BPM.__col_b2hR][0]
            off_b2vR = self.__offsets[BPM.__col_b2vR][0]

            for i_iter in range(0, len(data_in_unified_format)):

                if pos_in_offsets_time_array < n_offsets_intervals - 1 and \
                        data_in_unified_format[BPM.__col_time][i_iter] >= offsets_time[pos_in_offsets_time_array + 1]:
                    pos_in_offsets_time_array += 1
                    print("Applying offsets for time: " + str(offsets_time[pos_in_offsets_time_array]))
                    off_b1h = self.__offsets[BPM.__col_b1h][pos_in_offsets_time_array]
                    off_b1v = self.__offsets[BPM.__col_b1v][pos_in_offsets_time_array]
                    off_b2h = self.__offsets[BPM.__col_b2h][pos_in_offsets_time_array]
                    off_b2v = self.__offsets[BPM.__col_b2v][pos_in_offsets_time_array]

                    off_b1hL = self.__offsets[BPM.__col_b1hL][pos_in_offsets_time_array]
                    off_b1vL = self.__offsets[BPM.__col_b1vL][pos_in_offsets_time_array]
                    off_b2hL = self.__offsets[BPM.__col_b2hL][pos_in_offsets_time_array]
                    off_b2vL = self.__offsets[BPM.__col_b2vL][pos_in_offsets_time_array]
                    off_b1hR = self.__offsets[BPM.__col_b1hR][pos_in_offsets_time_array]
                    off_b1vR = self.__offsets[BPM.__col_b1vR][pos_in_offsets_time_array]
                    off_b2hR = self.__offsets[BPM.__col_b2hR][pos_in_offsets_time_array]
                    off_b2vR = self.__offsets[BPM.__col_b2vR][pos_in_offsets_time_array]

                value_b1hL = data_in_unified_format[BPM.__col_b1hL][i_iter]
                value_b1hR = data_in_unified_format[BPM.__col_b1hR][i_iter]
                value_b1vL = data_in_unified_format[BPM.__col_b1vL][i_iter]
                value_b1vR = data_in_unified_format[BPM.__col_b1vR][i_iter]

                value_b2hL = data_in_unified_format[BPM.__col_b2hL][i_iter]
                value_b2hR = data_in_unified_format[BPM.__col_b2hR][i_iter]
                value_b2vL = data_in_unified_format[BPM.__col_b2vL][i_iter]
                value_b2vR = data_in_unified_format[BPM.__col_b2vR][i_iter]

                mean_value_LR_b1h = (value_b1hL + value_b1hR) / 2
                mean_value_LR_b1v = (value_b1vL + value_b1vR) / 2
                mean_value_LR_b2h = (value_b2hL + value_b2hR) / 2
                mean_value_LR_b2v = (value_b2vL + value_b2vR) / 2

                col_b1h.append(mean_value_LR_b1h - off_b1h)
                col_b1v.append(mean_value_LR_b1v - off_b1v)
                col_b2h.append(mean_value_LR_b2h - off_b2h)
                col_b2v.append(mean_value_LR_b2v - off_b2v)

                if self.__do_deep_studies:
                    col_b1h_L.append(value_b1hL - off_b1hL)
                    col_b1v_L.append(value_b1vL - off_b1vL)
                    col_b2h_L.append(value_b2hL - off_b2hL)
                    col_b2v_L.append(value_b2vL - off_b2vL)

                    col_b1h_R.append(value_b1hR - off_b1hR)
                    col_b1v_R.append(value_b1vR - off_b1vR)
                    col_b2h_R.append(value_b2hR - off_b2hR)
                    col_b2v_R.append(value_b2vR - off_b2vR)

            data_in_unified_format[BPM.__col_b1h] = np.array(col_b1h)
            data_in_unified_format[BPM.__col_b1v] = np.array(col_b1v)
            data_in_unified_format[BPM.__col_b2h] = np.array(col_b2h)
            data_in_unified_format[BPM.__col_b2v] = np.array(col_b2v)

            if self.__do_deep_studies:
                data_in_unified_format[BPM.__col_b1hL] = np.array(col_b1h_L)
                data_in_unified_format[BPM.__col_b1vL] = np.array(col_b1v_L)
                data_in_unified_format[BPM.__col_b2hL] = np.array(col_b2h_L)
                data_in_unified_format[BPM.__col_b2vL] = np.array(col_b2v_L)

                data_in_unified_format[BPM.__col_b1hR] = np.array(col_b1h_R)
                data_in_unified_format[BPM.__col_b1vR] = np.array(col_b1v_R)
                data_in_unified_format[BPM.__col_b2hR] = np.array(col_b2h_R)
                data_in_unified_format[BPM.__col_b2vR] = np.array(col_b2v_R)

        return data_in_unified_format

    # It gets the difference between in_data specified cols and Nominal. If no matching to nominal cols dict is provided
    # then a default dict is constructed from the input cols taken the same ones for Nominal.
    def get_cols_diff(self, in_data, cols_names: list, cols_names_nominal: dict = None):
        # print(cols_names)
        if cols_names_nominal is None:
            cols_names_nominal = {}
            for i_col in cols_names:
                cols_names_nominal[i_col] = get_base_name_from_DOROS_LR(i_col)
        # print(cols_names_nominal)
        assert len(cols_names) == len(cols_names_nominal)
        for col_name in cols_names:
            nominal_col_name = cols_names_nominal[col_name]
            in_data[col_name + BPM.__ref_col_diff_suffix_label] = in_data[col_name] - \
                                                                  self.__nominal_data[nominal_col_name]

    def get_interpolation_to_nominal_time(self):
        in_data = self.__in_data_def_format
        self.__interpolation_to_nominal_time = pd.DataFrame()
        x = self.__nominal_data[BPM.__col_time]
        self.__interpolation_to_nominal_time[BPM.__col_time] = x
        if self.__do_deep_studies:
            col_names = BPM.__cols_to_be_position_scaled
        else:
            col_names = BPM.__ref_col_names
        for col_name in col_names:
            tck = interpolate.splrep(in_data[BPM.__col_time], in_data[col_name], s=0)
            y = interpolate.splev(x, tck, der=0)
            self.__interpolation_to_nominal_time[col_name] = np.array(y)

    def add_time_min_col(self, in_df):
        if BPM.__col_time in in_df.columns:
            time0 = self.__zero_time
            in_df[BPM.__col_time_min] = (in_df[BPM.__col_time] - time0) * BPM.__scale_time
        else:
            raise AssertionError("No time column available")

    def apply_length_scale_correction(self, base_cols: list, data_to_use=None):
        print("Applying length scale correction ...")
        if data_to_use is None:
            data_to_use = self.__interpolation_to_nominal_time
        # print(list(data_to_use))
        x_correction = setts.config_dict[self.name][self.__fill][setts.conf_label_length_scale][0]
        y_correction = setts.config_dict[self.name][self.__fill][setts.conf_label_length_scale][1]
        h_correction = x_correction
        v_correction = y_correction

        corrected_cols_dict = {}

        for i_col_name in base_cols:
            corrected_cols_dict[i_col_name] = []

        for index in range(0, len(data_to_use)):
            for i_col_name in base_cols:
                if is_H(i_col_name):
                    correction = h_correction
                elif is_V(i_col_name):
                    correction = v_correction
                else:
                    raise AssertionError()
                corrected_cols_dict[i_col_name].append(data_to_use[i_col_name][index] * correction)

        for i_col_name in base_cols:
            data_to_use[i_col_name + BPM.__ref_col_lscale_suffix_label] = \
                np.array(corrected_cols_dict[i_col_name])

    def read_deflection_correction_input(self):
        reading_cols_dict = {
            '#timestamp_start': BPM.__col_correction_ini_time,
            'timestamp_end': BPM.__col_correction_end_time,
            'correction_x': BPM.__col_correction_x,
            'correction_y': BPM.__col_correction_y
        }
        unit_factor = 1000.
        deflection_correction_data = pd.read_csv(self.__settings[setts.conf_label_beam_deflection],
                                                 index_col=False, usecols=list(reading_cols_dict))
        deflection_correction_data.rename(inplace=True, columns=reading_cols_dict)
        deflection_correction_data[BPM.__col_correction_x] *= unit_factor
        deflection_correction_data[BPM.__col_correction_y] *= unit_factor
        # print(deflection_correction_data)

        deflection_correction_data.sort_values(BPM.__col_correction_ini_time, inplace=True)
        deflection_correction_data.dropna(inplace=True)
        deflection_correction_data.reset_index(drop=True, inplace=True)

        # read correction geometrical factor from settings and set to 1. if not found
        x_gfactor = 2.
        y_gfactor = 2.
        if setts.conf_label_beam_deflection_gfactor in self.__settings_list:
            gfactors = self.__settings[setts.conf_label_beam_deflection_gfactor]
            x_gfactor = gfactors[0]
            y_gfactor = gfactors[1]
        else:
            print("WARNING: No geometrical factor was provided. Taking 2.0 as an estimated value")

        correction_gfactors = [x_gfactor, y_gfactor]
        return deflection_correction_data, correction_gfactors

    def read_orbit_drift_correction_input(self):
        reading_cols_dict = {
            '#timestamp_start': BPM.__col_correction_ini_time,
            'timestamp_end': BPM.__col_correction_end_time,
            'correction_x': BPM.__col_correction_x,
            'correction_y': BPM.__col_correction_y
        }
        unit_factor = 1000.
        orbit_drift_correction_data = pd.read_csv(self.__settings[setts.conf_label_beam_deflection],
                                                  index_col=False, usecols=list(reading_cols_dict))
        orbit_drift_correction_data.rename(inplace=True, columns=reading_cols_dict)
        orbit_drift_correction_data[BPM.__col_correction_x] *= unit_factor
        orbit_drift_correction_data[BPM.__col_correction_y] *= unit_factor

        orbit_drift_correction_data.dropna(inplace=True)
        orbit_drift_correction_data.sort_values(by=[BPM.__col_correction_ini_time])
        orbit_drift_correction_data.reset_index(drop=True, inplace=True)
        return orbit_drift_correction_data

    def apply_per_time_range_addition_correction(self, basic_cols, mode, data_to_use=None):
        if mode == "orbit drift":
            print("Applying orbit drift correction ...")
            correction_data = self.read_orbit_drift_correction_input()
            beam_factor = 0.5
            correction_gfactors = [1, 1]
            default_correction = 0.0
            corrected_suffix = BPM.__ref_col_orbit_drift_suffix_label
        elif mode == "deflection":
            print("Applying beam-beam deflection correction ...")
            correction_data, correction_gfactors = self.read_deflection_correction_input()
            beam_factor = 0.5
            default_correction = 0.0
            corrected_suffix = BPM.__ref_col_deflection_suffix_label
        else:
            raise AssertionError("Mode not implemented!!")

        x_gfactor = correction_gfactors[0]
        y_gfactor = correction_gfactors[1]

        if data_to_use is None:
            data_to_use = self.__interpolation_to_nominal_time

        corrected_cols_dict = {}

        for i_col_name in basic_cols:
            corrected_cols_dict[i_col_name] = []

        correction_time_index = 0
        ini_time = correction_data[BPM.__col_correction_ini_time][0]
        end_time = correction_data[BPM.__col_correction_end_time][0]

        for data_index in range(0, len(data_to_use)):
            data_time = data_to_use[BPM.__col_time][data_index]
            available_correction = False
            # print(data_time, ini_time, end_time)
            if data_time < ini_time or data_time > end_time:
                # print("Looking for correct index")
                time_matched_index = correction_data[
                    (correction_data[BPM.__col_correction_ini_time] <= data_time) &
                    (correction_data[BPM.__col_correction_end_time] >= data_time)].index.tolist()
                if len(time_matched_index) >= 1:
                    correction_time_index = time_matched_index[0]
                    ini_time = correction_data[BPM.__col_correction_ini_time][correction_time_index]
                    end_time = correction_data[BPM.__col_correction_end_time][correction_time_index]
                    available_correction = True
                # print(time_matched_index, ini_time, end_time)
            else:
                available_correction = True

            for i_col in basic_cols:
                uncorrected_val = data_to_use[i_col][data_index]
                if not available_correction:
                    # print("Not correction found")
                    correction = default_correction
                else:
                    if is_H(i_col):
                        correction_direction_label = BPM.__col_correction_x
                        gfactor = x_gfactor
                    elif is_V(i_col):
                        correction_direction_label = BPM.__col_correction_y
                        gfactor = y_gfactor
                    else:
                        raise AssertionError()
                    correction = \
                        correction_data[correction_direction_label][correction_time_index] * beam_factor * gfactor
                    if is_B1(i_col):
                        correction *= -1

                corrected_cols_dict[i_col].append(uncorrected_val+correction)

        for i_col_name in basic_cols:
            data_to_use[i_col_name + corrected_suffix] = \
                np.array(corrected_cols_dict[i_col_name])

    def apply_corrections_to_diff(self, basic_cols: list, apply_lscale: bool = False, apply_deflection: bool = False,
                                  apply_orbit_drift: bool = False):
        data_to_use = self.__interpolation_to_nominal_time

        col_names_diff_raw_data_to_raw_nominal = [basename + BPM.__ref_col_diff_suffix_label for basename in basic_cols]
        col_names_diff_lscale_data_to_lscale_nominal = [basename + BPM.__ref_col_lscale_suffix_label +
                                                        BPM.__ref_col_diff_suffix_label for basename in basic_cols]
        col_names_lscale_raw_data = [basename + BPM.__ref_col_lscale_suffix_label for basename in basic_cols]
        col_names_deflection_raw_data = [basename + BPM.__ref_col_deflection_suffix_label for basename in basic_cols]
        col_names_deflection_raw_nominal = [basename + BPM.__ref_col_deflection_suffix_label
                                            for basename in BPM.__ref_col_names]
        col_names_deflection_orbit_drift_raw_data = [basename + BPM.__ref_col_orbit_drift_suffix_label
                                                     for basename in col_names_deflection_raw_data]

        corrected_cols_suffix = ""
        if apply_lscale:
            corrected_cols_suffix += BPM.__ref_col_lscale_suffix_label
        if apply_deflection:
            corrected_cols_suffix += BPM.__ref_col_deflection_suffix_label

        if apply_deflection and not self.__is_nominal:
            self.apply_per_time_range_addition_correction(mode="deflection",
                                                          basic_cols=col_names_diff_raw_data_to_raw_nominal,
                                                          data_to_use=data_to_use)
            # print(list(data_to_use))
            self.__ref_col_diff_names_final_result = self.__ref_col_deflection_diff_names

        if apply_deflection and self.__is_nominal:
            self.apply_per_time_range_addition_correction(mode="deflection",
                                                          basic_cols=BPM.__ref_col_names,
                                                          data_to_use=data_to_use)
            # print(data_to_use)

        if apply_lscale:
            self.apply_length_scale_correction(base_cols=basic_cols, data_to_use=data_to_use)
            # print(list(data_to_use))
            if not self.__is_nominal:
                self.get_cols_diff(in_data=self.__interpolation_to_nominal_time, cols_names=col_names_lscale_raw_data)
                # print(list(data_to_use))
            self.__ref_col_diff_names_final_result = self.__ref_col_lscale_diff_names

        if apply_orbit_drift:
            self.apply_per_time_range_addition_correction(mode="orbit drift",
                                                          basic_cols=basic_cols, data_to_use=data_to_use)
            if not self.__is_nominal:
                self.apply_per_time_range_addition_correction(mode="orbit drift",
                                                              basic_cols=col_names_diff_raw_data_to_raw_nominal,
                                                              data_to_use=data_to_use)

        if apply_deflection and apply_orbit_drift and self.__is_nominal:
            self.apply_per_time_range_addition_correction(mode="orbit drift",
                                                          basic_cols=col_names_deflection_raw_nominal,
                                                          data_to_use=data_to_use)

        if apply_deflection and apply_lscale and self.__is_nominal:
            self.apply_length_scale_correction(base_cols=col_names_deflection_raw_nominal,
                                               data_to_use=data_to_use)

        if apply_deflection and apply_lscale and not self.__is_nominal:
            # case 2: BXdiff = (BXdoros - deflection) * dorosLS - BXnominal * nominalLS
            # BXdoros - deflection
            print("Applying combined correction: BXdiff = (BXdoros - deflection - OD) * dorosLS - BXnominal * nominalLS")
            self.apply_per_time_range_addition_correction(mode="deflection",
                                                          basic_cols=basic_cols,
                                                          data_to_use=data_to_use)
            input_cols_for_lscale = col_names_deflection_raw_data
            if apply_orbit_drift:
                self.apply_per_time_range_addition_correction(mode="orbit drift",
                                                              basic_cols=col_names_deflection_raw_data,
                                                              data_to_use=data_to_use)
                input_cols_for_lscale = col_names_deflection_orbit_drift_raw_data
                self.__ref_col_diff_names_final_result = self.__ref_col_deflection_and_orbit_and_lscale_diff_names
            else:
                self.__ref_col_diff_names_final_result = self.__ref_col_deflection_and_lscale_diff_names

            # (BXdoros - deflection) * dorosLS
            self.apply_length_scale_correction(base_cols=input_cols_for_lscale, data_to_use=data_to_use)

            # create correct matched Nominal cols (BXnominal * nominalLS) with respective (BXdoros - deflection) * dorosLS
            col_names_final_correction_data = []
            col_names_final_nominal_data = []
            col_names_final_correction_dict_nominal = {}
            if apply_orbit_drift:
                final_data_suffix = BPM.__ref_col_deflection_suffix_label + \
                                    BPM.__ref_col_orbit_drift_suffix_label + BPM.__ref_col_lscale_suffix_label
                final_nominal_suffix = BPM.__ref_col_lscale_suffix_label
            else:
                final_data_suffix = BPM.__ref_col_deflection_suffix_label + BPM.__ref_col_lscale_suffix_label
                final_nominal_suffix = BPM.__ref_col_lscale_suffix_label

            for basename in basic_cols:
                corrected_col_name = basename + final_data_suffix
                nominal_name_lscale = get_base_name_from_DOROS_LR(basename) + final_nominal_suffix
                col_names_final_correction_data.append(corrected_col_name)
                col_names_final_correction_dict_nominal[corrected_col_name] = nominal_name_lscale
                col_names_final_nominal_data.append(nominal_name_lscale)

            # print(col_names_deflection_lscale_raw_data_dict_nominal)
            # BXdiff = (BXdoros - deflection) * dorosLS - BXnominal * nominalLS
            self.get_cols_diff(in_data=data_to_use,
                               cols_names=col_names_final_correction_data,
                               cols_names_nominal=col_names_final_correction_dict_nominal)

        # print(self.__interpolation_to_nominal_time)

    def compute_offsets(self):
        print("Computing offsets ...")
        diffs = {}
        offsets = {}
        data_to_use = self.__in_data_def_format
        when_b1_h_is_zero = self.__nominal_data[BPM.__col_b1h] == 0.
        when_b1_v_is_zero = self.__nominal_data[BPM.__col_b1v] == 0.
        when_b2_h_is_zero = self.__nominal_data[BPM.__col_b2h] == 0.
        when_b2_v_is_zero = self.__nominal_data[BPM.__col_b2v] == 0.

        nominal_in_zero: pd.DataFrame = self.__nominal_data[when_b1_h_is_zero & when_b1_v_is_zero &
                                                            when_b2_h_is_zero & when_b2_v_is_zero].copy()

        mask = data_to_use[BPM.__col_time].isin(nominal_in_zero[BPM.__col_time])
        data_for_offsets = data_to_use[mask].copy()

        if len(data_for_offsets) != len(nominal_in_zero):
            # reduce nominal_in_zero to the available timestamps from data_for_offsets
            mask = nominal_in_zero[BPM.__col_time].isin(data_for_offsets[BPM.__col_time])
            nominal_in_zero = nominal_in_zero[mask].copy()

        # print(len(data_for_offsets), len(nominal_in_zero))
        assert len(data_for_offsets) == len(nominal_in_zero)

        nominal_in_zero.reset_index(inplace=True)
        data_for_offsets.reset_index(inplace=True)

        if self.name == BPM.__doros_name:
            raw_col_names = BPM.__ref_col_names_LR
        else:
            raise AssertionError("Detector not implemented in compute_offsets")

        for col_name in raw_col_names:
            diffs[col_name] = []
            self.__offsets[col_name] = []

        for i_iter in range(0, len(data_for_offsets)):
            # check that time matches
            assert data_for_offsets[BPM.__col_time][i_iter] == nominal_in_zero[BPM.__col_time][i_iter]
            for col_name in raw_col_names:
                nominal_col_name = get_base_name_from_DOROS_LR(col_name)
                diffs[col_name].append(data_for_offsets[col_name][i_iter] - nominal_in_zero[nominal_col_name][i_iter])

        diffs_df = pd.DataFrame.from_dict(diffs)
        diffs_df[BPM.__col_time] = nominal_in_zero[BPM.__col_time]

        self.__plt_plots["offsets_histos"] = diffs_df.hist(bins=50)[0][0].get_figure()

        diffs_df_split_time_list = []
        if setts.conf_label_offset_time in self.__settings_list:
            self.__offsets_time_split = self.__settings[setts.conf_label_offset_time]
            for index in range(0, len(self.__offsets_time_split)):
                if index == len(self.__offsets_time_split) - 1:
                    diffs_df_split_time_list.append(
                        diffs_df[diffs_df[BPM.__col_time] >= self.__offsets_time_split[index]])
                else:
                    diffs_df_split_time_list.append(
                        diffs_df[(diffs_df[BPM.__col_time] >= self.__offsets_time_split[index])
                                 & (diffs_df[BPM.__col_time] < self.__offsets_time_split[index + 1])])
            time_id = 0
            for i_df in diffs_df_split_time_list:
                self.__plt_plots["offsets_histos_" + str(self.__offsets_time_split[time_id])] = i_df.hist(bins=50)[0][
                    0].get_figure()
                for col_name in raw_col_names:
                    self.__offsets[col_name].append(i_df[col_name].mean())
                time_id += 1
        else:
            diffs_df_split_time_list = [diffs_df]
            for col_name in raw_col_names:
                self.__offsets[col_name] = diffs_df[col_name].mean()

    def get_offsets_from_config(self):
        offsets_values = setts.config_dict[self.name][self.__fill][setts.conf_label_offset_values]
        self.__offsets_time_split = self.__settings[setts.conf_label_offset_time]
        col_names = BPM.__ref_col_names
        for col_name in col_names:
            self.__offsets[col_name] = []
            for index in range(0, len(self.__offsets_time_split)):
                self.__offsets[col_name].append(offsets_values[index][setts.col_pos_offset_array[col_name]])

    def plot_detector_all_data(self, save_data_to_file=True):
        plot_name = self.name + "_all_data"
        print("plotting " + plot_name + "...")
        cols_to_plot = self.__cols_after_renaming
        cols_to_plot.remove(BPM.__col_time)

        self.__plt_plots[plot_name] = plotting.scatter_plot_from_pandas_frame(self.__in_data_def_format,
                                                                              x_data_label=BPM.__col_time,
                                                                              y_data_label=cols_to_plot,
                                                                              xlabel="time",
                                                                              ylabel="beam position", ncol_legend=4,
                                                                              plot_style='-',
                                                                              label_cms_status=False
                                                                              ).get_figure()

    def plot_detector_data(self, xrange=None, save_data_to_file=True, extra_name_suffix="", cols_to_plot=None,
                           data_to_use=None):
        plot_name = self.name + "_data" + extra_name_suffix
        x_min = None
        x_max = None
        if xrange:
            if xrange[0] == 0:
                x_min = xrange[0] + 0.0001
            else:
                x_min = xrange[0]
            x_max = xrange[1]
            plot_name += "_" + str(xrange[0]) + "_" + str(x_max)
        print("plotting " + plot_name + "...")

        if cols_to_plot is None:
            cols_to_plot = BPM.__ref_col_names

        if data_to_use is None:
            data_to_use = self.__in_data_def_format

        self.__plt_plots[plot_name] = plotting.scatter_plot_from_pandas_frame(data_to_use,
                                                                              x_data_label=BPM.__col_time_min,
                                                                              y_data_label=list(cols_to_plot),
                                                                              ymin=self.__y_range[0],
                                                                              ymax=self.__y_range[1],
                                                                              xmin=x_min, xmax=x_max,
                                                                              xlabel="time [min]",
                                                                              ylabel="beam position", ncol_legend=4,
                                                                              plot_style='-',
                                                                              label_cms_status=False
                                                                              ).get_figure()

    def plot_detector_data_diff_with_nominal(self, xrange=None, save_data_to_file=True, extra_name_suffix="",
                                             cols_to_plot=None, data_to_use=None, legend_labels=None):
        plot_name = self.name + "_data_diff_with_nominal" + extra_name_suffix
        x_min = None
        x_max = None
        if xrange:
            if xrange[0] == 0:
                x_min = xrange[0] + 0.0001
            else:
                x_min = xrange[0]
            x_max = xrange[1]
            plot_name += "_" + str(xrange[0]) + "_" + str(x_max)
        print("plotting " + plot_name + "...")
        if cols_to_plot is None:
            cols_to_plot = self.__ref_col_diff_names
        if data_to_use is None:
            data_to_use = self.__interpolation_to_nominal_time

        self.__plt_plots[plot_name] = plotting.scatter_plot_from_pandas_frame(data_to_use,
                                                                              x_data_label=BPM.__col_time_min,
                                                                              y_data_label=list(cols_to_plot),
                                                                              ymin=self.__y_diff_range[0],
                                                                              ymax=self.__y_diff_range[1],
                                                                              xmin=x_min, xmax=x_max,
                                                                              xlabel="time [min]",
                                                                              ylabel=self.name + " - Nominal diff." + " [" +
                                                                                     BPM.__distance_unit + "]",
                                                                              ncol_legend=4,
                                                                              plot_style='o',
                                                                              label_cms_status=False,
                                                                              legend_labels=legend_labels
                                                                              ).get_figure()

    def save_plots(self):
        print('\n\n Saving plots:')
        plotting.save_plots(self.__plt_plots, self.__output_dir, save_pickle=setts.save_figures_as_pickle)
        plotting.save_plots(self.__sns_plots, self.__output_dir, save_pickle=setts.save_figures_as_pickle)


def get_base_name_from_DOROS_LR(full_name: str):
    if "_L" in full_name:
        return full_name.replace("_L", "")
    elif "_R" in full_name:
        return full_name.replace("_R", "")
    else:
        return full_name


def check_string(input_string: str, to_check: str):
    if to_check in input_string:
        return True
    else:
        return False


def is_H(full_name: str):
    return check_string(full_name, "_H")


def is_V(full_name: str):
    return check_string(full_name, "_V")


def is_B1(full_name: str):
    return check_string(full_name, "B1_")


def is_B2(full_name: str):
    return check_string(full_name, "B2_")


def exclude_suffix_from_names(names_list: list, suffix: str):
    new_names = []
    for name in names_list:
        new_names.append(name.replace(suffix, ""))
    return new_names
