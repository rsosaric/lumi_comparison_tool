import pandas as pd
import settings_bpm as setts
import tools.plotting_tools as plotting
from scipy import interpolate
import numpy as np
import tools.lumi_tools as ltools

bpm_col_b1h = "B1_H"
bpm_col_b1v = "B1_V"
bpm_col_b2h = "B2_H"
bpm_col_b2v = "B2_V"
bpm_col_b1hR = "B1_H_R"
bpm_col_b1vR = "B1_V_R"
bpm_col_b2hR = "B2_H_R"
bpm_col_b2vR = "B2_V_R"
bpm_col_b1hL = "B1_H_L"
bpm_col_b1vL = "B1_V_L"
bpm_col_b2hL = "B2_H_L"
bpm_col_b2vL = "B2_V_L"

bpm_col_time = "timestamp"
bpm_col_time_min = "t [min]"


class BPM:
    __doros_name = "DOROS"
    __arcBPM_name = "arcBPM"
    __nominal_name = "Nominal"
    __allowed_detectors = (__nominal_name, __doros_name, __arcBPM_name)
    __col_b1h = bpm_col_b1h
    __col_b1v = bpm_col_b1v
    __col_b2h = bpm_col_b2h
    __col_b2v = bpm_col_b2v
    __col_b1hR = bpm_col_b1hR
    __col_b1vR = bpm_col_b1vR
    __col_b2hR = bpm_col_b2hR
    __col_b2vR = bpm_col_b2vR
    __col_b1hL = bpm_col_b1hL
    __col_b1vL = bpm_col_b1vL
    __col_b2hL = bpm_col_b2hL
    __col_b2vL = bpm_col_b2vL
    __col_H_diff = "H_diff"
    __col_V_diff = "V_diff"

    __col_time = bpm_col_time
    __col_time_min = bpm_col_time_min

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
    __ref_col_diff_suffix_label = '_diff_'
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
        __arcBPM_name: {"Timestamp": __col_time,
                        "B1_H_L": __col_b1hL,
                        "B1_V_L": __col_b1vL,
                        "B2_H_L": __col_b2hL,
                        "B2_V_L": __col_b2vL,
                        "B1_H_R": __col_b1hR,
                        "B1_V_R": __col_b1vR,
                        "B2_H_R": __col_b2hR,
                        "B2_V_R": __col_b2vR},
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
                 get_only_data: bool = False, apply_all_available_corrections=True, compute_orbit_drift=True,
                 nominal_data=None) -> None:
        self.name = name

        if self.name not in BPM.__allowed_detectors:
            raise AssertionError("Detector " + self.name + " no implemented")
        if self.name not in list(BPM.__col_names_to_read):
            raise AssertionError("Detector " + self.name + " data reading setting not available")

        self.__output_dir = "plots_bpm/" + str(fill) + "/" + self.name + "/"
        self.__plt_plots = {}
        self.__sns_plots = {}
        self.__nominal_data = None
        self.__data_in_zero_beam_position = None
        self.__path_to_data = []
        self.__offsets = {}
        self.__min_timestamp = None
        self.__max_timestamp = None
        self.__do_deep_studies = False
        self.__offsets_time_split = []
        self.__scans_info_dict = {}
        self.__scans_info_for_plotting = {}

        if self.name == BPM.__nominal_name:
            self.__is_nominal = True
            self.__interpolation_to_nominal_time = None
        else:
            self.__is_nominal = False
            self.__nominal_det = nominal_data
            self.__nominal_data = nominal_data.data
            self.__lscale_applied_in_nominal = self.__nominal_det.apply_lscale
            self.__orbit_drift_applied_in_nominal = self.__nominal_det.apply_orbit_drift
            self.__deflection_applied_in_nominal = self.__nominal_det.apply_deflection

            self.__interpolation_to_nominal_time = pd.DataFrame()

        # column names
        self.__ref_col_lscale_names = []
        self.__ref_col_lscale_names_LRAnalysing = []
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
        self.__ref_col_orbit_drift_and_lscale_diff_names = []
        self.__ref_col_orbit_drift_and_lscale_diff_names_LR = []
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

        self.__settings = setts.config_dict[self.name][fill]
        self.__settings_list = list(self.__settings)

        if fill:
            self.__fill = fill
            self.__path_to_data.append(self.getDataPathFromFillNumber(fill))
        elif data_file_name:
            self.__path_to_data.append(data_file_name)

        if setts.conf_label_deep_studies in self.__settings_list:
            if self.__settings[setts.conf_label_deep_studies]:
                self.__do_deep_studies = True
                print("Deep studies will be performed.")

        self.__y_range = None
        self.__y_diff_range = None

        if all(item in list(setts.config_dict[BPM.__nominal_name][self.__fill]) for item in
                [setts.conf_label_scans_time_windows, setts.conf_label_scans_names]):
            self.__scan_info_available = True
        else:
            self.__scan_info_available = False

        self.__y_range_orbit_drift = [setts.l_diff_orbit_drift_min_plot, setts.l_diff_orbit_drift_max_plot]

        if setts.conf_label_y_range in list(setts.config_dict[self.name][fill]):
            self.__y_range = setts.config_dict[self.name][fill][setts.conf_label_y_range]
        else:
            self.__y_range = [setts.l_min_plot, setts.l_max_plot]

        if setts.conf_label_y_diff_range in list(setts.config_dict[self.name][fill]):
            self.__y_diff_range = setts.config_dict[self.name][fill][setts.conf_label_y_diff_range]
        else:
            self.__y_diff_range = [setts.l_diff_zoom_min_plot, setts.l_diff_zoom_max_plot]

        if setts.conf_label_y_diff_smaller_range in list(setts.config_dict[self.name][fill]):
            self.__y_diff_smaller_range = setts.config_dict[self.name][fill][setts.conf_label_y_diff_smaller_range]
        else:
            self.__y_diff_smaller_range = self.__y_diff_range

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

        if setts.conf_label_exclusion_times in list(setts.config_dict[BPM.__nominal_name][self.__fill]):
            self.__times_to_exclude = \
                setts.config_dict[BPM.__nominal_name][self.__fill][setts.conf_label_exclusion_times]
        else:
            self.__times_to_exclude = []

        ltools.color_print("\n\nAnalysing " + self.name + " using data: " + str(self.__path_to_data), "green")

        if setts.conf_label_special_input_col_names in self.__settings_list:
            self.__rename_col_dict = self.__settings[setts.conf_label_special_input_col_names]
            print("     --->>> Using special input col names defined in settings. Please check if this is intended!!!")
        else:
            self.__rename_col_dict = BPM.__col_names_to_read[self.name]
        self.__cols_to_read_from_file = list(self.__rename_col_dict)
        self.__cols_after_renaming = [self.__rename_col_dict[x] for x in self.__cols_to_read_from_file]

        in_data = pd.DataFrame()
        n_file = 0
        for file_path in self.__path_to_data:
            if n_file == 0:
                in_data = pd.read_csv(file_path, usecols=self.__cols_to_read_from_file)
            else:
                in_data.append(pd.read_csv(file_path, usecols=self.__cols_to_read_from_file), ignore_index=True)
            n_file += 1

        in_data.rename(columns=self.__rename_col_dict, inplace=True)

        if setts.conf_label_timestamp_in_ms in self.__settings_list:
            if self.__settings[setts.conf_label_timestamp_in_ms]:
                in_data[BPM.__col_time] = in_data[BPM.__col_time] / 1000.

        in_data = in_data.astype({BPM.__col_time: int})

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

        if not self.__is_nominal:
            assert self.__zero_time == nominal_data.zero_time

        self.__in_data_def_format = self.covert_cols_to_default_format(in_data)

        self.__in_data_def_format = self.compute_extra_raw_avg_cols(self.__in_data_def_format)

        self.__in_data_def_format = self.apply_offsets(self.__in_data_def_format)

        self.__base_cols = list(BPM.__ref_col_names)

        if self.__scan_info_available:
            self.get_scans_data()

        if not self.__is_nominal:
            # (Method 1) use interpolated data to nominal times, maybe useful for data with low number of points
            # self.get_interpolation_to_nominal_time()
            # self.__processed_data = self.__interpolation_to_nominal_time

            # (Method 2) reduce Nominal and BPM to common timestamps
            common_available_times = np.intersect1d(self.__in_data_def_format[BPM.__col_time],
                                                    self.__nominal_data[BPM.__col_time])
            common_mask_nominal = self.__nominal_data[BPM.__col_time].isin(common_available_times)
            common_mask_bpm = self.__in_data_def_format[BPM.__col_time].isin(common_available_times)

            self.__nominal_data = self.__nominal_data[common_mask_nominal].copy()
            self.__processed_data = self.__in_data_def_format[common_mask_bpm].copy()

            self.__nominal_data.sort_values(BPM.__col_time, inplace=True)
            self.__processed_data.sort_values(BPM.__col_time, inplace=True)
            self.__nominal_data.reset_index(drop=True, inplace=True)
            self.__processed_data.reset_index(drop=True, inplace=True)

            if self.__do_deep_studies:
                self.__base_cols.extend(list(BPM.__ref_col_names_LR))
            self.get_cols_diff(in_data=self.__processed_data,
                               base_col_names=self.__base_cols)
            self.__ref_col_diff_names_final_result = self.__ref_col_diff_names
            self.add_time_min_col(self.__processed_data)

            if compute_orbit_drift:
                self.get_diffs_when_Nominal_is_zero()
                if not get_only_data:
                    self.plot_detector_data_diff_with_nominal(extra_name_suffix="_H_V_computed_orbit_drift",
                                                              cols_to_plot=[BPM.__col_H_diff, BPM.__col_V_diff],
                                                              data_to_use=self.__data_in_zero_beam_position,
                                                              yrange=self.__y_range_orbit_drift, marker_size=2.0,
                                                              leg_marker_scale=3,
                                                              plot_scan_info=True,
                                                              plot_scan_limits_lines=True,
                                                              scans_limits_to_use="time_range_minutes",
                                                              legend_labels=["H", "V"],
                                                              y_label="B2 - B1" +
                                                                      " [" + BPM.__distance_unit + "]")
                    self.plot_detector_data_diff_with_nominal(extra_name_suffix="_H_V_computed_orbit_drift_mod_limits",
                                                              cols_to_plot=[BPM.__col_H_diff, BPM.__col_V_diff],
                                                              data_to_use=self.__data_in_zero_beam_position,
                                                              yrange=self.__y_range_orbit_drift, marker_size=2.0,
                                                              leg_marker_scale=3,
                                                              plot_scan_info=True,
                                                              plot_scan_limits_lines=True,
                                                              scans_limits_to_use="mod_time_range_minutes",
                                                              legend_labels=["H", "V"],
                                                              y_label="B2 - B1" +
                                                                      " [" + BPM.__distance_unit + "]")
        else:
            self.__interpolation_to_nominal_time = self.__in_data_def_format
            self.__processed_data = self.__in_data_def_format

        if apply_all_available_corrections:
            self.apply_corrections(basic_cols=self.__base_cols,
                                   apply_lscale=self.__apply_lscale,
                                   apply_deflection=self.__apply_deflection,
                                   apply_orbit_drift=self.__apply_orbit_drift)

        # Plotting
        marker_size_for_sq_canvas_plots = 2.0
        use_squared_shape_for_final_plots = True
        if not get_only_data:
            if not self.__is_nominal:
                ltools.color_print(" \n\nFinal results columns: ", "green")
                ltools.color_print("    " + str(self.__ref_col_diff_names_final_result), "blue")
            self.plot_detector_data(plot_scan_info=self.__scan_info_available,
                                    plot_scan_limits_lines=self.__scan_info_available)

            if setts.conf_label_special_time_intervals in list(setts.config_dict[self.name][self.__fill]):
                for xrange in setts.config_dict[self.name][self.__fill][setts.conf_label_special_time_intervals]:
                    self.plot_detector_data(xrange=xrange)

            if self.name != BPM.__nominal_name:
                self.get_final_V_H_diff_cols()
                self.plot_detector_all_data()
                self.plot_detector_data_diff_with_nominal()
                self.plot_detector_data_diff_with_nominal(extra_name_suffix="_final",
                                                          cols_to_plot=self.__ref_col_diff_names_final_result,
                                                          data_to_use=self.__processed_data,
                                                          legend_labels=self.__ref_col_names)
                if self.__do_deep_studies:
                    self.plot_detector_data_diff_with_nominal(cols_to_plot=self.__ref_col_diff_names_LR,
                                                              data_to_use=self.__processed_data,
                                                              extra_name_suffix="_LR")
                    self.plot_detector_data_diff_with_nominal(cols_to_plot=self.__ref_col_diff_names_b1_LR,
                                                              data_to_use=self.__processed_data,
                                                              extra_name_suffix="_b1_LR")
                    self.plot_detector_data_diff_with_nominal(cols_to_plot=self.__ref_col_diff_names_b2_LR,
                                                              data_to_use=self.__processed_data,
                                                              extra_name_suffix="_b2_LR")
                if setts.conf_label_special_time_intervals in list(setts.config_dict[self.name][self.__fill]):
                    for xrange in setts.config_dict[self.name][self.__fill][setts.conf_label_special_time_intervals]:
                        self.plot_detector_data_diff_with_nominal(xrange=xrange,
                                                                  canvas_square_shape=use_squared_shape_for_final_plots,
                                                                  marker_size=marker_size_for_sq_canvas_plots,
                                                                  legend_labels=self.__ref_col_names)
                        self.plot_detector_data_diff_with_nominal(extra_name_suffix="_final",
                                                                  cols_to_plot=self.__ref_col_diff_names_final_result,
                                                                  data_to_use=self.__processed_data,
                                                                  xrange=xrange,
                                                                  marker_size=marker_size_for_sq_canvas_plots,
                                                                  canvas_square_shape=use_squared_shape_for_final_plots,
                                                                  legend_labels=self.__ref_col_names)

                if apply_all_available_corrections:
                    if self.__apply_lscale or self.__lscale_applied_in_nominal:
                        self.plot_detector_data_diff_with_nominal(extra_name_suffix="_scale_corr",
                                                                  cols_to_plot=self.__ref_col_lscale_diff_names,
                                                                  data_to_use=self.__processed_data,
                                                                  legend_labels=self.__ref_col_names)
                    if self.__apply_deflection or self.__deflection_applied_in_nominal:
                        # print(self.__ref_col_deflection_diff_names)
                        # print(list(self.__processed_data))
                        self.plot_detector_data_diff_with_nominal(extra_name_suffix="_deflection_corr",
                                                                  cols_to_plot=self.__ref_col_deflection_diff_names,
                                                                  data_to_use=self.__processed_data,
                                                                  legend_labels=self.__ref_col_names)
                        if self.apply_lscale:
                            self.plot_detector_data_diff_with_nominal(extra_name_suffix="_deflection_lscale_corr",
                                                                      cols_to_plot=self.__ref_col_deflection_and_lscale_diff_names,
                                                                      data_to_use=self.__processed_data,
                                                                      legend_labels=self.__ref_col_names)

                    if self.__apply_orbit_drift:
                        self.plot_detector_data_diff_with_nominal(extra_name_suffix="_orbit_drift",
                                                                  cols_to_plot=self.__ref_col_orbit_drift_diff_names,
                                                                  data_to_use=self.__processed_data,
                                                                  legend_labels=self.__ref_col_names)
                        if self.apply_lscale:
                            self.plot_detector_data_diff_with_nominal(extra_name_suffix="_orbit_drift_lscale",
                                                                      cols_to_plot=self.__ref_col_orbit_drift_and_lscale_diff_names,
                                                                      data_to_use=self.__processed_data,
                                                                      legend_labels=self.__ref_col_names)

                if len(self.__ref_col_diff_names_final_result) > 0:
                    if setts.conf_label_special_time_intervals in self.__settings_list:
                        for xrange in setts.config_dict[self.name][self.__fill][
                            setts.conf_label_special_time_intervals]:
                            self.plot_detector_data_diff_with_nominal(extra_name_suffix="_H_V_final",
                                                                      cols_to_plot=[BPM.__col_H_diff, BPM.__col_V_diff],
                                                                      data_to_use=self.__processed_data,
                                                                      legend_labels=["H", "V"],
                                                                      marker_size=marker_size_for_sq_canvas_plots,
                                                                      xrange=xrange,
                                                                      use_smaller_y_range=True,
                                                                      canvas_square_shape=use_squared_shape_for_final_plots,
                                                                      y_label="B2 - B1" +
                                                                              " [" + BPM.__distance_unit + "]")

                    self.plot_detector_data_diff_with_nominal(extra_name_suffix="_H_V_final",
                                                              cols_to_plot=[BPM.__col_H_diff, BPM.__col_V_diff],
                                                              data_to_use=self.__processed_data,
                                                              legend_labels=["H", "V"],
                                                              use_smaller_y_range=True,
                                                              y_label="B2 - B1" +
                                                                      " [" + BPM.__distance_unit + "]")

            self.save_plots()

            if not self.__is_nominal:
                #   Saving data to file
                data_to_save = self.__processed_data
                timing_cols = [BPM.__col_time, BPM.__col_time_min]
                data_cols = self.__ref_col_diff_names_final_result
                all_cols_to_save = timing_cols + data_cols
                final_hysteresis_cols = timing_cols + [BPM.__col_H_diff, BPM.__col_V_diff]

                # make sure to store data only in studied range
                ltools.color_print("\n\nSaving data in " + self.__output_dir + "plotted_data.csv", "green")
                ltools.color_print("Saving hysteresis data in " + self.__output_dir + str(self.__fill) +
                                   "_hysteresis.csv", "green")
                ltools.color_print("    Data saved only in range " + str(self.__y_diff_range), "blue")
                for i_col in data_cols:
                    data_to_save = data_to_save[(data_to_save[i_col] >= self.__y_diff_range[0]) &
                                                (data_to_save[i_col] <= self.__y_diff_range[1])]
                ltools.save_columns_from_pandas_to_file(data_to_save, all_cols_to_save, self.__output_dir +
                                                        "plotted_data.csv")
                ltools.save_columns_from_pandas_to_file(data_to_save, final_hysteresis_cols, self.__output_dir +
                                                        str(self.__fill) + "_hysteresis.csv")

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
        self.__ref_col_lscale_names = self.get_col_names(correction_suffix=BPM.__ref_col_lscale_suffix_label)
        self.__ref_col_lscale_names_LR = self.get_col_names(correction_suffix=BPM.__ref_col_lscale_suffix_label,
                                                            include_LR=True)
        self.__ref_col_lscale_names_b1_LR = self.get_col_names(correction_suffix=BPM.__ref_col_lscale_suffix_label,
                                                               include_LR=True, only_b1=True)
        self.__ref_col_lscale_names_b2_LR = self.get_col_names(correction_suffix=BPM.__ref_col_lscale_suffix_label,
                                                               include_LR=True, only_b2=True)
        self.__ref_col_orbit_drift_names = self.get_col_names(correction_suffix=BPM.__ref_col_lscale_suffix_label)
        self.__ref_col_orbit_drift_names_LR = self.get_col_names(
            correction_suffix=BPM.__ref_col_orbit_drift_suffix_label, include_LR=True)

        if self.__is_nominal:
            self.__ref_col_diff_names = []
            self.__ref_col_diff_names_LR = []
        else:
            self.__ref_col_diff_names = self.get_col_names(is_diff=True)
            self.__ref_col_diff_names_LR = self.get_col_names(is_diff=True, include_LR=True)
            self.__ref_col_diff_names_b1_LR = self.get_col_names(is_diff=True, include_LR=True, only_b1=True)
            self.__ref_col_diff_names_b2_LR = self.get_col_names(is_diff=True, include_LR=True, only_b2=True)

            if self.__lscale_applied_in_nominal:
                lscale_nominal_suffix = BPM.__ref_col_lscale_suffix_label
            else:
                lscale_nominal_suffix = ""
            self.__ref_col_lscale_diff_names = self.get_col_names(is_diff=True,
                                                                  correction_suffix=BPM.__ref_col_lscale_suffix_label,
                                                                  nominal_correction_suffix=lscale_nominal_suffix)
            self.__ref_col_lscale_diff_names_LR = self.get_col_names(is_diff=True, include_LR=True,
                                                                     correction_suffix=BPM.__ref_col_lscale_suffix_label,
                                                                     nominal_correction_suffix=lscale_nominal_suffix)

            self.__ref_col_lscale_diff_names_b1_LR = self.get_col_names(
                is_diff=True, include_LR=True, only_b1=True, correction_suffix=BPM.__ref_col_lscale_suffix_label,
                nominal_correction_suffix=lscale_nominal_suffix)
            self.__ref_col_lscale_diff_names_b2_LR = self.get_col_names(
                is_diff=True, include_LR=True, only_b2=True, correction_suffix=BPM.__ref_col_lscale_suffix_label,
                nominal_correction_suffix=lscale_nominal_suffix)

            # deflection is applied in Nominal
            self.__ref_col_deflection_diff_names = self.get_col_names(
                is_diff=True, nominal_correction_suffix=BPM.__ref_col_deflection_suffix_label)
            self.__ref_col_deflection_diff_names_LR = self.get_col_names(
                is_diff=True, nominal_correction_suffix=BPM.__ref_col_deflection_suffix_label, include_LR=True)
            self.__ref_col_deflection_diff_names_b1_LR = self.get_col_names(
                is_diff=True, nominal_correction_suffix=BPM.__ref_col_deflection_suffix_label, include_LR=True,
                only_b1=True)
            self.__ref_col_deflection_diff_names_b2_LR = self.get_col_names(
                is_diff=True, nominal_correction_suffix=BPM.__ref_col_deflection_suffix_label, include_LR=True,
                only_b2=True)

            # orbit drift is applied in DOROS, etc
            self.__ref_col_orbit_drift_diff_names = self.get_col_names(
                is_diff=True, correction_suffix=BPM.__ref_col_orbit_drift_suffix_label)
            self.__ref_col_orbit_drift_diff_names_LR = self.get_col_names(
                is_diff=True, correction_suffix=BPM.__ref_col_orbit_drift_suffix_label, include_LR=True)

            # combined corrections
            combined_correction_data_suffix = BPM.__ref_col_lscale_suffix_label
            combined_correction_nominal_suffix = BPM.__ref_col_deflection_suffix_label + \
                                                 BPM.__ref_col_lscale_suffix_label

            self.__ref_col_deflection_and_lscale_diff_names = self.get_col_names(
                is_diff=True, correction_suffix=combined_correction_data_suffix,
                nominal_correction_suffix=combined_correction_nominal_suffix)
            self.__ref_col_deflection_and_lscale_diff_names_LR = self.get_col_names(
                is_diff=True, correction_suffix=combined_correction_data_suffix,
                nominal_correction_suffix=combined_correction_nominal_suffix, include_LR=True)
            self.__ref_col_deflection_and_lscale_diff_names_b1_LR = self.get_col_names(
                is_diff=True, correction_suffix=combined_correction_data_suffix,
                nominal_correction_suffix=combined_correction_nominal_suffix, include_LR=True, only_b1=True)
            self.__ref_col_deflection_and_lscale_diff_names_b2_LR = self.get_col_names(
                is_diff=True, correction_suffix=combined_correction_data_suffix,
                nominal_correction_suffix=combined_correction_nominal_suffix, include_LR=True, only_b2=True)

            self.__ref_col_orbit_drift_and_lscale_diff_names = self.get_col_names(
                is_diff=True, correction_suffix=BPM.__ref_col_orbit_drift_suffix_label +
                                                BPM.__ref_col_lscale_suffix_label,
                nominal_correction_suffix=BPM.__ref_col_lscale_suffix_label)
            self.__ref_col_orbit_drift_and_lscale_diff_names_LR = self.get_col_names(
                is_diff=True, correction_suffix=BPM.__ref_col_orbit_drift_suffix_label +
                                                BPM.__ref_col_lscale_suffix_label,
                nominal_correction_suffix=BPM.__ref_col_lscale_suffix_label, include_LR=True)

            self.__ref_col_deflection_and_orbit_and_lscale_diff_names = self.get_col_names(
                is_diff=True,
                correction_suffix=combined_correction_data_suffix + BPM.__ref_col_orbit_drift_suffix_label,
                nominal_correction_suffix=combined_correction_nominal_suffix)
            self.__ref_col_deflection_and_orbit_and_lscale_diff_names_LR = self.get_col_names(
                is_diff=True,
                correction_suffix=combined_correction_data_suffix + BPM.__ref_col_orbit_drift_suffix_label,
                nominal_correction_suffix=combined_correction_nominal_suffix, include_LR=True)

            self.__ref_col_diff_names_final_result = self.__ref_col_diff_names

    def getDataPathFromFillNumber(self, fill: int):
        path_to_fill_data = setts.data_base_folder + self.name + "/Fill" + str(fill) + "/"
        if self.name in list(setts.config_dict):
            data_path = path_to_fill_data + setts.config_dict[self.name][fill][setts.conf_label_data_file_path]
        else:
            raise AssertionError("detector not configured in settings")

        return data_path

    def get_scans_data(self):
        if self.__is_nominal:
            names = self.__settings[setts.conf_label_scans_names]
            scans_times = self.__settings[setts.conf_label_scans_time_windows]
        else:
            names = setts.config_dict[BPM.__nominal_name][self.__fill][setts.conf_label_scans_names]
            scans_times = setts.config_dict[BPM.__nominal_name][self.__fill][setts.conf_label_scans_time_windows]

        assert len(names) == len(scans_times)
        assert sorted(scans_times) == scans_times

        modified_scans_times = []
        for i in range(0, len(names)):
            mod_low_limit = scans_times[i][0]
            mod_top_limit = scans_times[i][1]

            if i < len(names) - 1:
                delta_top = (scans_times[i+1][0] - scans_times[i][1])/2
            else:
                delta_top = setts.delta_scan_time_window

            if i == 0:
                delta_low = setts.delta_scan_time_window
            else:
                delta_low = (scans_times[i][0] - scans_times[i-1][1]) / 2

            # Check if modified times are inside excluded regions
            for excluded_range in self.__times_to_exclude:
                if excluded_range[0] <= mod_low_limit - delta_low <= excluded_range[1]:
                    delta_low = setts.delta_scan_time_window
                    break
            for excluded_range in self.__times_to_exclude:
                if excluded_range[0] <= mod_top_limit + delta_top <= excluded_range[1]:
                    delta_top = setts.delta_scan_time_window
                    break

            mod_low_limit -= delta_low
            mod_top_limit += delta_top

            times_i_min = np.array([scans_times[i][0] - self.__zero_time, scans_times[i][1] - self.__zero_time]) * \
                          BPM.__scale_time

            mod_times_i_min = np.array([mod_low_limit - self.__zero_time, mod_top_limit - self.__zero_time]) * \
                              BPM.__scale_time

            self.__scans_info_dict[names[i]] = {"time_range": scans_times[i],
                                                "time_range_minutes": times_i_min,
                                                "mod_time_range_minutes": mod_times_i_min,
                                                "mean_time": int((scans_times[i][0] + scans_times[i][1]) / 2),
                                                "mean_time_minutes": (times_i_min[0] + times_i_min[1]) / 2}

            self.__scans_info_for_plotting[names[i]] = (times_i_min[0] + times_i_min[1]) / 2

        print(self.__scans_info_for_plotting)

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
                    if setts.conf_label_unit_position_factor_to_um in self.__settings_list:
                        rescale_factor = self.__settings[setts.conf_label_unit_position_factor_to_um]
                        print(
                            "     -->> column " + i_col + " scaled using values from settings: " + str(rescale_factor))
                    else:
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

    def compute_extra_raw_avg_cols(self, in_data: pd.DataFrame):
        data_in_unified_format = in_data
        if self.name == BPM.__doros_name or self.name == BPM.__arcBPM_name:
            col_b1h = []
            col_b1v = []
            col_b2h = []
            col_b2v = []
            for i_iter in range(0, len(data_in_unified_format)):
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

                col_b1h.append(mean_value_LR_b1h)
                col_b1v.append(mean_value_LR_b1v)
                col_b2h.append(mean_value_LR_b2h)
                col_b2v.append(mean_value_LR_b2v)

            data_in_unified_format[BPM.__col_b1h] = np.array(col_b1h)
            data_in_unified_format[BPM.__col_b1v] = np.array(col_b1v)
            data_in_unified_format[BPM.__col_b2h] = np.array(col_b2h)
            data_in_unified_format[BPM.__col_b2v] = np.array(col_b2v)

        return data_in_unified_format

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

        if self.name == BPM.__doros_name or self.name == BPM.__arcBPM_name:
            raw_col_names = BPM.__ref_col_names_LR + BPM.__ref_col_names
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

    def apply_offsets(self, in_data: pd.DataFrame):
        data_in_unified_format = in_data

        if self.name == BPM.__doros_name or self.name == BPM.__arcBPM_name:
            offsets_computed = False
            offsets_for_LR_cols_available = False
            if setts.conf_label_compute_offsets in self.__settings_list or \
                    setts.conf_label_offset_values not in self.__settings_list:
                if self.__settings[setts.conf_label_compute_offsets] or \
                        setts.conf_label_offset_values not in self.__settings_list:
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

                value_b1h = data_in_unified_format[BPM.__col_b1h][i_iter]
                value_b1v = data_in_unified_format[BPM.__col_b1v][i_iter]
                value_b2h = data_in_unified_format[BPM.__col_b2h][i_iter]
                value_b2v = data_in_unified_format[BPM.__col_b2v][i_iter]

                col_b1h.append(value_b1h - off_b1h)
                col_b1v.append(value_b1v - off_b1v)
                col_b2h.append(value_b2h - off_b2h)
                col_b2v.append(value_b2v - off_b2v)

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

    # It gets the difference between in_data specified cols and Nominal with the respective correction labels
    def get_cols_diff(self, in_data, base_col_names: list, correction_suffix: str = '',
                      nominal_correction_suffix: str = '', input_nominal_data=None):
        if input_nominal_data is None:
            input_nominal_data = self.__nominal_data
        for i_col in base_col_names:
            col_name_nominal = get_base_name_from_DOROS_LR(i_col) + nominal_correction_suffix
            col_name_data = i_col + correction_suffix
            in_data[col_name_data + BPM.__ref_col_diff_suffix_label + col_name_nominal] = \
                in_data[col_name_data] - input_nominal_data[col_name_nominal]

    def get_col_names(self, correction_suffix: str = '', is_diff: bool = False,
                      nominal_correction_suffix: str = '', include_LR: bool = False,
                      only_b1: bool = False, only_b2: bool = False):

        result_col_names = []
        cols_names_nominal = {}

        if include_LR:
            if only_b1:
                base_cols = BPM.__cols_b1_names_LR
            elif only_b2:
                base_cols = BPM.__cols_b2_names_LR
            else:
                base_cols = BPM.__ref_col_names_LR
        else:
            if only_b1:
                base_cols = BPM.__cols_b1_names
            elif only_b2:
                base_cols = BPM.__cols_b2_names
            else:
                base_cols = BPM.__ref_col_names

        for i_col in base_cols:
            cols_names_nominal[i_col] = get_base_name_from_DOROS_LR(i_col)

        for i_col in base_cols:
            if is_diff:
                result_name = i_col + correction_suffix + BPM.__ref_col_diff_suffix_label + \
                              cols_names_nominal[i_col] + nominal_correction_suffix
            else:
                result_name = i_col + correction_suffix
            result_col_names.append(result_name)

        return result_col_names

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

    def get_diffs_when_Nominal_is_zero(self):
        data_to_use = self.__processed_data
        zero_pos_window = setts.zero_position_window

        nominal_in_zero: pd.DataFrame = self.__nominal_data[(self.__nominal_data[BPM.__col_b1h] < zero_pos_window) &
                                                            (self.__nominal_data[BPM.__col_b1h] > -1*zero_pos_window) &
                                                            (self.__nominal_data[BPM.__col_b1v] < zero_pos_window) &
                                                            (self.__nominal_data[BPM.__col_b1v] > -1*zero_pos_window) &
                                                            (self.__nominal_data[BPM.__col_b2h] < zero_pos_window) &
                                                            (self.__nominal_data[BPM.__col_b2h] > -1*zero_pos_window) &
                                                            (self.__nominal_data[BPM.__col_b2v] < zero_pos_window) &
                                                            (self.__nominal_data[BPM.__col_b2v] > -1*zero_pos_window)
                                                            ].copy()

        for time_range in self.__times_to_exclude:
            low_limit = time_range[0]
            top_limit = time_range[1]
            nominal_in_zero = nominal_in_zero.loc[~nominal_in_zero[BPM.__col_time].between(low_limit, top_limit)].copy()

        mask = data_to_use[BPM.__col_time].isin(nominal_in_zero[BPM.__col_time])
        self.__data_in_zero_beam_position = data_to_use[mask].copy()
        self.__data_in_zero_beam_position.reset_index(drop=True, inplace=True)

        # Getting H Raw differences
        self.__data_in_zero_beam_position[BPM.__col_H_diff] = self.__data_in_zero_beam_position[BPM.__col_b2h] - \
                                                              self.__data_in_zero_beam_position[BPM.__col_b1h]
        # Getting V Raw differences
        self.__data_in_zero_beam_position[BPM.__col_V_diff] = self.__data_in_zero_beam_position[BPM.__col_b2v] - \
                                                              self.__data_in_zero_beam_position[BPM.__col_b1v]

    def apply_length_scale_correction(self, base_cols: list, data_to_use=None):
        print("Applying length scale correction ...")
        if data_to_use is None:
            data_to_use = self.__processed_data
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
        unit_factor = 1.
        orbit_drift_correction_data = pd.read_csv(self.__settings[setts.conf_label_orbit_drift],
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

        if not self.__is_nominal:
            sign_b1 = -1
            sign_b2 = 1
        else:
            sign_b1 = 1
            sign_b2 = -1

        x_gfactor = correction_gfactors[0]
        y_gfactor = correction_gfactors[1]

        if data_to_use is None:
            data_to_use = self.__processed_data

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
                        correction *= sign_b1
                    else:
                        correction *= sign_b2

                corrected_cols_dict[i_col].append(uncorrected_val + correction)

        for i_col_name in basic_cols:
            data_to_use[i_col_name + corrected_suffix] = \
                np.array(corrected_cols_dict[i_col_name])

    def apply_corrections(self, basic_cols: list, apply_lscale: bool = False, apply_deflection: bool = False,
                          apply_orbit_drift: bool = False):

        data_to_use = self.__processed_data
        after_deflection_cols = [basename + BPM.__ref_col_deflection_suffix_label for basename in basic_cols]
        after_orbit_drift_cols = [basename + BPM.__ref_col_orbit_drift_suffix_label for basename in basic_cols]
        after_deflection_orbit_drift_cols = \
            [basename + BPM.__ref_col_orbit_drift_suffix_label for basename in after_deflection_cols]

        if apply_deflection:
            self.apply_per_time_range_addition_correction(mode="deflection",
                                                          basic_cols=basic_cols,
                                                          data_to_use=data_to_use)

        if not self.__is_nominal and (apply_deflection or self.__deflection_applied_in_nominal):
            correction_suffix_data_for_deflection_diff = ""
            correction_suffix_nominal_for_deflection_diff = ""
            if apply_deflection:
                correction_suffix_data_for_deflection_diff = BPM.__ref_col_deflection_suffix_label
            elif self.__deflection_applied_in_nominal:
                correction_suffix_nominal_for_deflection_diff = BPM.__ref_col_deflection_suffix_label

            self.get_cols_diff(in_data=data_to_use,
                               base_col_names=basic_cols,
                               correction_suffix=correction_suffix_data_for_deflection_diff,
                               nominal_correction_suffix=correction_suffix_nominal_for_deflection_diff)

        if apply_orbit_drift:
            # Applying only to data
            self.apply_per_time_range_addition_correction(mode="orbit drift",
                                                          basic_cols=basic_cols, data_to_use=data_to_use)
            self.get_cols_diff(in_data=data_to_use,
                               base_col_names=basic_cols,
                               correction_suffix=BPM.__ref_col_orbit_drift_suffix_label,
                               nominal_correction_suffix="")

        if apply_lscale:
            self.apply_length_scale_correction(base_cols=basic_cols, data_to_use=data_to_use)

            # Getting basic lscale difference
            if not self.__is_nominal:
                if self.__lscale_applied_in_nominal:
                    lscale_nominal_suffix_for_diff = BPM.__ref_col_lscale_suffix_label
                else:
                    lscale_nominal_suffix_for_diff = ""
                self.get_cols_diff(in_data=self.__processed_data,
                                   base_col_names=basic_cols,
                                   correction_suffix=BPM.__ref_col_lscale_suffix_label,
                                   nominal_correction_suffix=lscale_nominal_suffix_for_diff)

            if apply_deflection:
                self.apply_length_scale_correction(base_cols=after_deflection_cols, data_to_use=data_to_use)

            if not self.__is_nominal and (self.apply_deflection or self.__deflection_applied_in_nominal):
                self.get_cols_diff(in_data=self.__processed_data,
                                   base_col_names=basic_cols,
                                   correction_suffix=BPM.__ref_col_lscale_suffix_label,
                                   nominal_correction_suffix=BPM.__ref_col_deflection_suffix_label +
                                                             BPM.__ref_col_lscale_suffix_label)

            if apply_orbit_drift:
                self.apply_length_scale_correction(base_cols=after_orbit_drift_cols, data_to_use=data_to_use)
                if not self.__is_nominal:
                    self.get_cols_diff(in_data=self.__processed_data,
                                       base_col_names=basic_cols,
                                       correction_suffix=BPM.__ref_col_orbit_drift_suffix_label +
                                                         BPM.__ref_col_lscale_suffix_label,
                                       nominal_correction_suffix=BPM.__ref_col_lscale_suffix_label)

            if apply_deflection and apply_orbit_drift:
                self.apply_length_scale_correction(base_cols=after_deflection_orbit_drift_cols, data_to_use=data_to_use)

            if not self.__is_nominal:
                nominal_correction_suffix = ""
                data_correction_suffix = ""
                if self.__deflection_applied_in_nominal:
                    nominal_correction_suffix += BPM.__ref_col_deflection_suffix_label
                if self.__lscale_applied_in_nominal:
                    nominal_correction_suffix += BPM.__ref_col_lscale_suffix_label

                # this is the case for deflection only applied to Nominal
                if self.__apply_orbit_drift:
                    data_correction_suffix += BPM.__ref_col_orbit_drift_suffix_label
                data_correction_suffix += BPM.__ref_col_lscale_suffix_label

                # Getting B1diff = B1doros x dorosLS - (B1nominal + deflection) x nominalLS
                self.get_cols_diff(in_data=self.__processed_data,
                                   base_col_names=basic_cols,
                                   correction_suffix=data_correction_suffix,
                                   nominal_correction_suffix=nominal_correction_suffix)

        # Get final result column names
        nominal_suffix_for_final_result = ""
        data_suffix_for_final_result = ""

        if not self.__is_nominal:
            # Nominal
            if self.__deflection_applied_in_nominal:
                nominal_suffix_for_final_result += BPM.__ref_col_deflection_suffix_label
            if self.__lscale_applied_in_nominal:
                nominal_suffix_for_final_result += BPM.__ref_col_lscale_suffix_label
            # Data
            if apply_orbit_drift:
                data_suffix_for_final_result += BPM.__ref_col_orbit_drift_suffix_label
            if apply_lscale:
                data_suffix_for_final_result += BPM.__ref_col_lscale_suffix_label

            self.__ref_col_diff_names_final_result = self.get_col_names(
                correction_suffix=data_suffix_for_final_result,
                nominal_correction_suffix=nominal_suffix_for_final_result, is_diff=True)

    def get_final_V_H_diff_cols(self):
        if len(self.__ref_col_diff_names_final_result) > 0:
            cols_to_use = self.__ref_col_diff_names_final_result
        else:
            cols_to_use = BPM.__ref_col_names
        data_to_use = self.__processed_data

        matching_to_basename_dict = {}

        for base_name in BPM.__ref_col_names:
            for i_name in cols_to_use:
                if base_name in i_name:
                    matching_to_basename_dict[base_name] = i_name
                    break

        assert len(matching_to_basename_dict) == len(BPM.__ref_col_names)

        # Getting H differences
        data_to_use[BPM.__col_H_diff] = data_to_use[matching_to_basename_dict[BPM.__col_b2h]] - \
                                        data_to_use[matching_to_basename_dict[BPM.__col_b1h]]
        # Getting V differences
        data_to_use[BPM.__col_V_diff] = data_to_use[matching_to_basename_dict[BPM.__col_b2v]] - \
                                        data_to_use[matching_to_basename_dict[BPM.__col_b1v]]

    def plot_detector_all_data(self, save_data_to_file=True):
        plot_name = self.name + "_all_data"
        cols_to_plot = self.__cols_after_renaming
        cols_to_plot.remove(BPM.__col_time)

        ltools.color_print("plotting " + plot_name + "...", "green")

        self.__plt_plots[plot_name] = plotting.scatter_plot_from_pandas_frame(self.__in_data_def_format,
                                                                              x_data_label=BPM.__col_time,
                                                                              y_data_label=cols_to_plot,
                                                                              xlabel="time",
                                                                              ylabel="beam position", ncol_legend=4,
                                                                              plot_style='-',
                                                                              label_cms_status=False
                                                                              ).get_figure()

    def plot_detector_data(self, xrange=None, extra_name_suffix="", cols_to_plot=None,
                           data_to_use=None, plot_scan_info=False, plot_scan_limits_lines=False,
                           scans_limits_to_use="time_range_minutes"):
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

        if cols_to_plot is None:
            cols_to_plot = BPM.__ref_col_names

        if data_to_use is None:
            data_to_use = self.__in_data_def_format

        ltools.color_print("plotting " + plot_name + "...", "green")

        ymin = self.__y_range[0]
        ymax = self.__y_range[1]

        if plot_scan_info:
            if self.__scans_info_dict is None:
                print("     -->> Warning: trying to plot scan info but not properly configured! Plot name: " +
                      plot_name)
                draw_labels_pos_dict = None
                draw_lines_scans_limits = None
            else:
                draw_labels_pos_dict = {}
                draw_lines_scans_limits = []
                y_pos_scans = ymin + abs(ymin * setts.delta_pos_for_scan_labels)
                for i_text in list(self.__scans_info_for_plotting):
                    draw_labels_pos_dict[i_text] = [self.__scans_info_for_plotting[i_text], y_pos_scans]
                    draw_lines_scans_limits.extend(self.__scans_info_dict[i_text][scans_limits_to_use])
        else:
            draw_labels_pos_dict = None
            draw_lines_scans_limits = None

        if not plot_scan_limits_lines:
            draw_lines_scans_limits = None

        self.__plt_plots[plot_name] = plotting.scatter_plot_from_pandas_frame(data_to_use,
                                                                              x_data_label=BPM.__col_time_min,
                                                                              y_data_label=list(cols_to_plot),
                                                                              ymin=ymin,
                                                                              ymax=ymax,
                                                                              xmin=x_min, xmax=x_max,
                                                                              xlabel="time [min]",
                                                                              ylabel="beam position", ncol_legend=4,
                                                                              plot_style='-',
                                                                              label_cms_status=False,
                                                                              draw_labels_pos_dict=draw_labels_pos_dict,
                                                                              draw_vertical_line_pos=draw_lines_scans_limits
                                                                              ).get_figure()

    def plot_detector_data_diff_with_nominal(self, xrange=None, yrange=None, extra_name_suffix="",
                                             cols_to_plot=None, data_to_use=None, legend_labels=None,
                                             y_label=None, use_smaller_y_range=False, canvas_square_shape=False,
                                             marker_size=1.0, leg_marker_scale=7, plot_scan_info=False,
                                             plot_scan_limits_lines=False, scans_limits_to_use="time_range_minutes"):
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

        if y_label is None:
            y_label = self.name + " - Nominal diff." + " [" + BPM.__distance_unit + "]"
        if cols_to_plot is None:
            cols_to_plot = self.__ref_col_diff_names
        if data_to_use is None:
            data_to_use = self.__processed_data

        ltools.color_print("plotting " + plot_name, 'green')
        print("     using columns: " + str(cols_to_plot))

        if yrange is None:
            if use_smaller_y_range:
                ymin = self.__y_diff_smaller_range[0]
                ymax = self.__y_diff_smaller_range[1]
            else:
                ymin = self.__y_diff_range[0]
                ymax = self.__y_diff_range[1]
        else:
            ymin = yrange[0]
            ymax = yrange[1]

        canvas_shape = 'nsq'
        if canvas_square_shape:
            canvas_shape = 'sq'

        if plot_scan_info:
            if self.__scans_info_dict is None:
                print("     -->> Warning: trying to plot scan info but not properly configured! Plot name: " +
                      plot_name)
                draw_labels_pos_dict = None
                draw_lines_scans_limits = None
            else:
                draw_labels_pos_dict = {}
                draw_lines_scans_limits = []
                y_pos_scans = ymin + abs(ymin * setts.delta_pos_for_scan_labels)
                for i_text in list(self.__scans_info_for_plotting):
                    draw_labels_pos_dict[i_text] = [self.__scans_info_for_plotting[i_text], y_pos_scans]
                    draw_lines_scans_limits.extend(self.__scans_info_dict[i_text][scans_limits_to_use])
        else:
            draw_labels_pos_dict = None
            draw_lines_scans_limits = None

        if not plot_scan_limits_lines:
            draw_lines_scans_limits = None

        self.__plt_plots[plot_name] = plotting.scatter_plot_from_pandas_frame(data_to_use,
                                                                              x_data_label=BPM.__col_time_min,
                                                                              y_data_label=list(cols_to_plot),
                                                                              ymin=ymin,
                                                                              ymax=ymax,
                                                                              xmin=x_min, xmax=x_max,
                                                                              xlabel="time [min]",
                                                                              ylabel=y_label,
                                                                              ncol_legend=4,
                                                                              plot_style='o',
                                                                              leg_marker_sc=leg_marker_scale,
                                                                              label_cms_status=False,
                                                                              legend_labels=legend_labels,
                                                                              fig_size_shape=canvas_shape,
                                                                              marker_size=marker_size,
                                                                              draw_labels_pos_dict=draw_labels_pos_dict,
                                                                              draw_vertical_line_pos=draw_lines_scans_limits
                                                                              ).get_figure()

    def save_plots(self):
        ltools.color_print('\n\n Saving plots:', "green")
        plotting.save_plots(self.__plt_plots, self.__output_dir, save_pickle=setts.save_figures_as_pickle)
        plotting.save_plots(self.__sns_plots, self.__output_dir, save_pickle=setts.save_figures_as_pickle)

    @property
    def apply_lscale(self):
        return self.__apply_lscale

    @property
    def apply_orbit_drift(self):
        return self.__apply_orbit_drift

    @property
    def apply_deflection(self):
        return self.__apply_deflection

    @property
    def zero_time(self):
        return self.__zero_time


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
