import pandas as pd
import settings_bpm as setts
import settings as def_setts
import tools.plotting_tools as plotting
import tools.live_plotting_tools as live_plotting
from scipy import interpolate
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
# import json
import numpy as np
import tools.lumi_tools as ltools
from lmfit import Model

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
    __col_H_sum = "H_sum"
    __col_V_sum = "V_sum"

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
    __cols_H_names = (__col_b1h, __col_b2h)
    __cols_V_names = (__col_b1v, __col_b2v)
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
                 get_only_data: bool = False, apply_all_available_corrections=True, compute_orbit_drift=False,
                 nominal_data=None, get_orbit_drift_plots=False) -> None:
        self.name = name
        self.__debugging = setts.debugging
        if self.name not in BPM.__allowed_detectors:
            raise AssertionError("Detector " + self.name + " no implemented")
        if self.name not in list(BPM.__col_names_to_read):
            raise AssertionError("Detector " + self.name + " data reading setting not available")

        self.__output_dir = "plots_bpm/" + str(fill) + "/" + self.name + "/"
        self.__output_per_scan_studies = "per_scan_studies/"
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
        self.__scans_timestamp_limits_dict = {}
        self.__scans_info_for_plotting = {}
        self.__scans_with_inverted_beam_sign = []
        self.__orbit_drifts_per_scan = {}
        self.__data_per_scan_dict = {}
        self.__data_per_scan_nominal_dict = {}
        self.__lscale_per_beam_epsilon_dict = {}
        self.__lscale_per_beam_lscale_slope_dict = {}
        # computed mean values
        self.__lscale_per_beam_slope_b1_x = None
        self.__lscale_per_beam_slope_b1_y = None
        self.__lscale_per_beam_slope_b2_x = None
        self.__lscale_per_beam_slope_b2_y = None
        self.__lscale_per_beam_epsilon_b1_x = None
        self.__lscale_per_beam_epsilon_b1_y = None
        self.__lscale_per_beam_epsilon_b2_x = None
        self.__lscale_per_beam_epsilon_b2_y = None
        self.__lscale_per_beam_slope_b1_x_err = None
        self.__lscale_per_beam_slope_b1_y_err = None
        self.__lscale_per_beam_slope_b2_x_err = None
        self.__lscale_per_beam_slope_b2_y_err = None
        self.__lscale_per_beam_epsilon_b1_x_err = None
        self.__lscale_per_beam_epsilon_b1_y_err = None
        self.__lscale_per_beam_epsilon_b2_x_err = None
        self.__lscale_per_beam_epsilon_b2_y_err = None
        self.__epsilon_x_in_plotting_structure_dict = {}
        self.__epsilon_y_in_plotting_structure_dict = {}
        self.__slope_x_in_plotting_structure_dict = {}
        self.__slope_y_in_plotting_structure_dict = {}
        self.__loaded_epsilon_correction_dict = {}

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
        self.__nominal_suffix_for_final_result = ""
        self.__data_suffix_for_final_result = ""

        self.fill_col_names()

        self.__settings = setts.config_dict[self.name][fill]
        self.__settings_list = list(self.__settings)

        if self.__is_nominal:
            self.__nominal_settings = self.__settings
            self.__nominal_settings_list = self.__settings_list
        else:
            self.__nominal_settings = self.__nominal_det.settings
            self.__nominal_settings_list = self.__nominal_det.settings_list

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

        if setts.conf_label_orbit_drift_plot_range in self.__settings_list:
            self.__y_range_orbit_drift = self.__settings[setts.conf_label_orbit_drift_plot_range]
        else:
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

        ltools.check_and_create_folder(self.__output_dir)

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

        self.__in_data_def_format = self.convert_cols_to_default_format(in_data)

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
            self.__ref_col_diff_names_final_result_b1 = self.get_col_names(is_diff=True, only_b1=True)
            self.__ref_col_diff_names_final_result_b2 = self.get_col_names(is_diff=True, only_b2=True)
            self.__ref_col_diff_names_final_result_H = self.get_col_names(is_diff=True, only_H=True)
            self.__ref_col_diff_names_final_result_V = self.get_col_names(is_diff=True, only_V=True)

            self.add_time_min_col(self.__processed_data)

            if compute_orbit_drift:
                orbit_drift_computed = False
                self.get_diffs_when_nominal_is_zero()
                try:
                    self.compute_orbit_drifts_per_limit_scan_points()
                    orbit_drift_computed = True
                except:
                    ltools.color_print("\n ** WARNING ** : Orbit drifts not computed \n", 'yellow')
                    orbit_drift_computed = False

                if not get_only_data or get_orbit_drift_plots:
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
                    if orbit_drift_computed:
                        self.plot_orbit_drift_result()
        else:
            self.__interpolation_to_nominal_time = self.__in_data_def_format
            self.__processed_data = self.__in_data_def_format

        # #  ### getting initial diffs
        # # Getting H differences
        # self.__processed_data[BPM.__col_H_diff + "_ini"] = \
        #     self.__processed_data[self.get_single_diff_name(BPM.__col_b2h)] - \
        #     self.__processed_data[self.get_single_diff_name(BPM.__col_b1h)]
        # # Getting V differences
        # self.__processed_data[BPM.__col_V_diff + "_ini"] = \
        #     self.__processed_data[self.get_single_diff_name(BPM.__col_b2v)] - \
        #     self.__processed_data[self.get_single_diff_name(BPM.__col_b1v)]

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
            if def_setts.show_live_plots:
                self.plot_live_detector_data()

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
                                                          legend_labels=self.__ref_col_names,
                                                          plot_scan_info=self.__scan_info_available,
                                                          plot_scan_limits_lines=self.__scan_info_available
                                                          )
                self.plot_detector_data_diff_with_nominal(extra_name_suffix="_final_b1",
                                                          cols_to_plot=self.__ref_col_diff_names_final_result_b1,
                                                          data_to_use=self.__processed_data,
                                                          legend_labels=self.__cols_b1_names)
                self.plot_detector_data_diff_with_nominal(extra_name_suffix="_final_b2",
                                                          cols_to_plot=self.__ref_col_diff_names_final_result_b2,
                                                          data_to_use=self.__processed_data,
                                                          legend_labels=self.__cols_b2_names)
                self.plot_detector_data_diff_with_nominal(extra_name_suffix="_per_beam_final_H",
                                                          cols_to_plot=self.__ref_col_diff_names_final_result_H,
                                                          data_to_use=self.__processed_data,
                                                          legend_labels=self.__cols_H_names)
                self.plot_detector_data_diff_with_nominal(extra_name_suffix="_per_beam_final_V",
                                                          cols_to_plot=self.__ref_col_diff_names_final_result_V,
                                                          data_to_use=self.__processed_data,
                                                          legend_labels=self.__cols_V_names)
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
                                                              plot_scan_info=self.__scan_info_available,
                                                              plot_scan_limits_lines=self.__scan_info_available,
                                                              y_label="B2 - B1" +
                                                                      " [" + BPM.__distance_unit + "]")
                    self.plot_detector_data_diff_with_nominal(extra_name_suffix="_H_V_sum_final",
                                                              cols_to_plot=[BPM.__col_H_sum, BPM.__col_V_sum],
                                                              data_to_use=self.__processed_data,
                                                              legend_labels=["H", "V"],
                                                              y_label="B2 + B1" +
                                                                      " [" + BPM.__distance_unit + "]")

                # per scan studies
                self.do_per_scan_studies()

            self.save_plots()

            if not self.__is_nominal:
                #   Saving data to file
                data_to_save = self.__processed_data
                timing_cols = [BPM.__col_time, BPM.__col_time_min]
                data_cols = self.__ref_col_diff_names_final_result
                all_cols_to_save = timing_cols + data_cols
                raw_data_cols = timing_cols + list(BPM.__ref_col_names) + self.__ref_col_diff_names
                final_hysteresis_cols = timing_cols + [BPM.__col_H_diff, BPM.__col_V_diff]

                # make sure to store data only in studied range
                ltools.color_print("\n\n (OUTPUT FILE) Saving data in " + self.__output_dir + "plotted_data.csv",
                                   "yellow")
                ltools.color_print(" (OUTPUT FILE) Saving data with only offsets in " + self.__output_dir +
                                   "_only_offsets_data.csv", "yellow")
                ltools.color_print(" (OUTPUT FILE) Saving hysteresis data in " + self.__output_dir + str(self.__fill) +
                                   "_hysteresis.csv", "yellow")
                # ltools.color_print("    Data saved only in range " + str(self.__y_diff_range), "blue")
                # for i_col in data_cols:
                #     data_to_save = data_to_save[(data_to_save[i_col] >= self.__y_diff_range[0]) &
                #                                 (data_to_save[i_col] <= self.__y_diff_range[1])]
                ltools.save_columns_from_pandas_to_file(data_to_save, all_cols_to_save, self.__output_dir +
                                                        self.name + "_" + str(self.__fill) + "_plotted_data.csv")
                ltools.save_columns_from_pandas_to_file(data_to_save, raw_data_cols, self.__output_dir +
                                                        self.name + "_" + str(self.__fill) + "_only_offsets_data.csv")
                ltools.save_columns_from_pandas_to_file(data_to_save, final_hysteresis_cols, self.__output_dir +
                                                        self.name + "_" + str(self.__fill) + "_hysteresis.csv")

                # ltools.save_columns_from_pandas_to_file(data_to_save, self.__output_dir +
                #                                         "plotted_data_all.csv")

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
                delta_top = (scans_times[i + 1][0] - scans_times[i][1]) / 2
            else:
                delta_top = setts.delta_scan_time_window

            if i == 0:
                delta_low = setts.delta_scan_time_window
            else:
                delta_low = (scans_times[i][0] - scans_times[i - 1][1]) / 2

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
            mod_low_limit = int(mod_low_limit)
            mod_top_limit = int(mod_top_limit)

            times_i_min = np.array([scans_times[i][0] - self.__zero_time, scans_times[i][1] - self.__zero_time]) * \
                          BPM.__scale_time

            mod_times_i_min = np.array([mod_low_limit - self.__zero_time, mod_top_limit - self.__zero_time]) * \
                              BPM.__scale_time

            self.__scans_info_dict[names[i]] = {"time_range": scans_times[i],
                                                "time_range_minutes": times_i_min,
                                                "mod_time_range_minutes": mod_times_i_min,
                                                "mod_time_range": [mod_low_limit, mod_top_limit],
                                                "mean_time": int((scans_times[i][0] + scans_times[i][1]) / 2),
                                                "mean_time_minutes": (times_i_min[0] + times_i_min[1]) / 2}
            self.__scans_timestamp_limits_dict[names[i]] = scans_times[i]

            self.__scans_info_for_plotting[names[i]] = (times_i_min[0] + times_i_min[1]) / 2

        # print(self.__scans_info_for_plotting)
        ltools.color_print("\n ====>>> Scan info loaded from settings \n", "blue")

        # Reading info about special beam sign for specific scans
        if setts.conf_label_scans_with_inverted_beam_sign in self.__nominal_settings:
            self.__scans_with_inverted_beam_sign = \
                self.__nominal_settings[setts.conf_label_scans_with_inverted_beam_sign]
            ltools.color_print(" ====>>> Special sign is set for scans " + str(self.__scans_with_inverted_beam_sign)
                               + "\n", "blue")

    def convert_cols_to_default_format(self, in_data: pd.DataFrame, rescale: bool = True):
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

    def compute_offsets(self, mode="only after optimization"):
        print("Computing offsets ...")
        if setts.conf_label_offset_time in self.__settings_list:
            self.__offsets_time_split = self.__settings[setts.conf_label_offset_time]
        else:
            raise AssertionError("optimization scans times must be defined in settings")
        available_modes = ("only after optimization", "all head-on")
        if mode not in available_modes:
            raise AssertionError("mode not correct")
        mode_only_short_after_opt = False
        mode_all_head_on = False
        if mode == "only after optimization":
            mode_only_short_after_opt = True
            if setts.conf_label_offset_time not in self.__settings_list:
                raise AssertionError(
                    "The optimization times must be defined for using the _only after optimization_ method")
        elif mode == "all head-on":
            mode_all_head_on = True

        if setts.conf_label_time_window_for_offsets in self.__settings_list:
            time_window_for_offsets = self.__settings[setts.conf_label_time_window_for_offsets]
        else:
            time_window_for_offsets = setts.time_window_for_offsets

        diffs = {}
        offsets = {}
        data_to_use = self.__in_data_def_format
        when_b1_h_is_zero = self.__nominal_data[BPM.__col_b1h] == 0.
        when_b1_v_is_zero = self.__nominal_data[BPM.__col_b1v] == 0.
        when_b2_h_is_zero = self.__nominal_data[BPM.__col_b2h] == 0.
        when_b2_v_is_zero = self.__nominal_data[BPM.__col_b2v] == 0.

        zero_pos_window = setts.zero_position_window_for_offsets

        nominal_in_zero: pd.DataFrame = self.__nominal_data[(self.__nominal_data[BPM.__col_b1h] < zero_pos_window) &
                                                            (self.__nominal_data[
                                                                 BPM.__col_b1h] > -1 * zero_pos_window) &
                                                            (self.__nominal_data[BPM.__col_b1v] < zero_pos_window) &
                                                            (self.__nominal_data[
                                                                 BPM.__col_b1v] > -1 * zero_pos_window) &
                                                            (self.__nominal_data[BPM.__col_b2h] < zero_pos_window) &
                                                            (self.__nominal_data[
                                                                 BPM.__col_b2h] > -1 * zero_pos_window) &
                                                            (self.__nominal_data[BPM.__col_b2v] < zero_pos_window) &
                                                            (self.__nominal_data[BPM.__col_b2v] > -1 * zero_pos_window)
                                                            ].copy()

        mask = data_to_use[BPM.__col_time].isin(nominal_in_zero[BPM.__col_time])
        data_for_offsets = data_to_use[mask].copy()

        offsets_limits = []
        for i_time in self.__offsets_time_split:
            offsets_limits.append([i_time, i_time + time_window_for_offsets])

        plot_data_for_offsets = plotting.scatter_plot_from_pandas_frame(data_for_offsets,
                                                                        x_data_label=BPM.__col_time,
                                                                        y_data_label=list(BPM.__ref_col_names),
                                                                        xlabel="time [timestamp]",
                                                                        ylabel="unscaled positions",
                                                                        ncol_legend=4,
                                                                        label_cms_status=False,
                                                                        draw_vertical_line_pos=self.__settings[
                                                                            setts.conf_label_offset_time],
                                                                        draw_vertical_bands_array_pos=offsets_limits,
                                                                        title="Fill" + str(self.__fill),
                                                                        title_loc="right"
                                                                        ).get_figure()
        plot_nominal_data_for_offsets = plotting.scatter_plot_from_pandas_frame(nominal_in_zero,
                                                                                x_data_label=BPM.__col_time,
                                                                                y_data_label=list(BPM.__ref_col_names),
                                                                                xlabel="time [timestamp]",
                                                                                ylabel="unscaled positions",
                                                                                ncol_legend=4,
                                                                                label_cms_status=False,
                                                                                draw_vertical_line_pos=self.__settings[
                                                                                    setts.conf_label_offset_time],
                                                                                draw_vertical_bands_array_pos=offsets_limits,
                                                                                title="Fill" + str(self.__fill),
                                                                                title_loc="right"
                                                                                ).get_figure()

        plot_data_for_offsets.savefig(self.__output_dir + "offsets_debug_data_for_offsets.png")
        plot_nominal_data_for_offsets.savefig(self.__output_dir + "offsets_debug_nominal_data_for_offsets.png")

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

        if mode_all_head_on:
            if setts.conf_label_offset_time in self.__settings_list:
                self.__offsets_time_split = self.__settings[setts.conf_label_offset_time]
                for index in range(0, len(self.__offsets_time_split)):
                    this_opt_scan_time = self.__offsets_time_split[index]
                    if index == len(self.__offsets_time_split) - 1:
                        diffs_df_split_time_list.append(
                            diffs_df[diffs_df[BPM.__col_time] >= this_opt_scan_time])
                    else:
                        diffs_df_split_time_list.append(
                            diffs_df[(diffs_df[BPM.__col_time] >= this_opt_scan_time)
                                     & (diffs_df[BPM.__col_time] < self.__offsets_time_split[index + 1])])
                time_id = 0
                for i_df in diffs_df_split_time_list:
                    self.__plt_plots["offsets_histos_" + str(self.__offsets_time_split[time_id])] = \
                    i_df.hist(bins=50)[0][
                        0].get_figure()
                    for col_name in raw_col_names:
                        self.__offsets[col_name].append(i_df[col_name].mean())
                    time_id += 1
            else:
                diffs_df_split_time_list = [diffs_df]
                for col_name in raw_col_names:
                    self.__offsets[col_name] = diffs_df[col_name].mean()

        elif mode_only_short_after_opt:
            self.__offsets_time_split = self.__settings[setts.conf_label_offset_time]
            for index in range(0, len(self.__offsets_time_split)):
                this_opt_scan_time = self.__offsets_time_split[index]
                end_of_time_window = this_opt_scan_time + time_window_for_offsets

                diffs_df_split_time_list.append(
                    diffs_df[(diffs_df[BPM.__col_time] >= this_opt_scan_time)
                             & (diffs_df[BPM.__col_time] < end_of_time_window)])
            time_id = 0
            for i_df in diffs_df_split_time_list:
                self.__plt_plots["offsets_histos_" + str(self.__offsets_time_split[time_id])] = i_df.hist(bins=50)[0][
                    0].get_figure()
                for col_name in raw_col_names:
                    self.__offsets[col_name].append(i_df[col_name].mean())
                time_id += 1

        # save computed offsets into file
        offsets_in_plain_list_format_short = []
        for time_index in range(0, len(self.__offsets[BPM.__ref_col_names[0]])):
            col_offsets_in_time = []
            for col_name in BPM.__ref_col_names:
                col_offsets_in_time.append(self.__offsets[col_name][time_index])
            offsets_in_plain_list_format_short.append(col_offsets_in_time)

        with open(self.__output_dir + self.name + "_offsets_all.txt", 'w') as file:
            file.write('Timestamps: ' + str(self.__offsets_time_split) + "\n")
            file.write('OffValues dict: ' + str(self.__offsets) + "\n")
        with open(self.__output_dir + self.name + "_offsets.txt", 'w') as file:
            file.write('Timestamps: ' + str(self.__offsets_time_split) + "\n")
            file.write('OffValues: ' + str(offsets_in_plain_list_format_short) + "\n")
            file.write('Col Names: ' + str(list(BPM.__ref_col_names)) + "\n")

        ltools.color_print("\n (OUTPUT FILE) Computed offsets saved in " +
                           self.__output_dir + self.name + "_offsets.txt \n", "yellow")

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
                    if setts.conf_label_compute_offsets_method in self.__settings_list:
                        self.compute_offsets(mode=self.__settings[setts.conf_label_compute_offsets_method])
                    else:
                        self.compute_offsets()
                    offsets_computed = True
            if setts.conf_label_offset_values in self.__settings_list:
                self.get_offsets_from_config()
                ltools.color_print("\n ====>>> Using offsets defined in settings \n", "blue")

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

    def get_single_diff_name(self, base_col_name: str, correction_suffix: str = '',
                             nominal_correction_suffix: str = ''):
        result_name = base_col_name + correction_suffix + BPM.__ref_col_diff_suffix_label + \
                      base_col_name + nominal_correction_suffix
        return result_name

    def get_col_names(self, correction_suffix: str = '', is_diff: bool = False,
                      nominal_correction_suffix: str = '', include_LR: bool = False,
                      only_b1: bool = False, only_b2: bool = False, only_H: bool = False, only_V: bool = False):

        result_col_names = []
        cols_names_nominal = {}

        if include_LR:
            if only_b1:
                base_cols = BPM.__cols_b1_names_LR
            elif only_b2:
                base_cols = BPM.__cols_b2_names_LR
            elif only_H:
                base_cols = BPM.__cols_H_names
            elif only_V:
                base_cols = BPM.__cols_V_names
            else:
                base_cols = BPM.__ref_col_names_LR
        else:
            if only_b1:
                base_cols = BPM.__cols_b1_names
            elif only_b2:
                base_cols = BPM.__cols_b2_names
            elif only_H:
                base_cols = BPM.__cols_H_names
            elif only_V:
                base_cols = BPM.__cols_V_names
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

    def get_final_col_name_data(self, col_name: str):
        return col_name + self.__data_suffix_for_final_result

    def get_final_col_name_nominal(self, col_name: str):
        return col_name + self.__nominal_suffix_for_final_result

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

    def get_diffs_when_nominal_is_zero(self):
        data_to_use = self.__processed_data
        zero_pos_window = setts.zero_position_window

        nominal_in_zero: pd.DataFrame = self.__nominal_data[(self.__nominal_data[BPM.__col_b1h] < zero_pos_window) &
                                                            (self.__nominal_data[
                                                                 BPM.__col_b1h] > -1 * zero_pos_window) &
                                                            (self.__nominal_data[BPM.__col_b1v] < zero_pos_window) &
                                                            (self.__nominal_data[
                                                                 BPM.__col_b1v] > -1 * zero_pos_window) &
                                                            (self.__nominal_data[BPM.__col_b2h] < zero_pos_window) &
                                                            (self.__nominal_data[
                                                                 BPM.__col_b2h] > -1 * zero_pos_window) &
                                                            (self.__nominal_data[BPM.__col_b2v] < zero_pos_window) &
                                                            (self.__nominal_data[BPM.__col_b2v] > -1 * zero_pos_window)
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

    def compute_orbit_drifts_per_limit_scan_points(self):
        assert len(self.__data_in_zero_beam_position) != 0

        for scan_name in list(self.__scans_info_dict):
            limit_times = self.__scans_info_dict[scan_name]["time_range"]
            limit_times_min = self.__scans_info_dict[scan_name]["time_range_minutes"]
            mod_limit_times = self.__scans_info_dict[scan_name]["mod_time_range"]
            scan_middle_time = self.__scans_info_dict[scan_name]["mean_time"]
            scan_middle_time_min = self.__scans_info_dict[scan_name]["mean_time_minutes"]
            time_window = setts.delta_scan_time_window
            ini_scan_time = limit_times[0]
            end_scan_time = limit_times[1]
            ini_scan_mod_time = mod_limit_times[0]
            end_scan_mod_time = mod_limit_times[1]

            orbit_drift_fit_out_dir = self.__output_dir + "/orbit_drift_special_plots/"

            data_inside_mod_time = \
                self.__data_in_zero_beam_position[
                    (self.__data_in_zero_beam_position[BPM.__col_time] >= mod_limit_times[0]) &
                    (self.__data_in_zero_beam_position[BPM.__col_time] <= mod_limit_times[1])].copy()

            ini_scan_time_window = [ini_scan_mod_time, ini_scan_time + time_window]
            end_scan_time_window = [end_scan_time - time_window, end_scan_mod_time]
            middle_scan_time_window = [scan_middle_time - time_window, scan_middle_time + time_window]
            ini_scan_points = data_inside_mod_time[(data_inside_mod_time[BPM.__col_time] >= ini_scan_time_window[0]) &
                                                   (data_inside_mod_time[BPM.__col_time] <= ini_scan_time_window[
                                                       1])].copy()
            end_scan_points = data_inside_mod_time[(data_inside_mod_time[BPM.__col_time] >= end_scan_time_window[0]) &
                                                   (data_inside_mod_time[BPM.__col_time] <= end_scan_time_window[
                                                       1])].copy()
            middle_scan_points = data_inside_mod_time[
                (data_inside_mod_time[BPM.__col_time] >= middle_scan_time_window[0]) &
                (data_inside_mod_time[BPM.__col_time] <= middle_scan_time_window[1])].copy()

            region_bands = [ini_scan_time_window, middle_scan_time_window, end_scan_time_window]
            all_data_inside_scan = plotting.scatter_plot_from_pandas_frame(data_inside_mod_time,
                                                                           x_data_label=BPM.__col_time,
                                                                           y_data_label=[BPM.__col_H_diff,
                                                                                         BPM.__col_V_diff],
                                                                           xlabel="time [min]",
                                                                           ylabel="positions",
                                                                           ncol_legend=4,
                                                                           label_cms_status=False,
                                                                           draw_vertical_line_pos=[ini_scan_mod_time,
                                                                                                   scan_middle_time,
                                                                                                   end_scan_mod_time],
                                                                           draw_vertical_bands_array_pos=region_bands,
                                                                           draw_vertical_bands_array_alpha=0.2,
                                                                           marker_size=2.0,
                                                                           title=scan_name + " scan",
                                                                           title_loc="right"
                                                                           ).get_figure()

            ltools.check_and_create_folder(orbit_drift_fit_out_dir, creation_info=False)
            all_data_inside_scan.savefig(orbit_drift_fit_out_dir + scan_name + "_scan_all_data.png")

            od_H_in_low_limit = ini_scan_points[BPM.__col_H_diff].mean()
            od_H_in_middle = middle_scan_points[BPM.__col_H_diff].mean()
            od_H_in_top_limit = end_scan_points[BPM.__col_H_diff].mean()

            od_V_in_low_limit = ini_scan_points[BPM.__col_V_diff].mean()
            od_V_in_middle = middle_scan_points[BPM.__col_V_diff].mean()
            od_V_in_top_limit = end_scan_points[BPM.__col_V_diff].mean()

            n_middle_scan_points_H = len(middle_scan_points[BPM.__col_H_diff])
            n_middle_scan_points_V = len(middle_scan_points[BPM.__col_V_diff])

            if "L" in scan_name and (n_middle_scan_points_H == 0 or n_middle_scan_points_V == 0):
                if n_middle_scan_points_H == 0:
                    od_H_in_middle = (od_H_in_low_limit + od_H_in_top_limit) / 2
                if n_middle_scan_points_V == 0:
                    od_V_in_middle = (od_V_in_low_limit + od_V_in_top_limit) / 2

            self.__orbit_drifts_per_scan[scan_name] = {
                "TimeWindows": [limit_times[0], scan_middle_time, limit_times[1]],
                "TimeWindows_min": [limit_times_min[0], scan_middle_time_min, limit_times_min[1]],
                "OrbitDrifts_X": [od_H_in_low_limit, od_H_in_middle, od_H_in_top_limit],
                "OrbitDrifts_Y": [od_V_in_low_limit, od_V_in_middle, od_V_in_top_limit]
            }

        # print(self.__orbit_drifts_per_scan)
        self.store_orbit_drift_in_file()

    def store_orbit_drift_in_file(self):
        scans = list(self.__orbit_drifts_per_scan)
        times = []
        od_x = []
        od_y = []

        for scan in scans:
            times.append(self.__orbit_drifts_per_scan[scan]["TimeWindows"])
            od_x.append(self.__orbit_drifts_per_scan[scan]["OrbitDrifts_X"])
            od_y.append(self.__orbit_drifts_per_scan[scan]["OrbitDrifts_Y"])

        with open(self.__output_dir + "OrbitDrifts_" + self.name + ".json", 'w') as file:
            file.write("{\n")
            file.write('"Names": ' + str(scans) + ",\n")
            file.write('"TimeWindows": ' + str(times) + ",\n")
            file.write('"OrbitDrifts_X": ' + str(od_x) + ",\n")
            file.write('"OrbitDrifts_Y": ' + str(od_y) + "\n")
            file.write("}")

        ltools.color_print("\n (OUTPUT FILE) Orbit drifts saved in " +
                           self.__output_dir + "OrbitDrifts_" + self.name + ".json \n", "yellow")

    def apply_length_scale_correction(self, base_cols: list, data_to_use=None):
        print("Applying length scale correction to cols: " + str(base_cols))
        if data_to_use is None:
            data_to_use = self.__processed_data
        # print(list(data_to_use))

        per_scan_correction_is_needed = False

        x_correction = setts.config_dict[self.name][self.__fill][setts.conf_label_length_scale][0]
        y_correction = setts.config_dict[self.name][self.__fill][setts.conf_label_length_scale][1]

        # Loading overall correction lenght scale per beam residual correction
        m_b1_x = 0.0
        m_b2_x = 0.0
        m_b1_y = 0.0
        m_b2_y = 0.0
        if setts.conf_label_length_scale_per_beam_epsilon in self.__settings_list:
            per_beam_correction = self.__settings[setts.conf_label_length_scale_per_beam_epsilon]
            m_b1_x = per_beam_correction[0][0]
            m_b2_x = per_beam_correction[0][1]
            m_b1_y = per_beam_correction[1][0]
            m_b2_y = per_beam_correction[1][1]

        # Loading overall correction lenght scale per beam residual correction
        elif setts.conf_label_apply_length_scale_per_beam_epsilon_from_file in self.__settings_list:
            if self.__settings[setts.conf_label_apply_length_scale_per_beam_epsilon_from_file]:
                if ltools.check_file_existence(self.__output_dir + self.__output_per_scan_studies +
                                               "lscale_per_beam_epsilon.json"):
                    self.read_length_scale_per_beam_epsilon_file()
                    per_scan_correction_is_needed = True
                else:
                    ltools.color_print("--->>> lscale per beam epsilon is enable but not file found. "
                                       "The file will be produce and then it would be needed to rerun and check "
                                       "that this message it is not printed again", "red")

        h_correction = x_correction
        v_correction = y_correction

        corrected_cols_dict = {}

        for i_col_name in base_cols:
            corrected_cols_dict[i_col_name] = []

        if not per_scan_correction_is_needed:
            for index in range(0, len(data_to_use)):
                for i_col_name in base_cols:
                    if is_H(i_col_name):
                        correction = h_correction
                        if is_B1(i_col_name):
                            epsilon = m_b1_x
                        else:
                            epsilon = m_b2_x
                    elif is_V(i_col_name):
                        correction = v_correction
                        if is_B1(i_col_name):
                            epsilon = m_b1_y
                        else:
                            epsilon = m_b2_y
                    else:
                        raise AssertionError()
                    correction -= epsilon
                    corrected_cols_dict[i_col_name].append(data_to_use[i_col_name][index] * correction)
        else:
            for index in range(0, len(data_to_use)):
                for i_col_name in base_cols:
                    data_time = data_to_use[BPM.__col_time][index]
                    scan_name = self.get_scan_name_from_timestamp(data_time)
                    correction = self.get_lscale_correction_value_per_scan(i_col_name, scan_name,
                                                                           h_correction, v_correction)
                    corrected_cols_dict[i_col_name].append(data_to_use[i_col_name][index] * correction)

        for i_col_name in base_cols:
            data_to_use[i_col_name + BPM.__ref_col_lscale_suffix_label] = \
                np.array(corrected_cols_dict[i_col_name])

    def read_length_scale_per_beam_epsilon_file(self):
        file_location = self.__output_dir + self.__output_per_scan_studies + "lscale_per_beam_epsilon.json"
        self.__loaded_epsilon_correction_dict = ltools.load_json_as_dict(file_location)
        # print(self.__loaded_epsilon_correction_dict)
        ltools.color_print("--->>> lscale_per_beam_epsilon.json file for per beam lscale correction has been loaded",
                           color="blue")

    def get_lscale_correction_value_per_scan(self, i_col_name, scan_name, h_correction, v_correction):
        base_col_name = get_base_name(i_col_name)
        if scan_name != "":
            if (is_H(i_col_name) and "X" in scan_name) or (is_V(i_col_name) and "Y" in scan_name):
                epsilon = self.__loaded_epsilon_correction_dict[scan_name][base_col_name]
            # elif is_H(i_col_name) and "Y" in scan_name:
            #     epsilon = self.__loaded_epsilon_correction_dict[scan_name.replace("Y", "X")][base_col_name]
            # elif is_V(i_col_name) and "X" in scan_name:
            #     epsilon = self.__loaded_epsilon_correction_dict[scan_name.replace("X", "Y")][base_col_name]
            else:
                epsilon = 0.0
        else:
            epsilon = 0.0

        if is_H(i_col_name):
            correction = h_correction
        elif is_V(i_col_name):
            correction = v_correction
        else:
            raise AssertionError()

        correction -= epsilon

        return correction

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

    def get_beam_sign_for_column_name(self, col_name, scan_name=""):
        if not self.__is_nominal:
            sign_b1 = -1
            sign_b2 = 1
        else:
            sign_b1 = 1
            sign_b2 = -1

        if scan_name in self.__scans_with_inverted_beam_sign:
            print("Inverting sign for scan_name")
            sign_b1 *= -1
            sign_b2 *= -1

        if is_B1(col_name):
            col_sign = sign_b1
        else:
            col_sign = sign_b2

        return col_sign

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

        # if not self.__is_nominal:
        #     sign_b1 = -1
        #     sign_b2 = 1
        # else:
        #     sign_b1 = 1
        #     sign_b2 = -1

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
            scan_name = self.get_scan_name_from_timestamp(data_time)
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
                    # if is_B1(i_col):
                    #     correction *= sign_b1
                    # else:
                    #     correction *= sign_b2
                    correction *= self.get_beam_sign_for_column_name(i_col, scan_name)

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

        if not self.__is_nominal:
            # Nominal
            if self.__deflection_applied_in_nominal:
                self.__nominal_suffix_for_final_result += BPM.__ref_col_deflection_suffix_label
            if self.__lscale_applied_in_nominal:
                self.__nominal_suffix_for_final_result += BPM.__ref_col_lscale_suffix_label
            # Data
            if apply_orbit_drift:
                self.__data_suffix_for_final_result += BPM.__ref_col_orbit_drift_suffix_label
            if apply_lscale:
                self.__data_suffix_for_final_result += BPM.__ref_col_lscale_suffix_label

            self.__ref_col_diff_names_final_result = self.get_col_names(
                correction_suffix=self.__data_suffix_for_final_result,
                nominal_correction_suffix=self.__nominal_suffix_for_final_result, is_diff=True)
            self.__ref_col_diff_names_final_result_b1 = self.get_col_names(
                correction_suffix=self.__data_suffix_for_final_result,
                nominal_correction_suffix=self.__nominal_suffix_for_final_result, is_diff=True, only_b1=True)
            self.__ref_col_diff_names_final_result_b2 = self.get_col_names(
                correction_suffix=self.__data_suffix_for_final_result,
                nominal_correction_suffix=self.__nominal_suffix_for_final_result, is_diff=True, only_b2=True)
            self.__ref_col_diff_names_final_result_H = self.get_col_names(
                correction_suffix=self.__data_suffix_for_final_result,
                nominal_correction_suffix=self.__nominal_suffix_for_final_result, is_diff=True, only_H=True)
            self.__ref_col_diff_names_final_result_V = self.get_col_names(
                correction_suffix=self.__data_suffix_for_final_result,
                nominal_correction_suffix=self.__nominal_suffix_for_final_result, is_diff=True, only_V=True)

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

        # Getting H differences
        data_to_use[BPM.__col_H_sum] = data_to_use[matching_to_basename_dict[BPM.__col_b2h]] + \
                                       data_to_use[matching_to_basename_dict[BPM.__col_b1h]]
        # Getting V differences
        data_to_use[BPM.__col_V_sum] = data_to_use[matching_to_basename_dict[BPM.__col_b2v]] + \
                                       data_to_use[matching_to_basename_dict[BPM.__col_b1v]]

    def get_scan_name_from_timestamp(self, timestamp):
        scan_name = ""
        all_scans = list(self.__scans_timestamp_limits_dict)
        for i_scan in all_scans:
            i_limits = self.__scans_timestamp_limits_dict[i_scan]
            if i_limits[0] <= timestamp <= i_limits[1]:
                scan_name = i_scan
                break
        return scan_name

    def do_per_scan_studies(self):
        if self.check_for_boolean_setting(setts.conf_compute_length_scale_epsilon_fits):
            scan_names = list(self.__scans_timestamp_limits_dict)

            self.__data_per_scan_dict = ltools.split_pandas_dataset_from_col_values_ranges(
                self.__processed_data, BPM.__col_time, self.__scans_timestamp_limits_dict)

            self.__data_per_scan_nominal_dict = ltools.split_pandas_dataset_from_col_values_ranges(
                self.__nominal_data, BPM.__col_time, self.__scans_timestamp_limits_dict)

            ltools.check_and_create_folder(self.__output_dir + self.__output_per_scan_studies, creation_info=False)

            lin_model = Model(ltools.lin_func)

            for scan_name in scan_names:
                self.__lscale_per_beam_epsilon_dict[scan_name] = {}
                self.__lscale_per_beam_lscale_slope_dict[scan_name] = {}

                # result plots
                scan_title_label = "(scan " + scan_name + ") "
                fill_and_scan_label = str(self.__fill) + "(scan " + scan_name + ")"
                scan_data = self.__data_per_scan_dict[scan_name]
                scan_nominal = self.__data_per_scan_nominal_dict[scan_name]

                ltools.check_and_create_folder(self.__output_dir + self.__output_per_scan_studies + scan_name + "/",
                                               creation_info=False)
                self.plot_detector_data_diff_with_nominal(extra_name_suffix="_" + scan_name + "_final",
                                                          extra_name_prefix=self.__output_per_scan_studies +
                                                                            scan_name + "/",
                                                          cols_to_plot=self.__ref_col_diff_names_final_result,
                                                          data_to_use=scan_data,
                                                          canvas_square_shape=True,
                                                          marker_size=3.0,
                                                          leg_marker_scale=4,
                                                          extra_title=scan_title_label,
                                                          legend_labels=self.__ref_col_names)
                self.plot_detector_data_diff_with_nominal(extra_name_suffix="_" + scan_name + "_no_corrections",
                                                          extra_name_prefix=self.__output_per_scan_studies +
                                                                            scan_name + "/",
                                                          cols_to_plot=self.__ref_col_diff_names,
                                                          data_to_use=scan_data,
                                                          canvas_square_shape=True,
                                                          marker_size=3.0,
                                                          leg_marker_scale=4,
                                                          extra_title=scan_title_label)

                self.plot_detector_data_diff_with_nominal(extra_name_suffix="_" + scan_name + "_H_V_final",
                                                          extra_name_prefix=self.__output_per_scan_studies +
                                                                            scan_name + "/",
                                                          cols_to_plot=[BPM.__col_H_diff, BPM.__col_V_diff],
                                                          data_to_use=scan_data,
                                                          legend_labels=["H", "V"],
                                                          canvas_square_shape=True,
                                                          marker_size=3.0,
                                                          leg_marker_scale=4,
                                                          # automatic_y_range=True,
                                                          use_smaller_y_range=True,
                                                          extra_title=scan_title_label,
                                                          y_label="B2 - B1" +
                                                                  " [" + BPM.__distance_unit + "]")
                # self.plot_detector_data_diff_with_nominal(extra_name_suffix="_" + scan_name + "_H_V_no_corrections",
                #                                           extra_name_prefix=self.__output_per_scan_studies +
                #                                                             scan_name + "/",
                #                                           cols_to_plot=[BPM.__col_H_diff + "_ini",
                #                                                         BPM.__col_V_diff + "_ini"],
                #                                           data_to_use=scan_data,
                #                                           legend_labels=["H", "V"],
                #                                           canvas_square_shape=True,
                #                                           marker_size=3.0,
                #                                           leg_marker_scale=4,
                #                                           yrange=[-40, 40],
                #                                           extra_title=scan_title_label,
                #                                           y_label="B2 - B1" +
                #                                                   " [" + BPM.__distance_unit + "]")

                self.plot_detector_data_diff_with_nominal(extra_name_suffix="_" + scan_name + "_H_V_sum_final",
                                                          extra_name_prefix=self.__output_per_scan_studies +
                                                                            scan_name + "/",
                                                          cols_to_plot=[BPM.__col_H_sum, BPM.__col_V_sum],
                                                          data_to_use=scan_data,
                                                          legend_labels=["H", "V"],
                                                          canvas_square_shape=True,
                                                          marker_size=3.0,
                                                          leg_marker_scale=4,
                                                          extra_title=scan_title_label,
                                                          y_label="B2 + B1" +
                                                                  " [" + BPM.__distance_unit + "]")
                self.plot_detector_data_diff_with_nominal(extra_name_suffix="_" + scan_name + "_per_beam_final_H",
                                                          extra_name_prefix=self.__output_per_scan_studies +
                                                                            scan_name + "/",
                                                          cols_to_plot=self.__ref_col_diff_names_final_result_H,
                                                          data_to_use=scan_data,
                                                          canvas_square_shape=True,
                                                          marker_size=3.0,
                                                          leg_marker_scale=4,
                                                          extra_title=scan_title_label,
                                                          legend_labels=self.__cols_H_names)
                self.plot_detector_data_diff_with_nominal(extra_name_suffix="_" + scan_name + "_per_beam_final_V",
                                                          extra_name_prefix=self.__output_per_scan_studies +
                                                                            scan_name + "/",
                                                          cols_to_plot=self.__ref_col_diff_names_final_result_V,
                                                          data_to_use=scan_data,
                                                          canvas_square_shape=True,
                                                          marker_size=3.0,
                                                          leg_marker_scale=4,
                                                          extra_title=scan_title_label,
                                                          legend_labels=self.__cols_V_names)

                # if self.apply_lscale:
                #     self.plot_detector_data_diff_with_nominal(extra_name_suffix="_scale_corr",
                #                                               cols_to_plot=self.__ref_col_lscale_diff_names,
                #                                               cols_to_plot=self.__ref_col_deflection_and_lscale_diff_names,
                #                                               data_to_use=scan_data,
                #                                               canvas_square_shape=True,
                #                                               marker_size=3.0,
                #                                               leg_marker_scale=4,
                #                                               extra_title=scan_title_label,
                #                                               legend_labels=self.__ref_col_names)
                #     if self.__apply_deflection:
                #         self.plot_detector_data_diff_with_nominal(extra_name_suffix="_deflection_lscale_corr",
                #                                                   extra_name_prefix=self.__output_per_scan_studies +
                #                                                                     scan_name + "/",
                #                                                   cols_to_plot=self.__ref_col_deflection_and_lscale_diff_names,
                #                                                   data_to_use=scan_data,
                #                                                   canvas_square_shape=True,
                #                                                   marker_size=3.0,
                #                                                   leg_marker_scale=4,
                #                                                   extra_title=scan_title_label,
                #                                                   legend_labels=self.__ref_col_names)
                #     if self.__apply_orbit_drift:
                #         self.plot_detector_data_diff_with_nominal(extra_name_suffix="_orbit_drift",
                #                                                   cols_to_plot=self.__ref_col_orbit_drift_diff_names,
                #                                                   data_to_use=self.__processed_data,
                #                                                   legend_labels=self.__ref_col_names)


                # length scale epsilon fitting
                for i_col in BPM.__ref_col_names:
                    if self.check_if_scan_fit_is_meaningful_for_col_name(i_col, scan_name):
                        i_col_corrected_data = self.get_final_col_name_data(i_col)
                        i_col_corrected_nominal = self.get_final_col_name_nominal(i_col)
                        x = scan_nominal[i_col_corrected_nominal]
                        y = scan_data[i_col_corrected_data]

                        result_fit = lin_model.fit(y, x=x, a=1, b=1)

                        slope = result_fit.params['a'].value
                        # slope_err = result_fit.params['a'].stderr
                        # intercept = result_fit.params['b'].value
                        # intercept_err = result_fit.params['b'].stderr
                        # chi2 = result_fit.redchi

                        plot_name = self.__output_per_scan_studies + scan_name + "/" + i_col + \
                                    "_data_vs_nominal_linear_fit"
                        self.__plt_plots[plot_name] = plotting.plot_from_lmfit(fitted_model=result_fit,
                                                                               xlabel="Nominal " + i_col,
                                                                               ylabel=self.name + " " + i_col,
                                                                               title=fill_and_scan_label,
                                                                               title_loc="right")

                        # saving results into dict
                        self.__lscale_per_beam_lscale_slope_dict[scan_name][i_col] = slope
                        self.__lscale_per_beam_epsilon_dict[scan_name][i_col] = slope - 1
            self.process_epsilon_fit_results()

            # Saving dict into json file:
            output_json_path = self.__output_dir + self.__output_per_scan_studies + "lscale_per_beam_epsilon.json"
            output_slope_json_path = self.__output_dir + self.__output_per_scan_studies + "lscale_per_beam_slope.json"
            if not ltools.check_file_existence(output_json_path):
                ltools.color_print("\n ===>> Producing lscale_per_beam_epsilon.json because such file was not found",
                                   "yellow")
                ltools.save_dict_as_json(self.__lscale_per_beam_epsilon_dict, output_json_path)
                ltools.save_dict_as_json(self.__lscale_per_beam_lscale_slope_dict, output_slope_json_path)
                self.plot_per_beam_epsilon_fit_results(self.__output_per_scan_studies,
                                                       extra_file_suffix="_original")
            else:
                ltools.color_print("\n ===>> lscale_per_beam_epsilon.json not saved because such file was found. "
                                   "Please delete it if you want to produce a new file ;)",
                                   "green")
                self.plot_per_beam_epsilon_fit_results(self.__output_per_scan_studies)

    def check_if_scan_fit_is_meaningful_for_col_name(self, col_name, scan_name):
        # known_scan_labels_for_x = ["X1", "X2", "X3", "X4", "X5",
        #                            "L1X", "L2X"]
        # known_scan_labels_for_y = ["Y1", "Y2", "Y3", "Y4", "Y5",
        #                            "L1Y", "L2Y"]
        signature_label_x = "X"
        signature_label_y = "Y"

        is_meaningful = False
        if signature_label_x in scan_name:
            if is_H(col_name):
                is_meaningful = True
        elif signature_label_y in scan_name:
            if is_V(col_name):
                is_meaningful = True
        else:
            raise AssertionError("Scan direction not found!")

        return is_meaningful

    def plot_detector_all_data(self, save_data_to_file=True):
        plot_name = self.name + "_all_data"
        cols_to_plot = self.__cols_after_renaming
        cols_to_plot.remove(BPM.__col_time)

        if self.__debugging:
            ltools.color_print("plotting " + plot_name + "...", "green")

        ylabel = "beam position" + " [" + BPM.__distance_unit + "]"

        self.__plt_plots[plot_name] = plotting.scatter_plot_from_pandas_frame(self.__in_data_def_format,
                                                                              x_data_label=BPM.__col_time,
                                                                              y_data_label=cols_to_plot,
                                                                              xlabel="time",
                                                                              ylabel=ylabel, ncol_legend=4,
                                                                              plot_style='-',
                                                                              label_cms_status=False
                                                                              ).get_figure()
        self.__plt_plots[plot_name + "_min"] = plotting.scatter_plot_from_pandas_frame(self.__in_data_def_format,
                                                                                       x_data_label=BPM.__col_time_min,
                                                                                       y_data_label=cols_to_plot,
                                                                                       xlabel="time",
                                                                                       ylabel=ylabel, ncol_legend=4,
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

        if self.__debugging:
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

        ylabel = "beam position" + " [" + BPM.__distance_unit + "]"

        self.__plt_plots[plot_name] = plotting.scatter_plot_from_pandas_frame(data_to_use,
                                                                              x_data_label=BPM.__col_time_min,
                                                                              y_data_label=list(cols_to_plot),
                                                                              ymin=ymin,
                                                                              ymax=ymax,
                                                                              xmin=x_min, xmax=x_max,
                                                                              xlabel="time [min]",
                                                                              ylabel=ylabel,
                                                                              ncol_legend=4,
                                                                              plot_style='-',
                                                                              label_cms_status=False,
                                                                              draw_labels_pos_dict=draw_labels_pos_dict,
                                                                              draw_vertical_line_pos=draw_lines_scans_limits,
                                                                              title="Fill" + str(self.__fill),
                                                                              title_loc="right"
                                                                              ).get_figure()

    def plot_live_detector_data(self, xrange=None, extra_name_suffix="", cols_to_plot=None,
                                data_to_use=None, plot_scan_info=False, plot_scan_limits_lines=False,
                                scans_limits_to_use="time_range_minutes"):
        if data_to_use is None:
            data_to_use = self.__in_data_def_format

        if cols_to_plot is None:
            cols_to_plot = BPM.__ref_col_names

        ylabel = "beam position" + " [um]"
        xlabel = "time [min]"
        live_plotting.live_line_from_pandas(data_frame=data_to_use,
                                            x_data_label=BPM.__col_time_min,
                                            y_data_label=list(cols_to_plot),
                                            show_also_info_in=[BPM.__col_time],
                                            xlabel=xlabel, ylabel=ylabel
                                            )
        live_plotting.live_scatter_from_pandas(data_frame=data_to_use,
                                               x_data_label=BPM.__col_time_min,
                                               y_data_label=list(cols_to_plot),
                                               show_also_info_in=[BPM.__col_time],
                                               xlabel=xlabel, ylabel=ylabel)

    def plot_detector_data_diff_with_nominal(self, xrange=None, yrange=None,
                                             extra_name_suffix="", extra_name_prefix="",
                                             extra_title="",
                                             cols_to_plot=None, data_to_use=None, legend_labels=None,
                                             y_label=None, use_smaller_y_range=False, canvas_square_shape=False,
                                             automatic_y_range=False,
                                             marker_size=1.0, leg_marker_scale=7, plot_scan_info=False,
                                             plot_scan_limits_lines=False, scans_limits_to_use="time_range_minutes"):
        plot_name = extra_name_prefix + self.name + "_data_diff_with_nominal" + extra_name_suffix
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

        if self.__debugging:
            ltools.color_print("plotting " + plot_name, 'green')
            print("     using columns: " + str(cols_to_plot))

        if automatic_y_range:
            ymin = None
            ymax = None
        elif yrange is None:
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
        title_text = extra_title + "Fill" + str(self.__fill)
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
                                                                              draw_vertical_line_pos=draw_lines_scans_limits,
                                                                              title=title_text,
                                                                              title_loc="right"
                                                                              ).get_figure()

    def process_epsilon_fit_results(self):
        input_dict = self.__lscale_per_beam_epsilon_dict
        input_slope_dict = self.__lscale_per_beam_lscale_slope_dict
        scan_names = list(input_dict)
        scan_names_x = ltools.filter_string_list_with_substring(scan_names, "X")
        scan_names_y = ltools.filter_string_list_with_substring(scan_names, "Y")
        assert len(scan_names) == len(scan_names_x) + len(scan_names_y)

        slopes_x_b1 = []
        slopes_y_b1 = []
        slopes_x_b2 = []
        slopes_y_b2 = []

        epsi_x_b1 = []
        epsi_y_b1 = []
        epsi_x_b2 = []
        epsi_y_b2 = []

        if setts.conf_exclusion_from_length_scale_per_beam_epsilon in self.__settings_list:
            exclude_from_mean = self.__settings[setts.conf_exclusion_from_length_scale_per_beam_epsilon]
        else:
            exclude_from_mean = {}
        scans_with_exclusions = list(exclude_from_mean)

        for i_col in BPM.__ref_col_names:
            beam_name = get_beam_name(i_col)
            if is_H(i_col):
                scan_names_to_loop = scan_names_x
                dict_to_fill = self.__epsilon_x_in_plotting_structure_dict
                slope_dict_to_fill = self.__slope_x_in_plotting_structure_dict
                if is_B1(i_col):
                    slopes_to_fill = slopes_x_b1
                    epsi_to_fill = epsi_x_b1
                else:
                    slopes_to_fill = slopes_x_b2
                    epsi_to_fill = epsi_x_b2
            elif is_V(i_col):
                scan_names_to_loop = scan_names_y
                if is_B1(i_col):
                    slopes_to_fill = slopes_y_b1
                    epsi_to_fill = epsi_y_b1
                else:
                    slopes_to_fill = slopes_y_b2
                    epsi_to_fill = epsi_y_b2
                dict_to_fill = self.__epsilon_y_in_plotting_structure_dict
                slope_dict_to_fill = self.__slope_y_in_plotting_structure_dict
            else:
                raise AssertionError("Something wrong with columns looping")

            i_beam = get_beam_name(i_col)
            dict_to_fill[i_beam] = []
            slope_dict_to_fill[i_beam] = []
            for i_scan in scan_names_to_loop:
                dict_to_fill[i_beam].append(input_dict[i_scan][i_col])
                slope_dict_to_fill[i_beam].append(input_slope_dict[i_scan][i_col])
                if i_scan not in scans_with_exclusions or beam_name not in exclude_from_mean[i_scan]:
                    epsi_to_fill.append(input_dict[i_scan][i_col])
                    slopes_to_fill.append(input_slope_dict[i_scan][i_col])
                else:
                    print(i_scan, i_col, "excluded from per beam length scale mean")

        # print(slopes_x_b1, slopes_y_b1, slopes_x_b2, slopes_y_b2)
        slopes_x_b1 = np.array(slopes_x_b1)
        slopes_y_b1 = np.array(slopes_y_b1)
        slopes_x_b2 = np.array(slopes_x_b2)
        slopes_y_b2 = np.array(slopes_y_b2)

        epsi_x_b1 = np.array(epsi_x_b1)
        epsi_y_b1 = np.array(epsi_y_b1)
        epsi_x_b2 = np.array(epsi_x_b2)
        epsi_y_b2 = np.array(epsi_y_b2)

        self.__lscale_per_beam_slope_b1_x = slopes_x_b1.mean()
        self.__lscale_per_beam_slope_b1_y = slopes_y_b1.mean()
        self.__lscale_per_beam_slope_b2_x = slopes_x_b2.mean()
        self.__lscale_per_beam_slope_b2_y = slopes_y_b2.mean()

        self.__lscale_per_beam_epsilon_b1_x = epsi_x_b1.mean()
        self.__lscale_per_beam_epsilon_b1_y = epsi_y_b1.mean()
        self.__lscale_per_beam_epsilon_b2_x = epsi_x_b2.mean()
        self.__lscale_per_beam_epsilon_b2_y = epsi_y_b2.mean()

        self.__lscale_per_beam_slope_b1_x_err = slopes_x_b1.std()/len(slopes_x_b1)
        self.__lscale_per_beam_slope_b1_y_err = slopes_y_b1.std()/len(slopes_y_b1)
        self.__lscale_per_beam_slope_b2_x_err = slopes_x_b2.std()/len(slopes_x_b2)
        self.__lscale_per_beam_slope_b2_y_err = slopes_y_b2.std()/len(slopes_y_b2)

        self.__lscale_per_beam_epsilon_b1_x_err = epsi_x_b1.std()/len(epsi_x_b1)
        self.__lscale_per_beam_epsilon_b1_y_err = epsi_y_b1.std()/len(epsi_y_b1)
        self.__lscale_per_beam_epsilon_b2_x_err = epsi_x_b2.std()/len(epsi_x_b2)
        self.__lscale_per_beam_epsilon_b2_y_err = epsi_y_b2.std()/len(epsi_y_b2)

    def plot_per_beam_epsilon_fit_results(self, output_folder, extra_file_suffix=""):
        input_dict = self.__lscale_per_beam_epsilon_dict
        if len(input_dict) == 0:
            ltools.color_print("self.__lscale_per_beam_epsilon_dict is empty")
            raise AssertionError("self.__lscale_per_beam_epsilon_dict is empty")
        plots_base_name = "epsilon_fit_values_by_scans"
        plot_name_x = output_folder + "x_" + plots_base_name + extra_file_suffix
        plot_name_y = output_folder + "y_" + plots_base_name + extra_file_suffix
        scan_names = list(input_dict)
        scan_names_x = ltools.filter_string_list_with_substring(scan_names, "X")
        scan_names_y = ltools.filter_string_list_with_substring(scan_names, "Y")
        top_right_text = "Fill" + str(self.__fill)

        if setts.conf_limits_for_length_scale_per_beam_epsilon in self.__settings_list:
            slope_min = self.__settings[setts.conf_limits_for_length_scale_per_beam_epsilon][0]
            slope_max = self.__settings[setts.conf_limits_for_length_scale_per_beam_epsilon][1]
        else:
            slope_min = -0.05
            slope_max = 0.05

        summary_text_x = "B1 mean (error): " + str(float("{0:.4f}".format(self.__lscale_per_beam_epsilon_b1_x))) + " (" \
                         + str(float("{0:.4f}".format(self.__lscale_per_beam_epsilon_b1_x_err))) + ")\n" + \
                         "B2 mean (error): " + str(float("{0:.4f}".format(self.__lscale_per_beam_epsilon_b2_x))) + " (" \
                         + str(float("{0:.4f}".format(self.__lscale_per_beam_epsilon_b2_x_err))) + ")"
        self.__plt_plots[plot_name_x] = plotting.plot_scatter_from_dict(self.__epsilon_x_in_plotting_structure_dict,
                                                                        xlabel="Scans for X",
                                                                        ylabel=self.name + "/Nominal slope residual",
                                                                        new_xticks=scan_names_x,
                                                                        title=top_right_text,
                                                                        title_loc='right',
                                                                        ymin=slope_min,
                                                                        ymax=slope_max,
                                                                        summary_text=summary_text_x,
                                                                        h_line_1=self.__lscale_per_beam_epsilon_b1_x,
                                                                        h_line_1_err=self.__lscale_per_beam_epsilon_b1_x_err,
                                                                        h_line_2=self.__lscale_per_beam_epsilon_b2_x,
                                                                        h_line_2_err=self.__lscale_per_beam_epsilon_b2_x_err
                                                                        )
        summary_text_y = "B1 mean (error): " + str(float("{0:.4f}".format(self.__lscale_per_beam_epsilon_b1_y))) + " (" \
                         + str(float("{0:.4f}".format(self.__lscale_per_beam_epsilon_b1_y_err))) + ")\n" + \
                         "B2 mean (error): " + str(float("{0:.4f}".format(self.__lscale_per_beam_epsilon_b2_y))) + " (" \
                         + str(float("{0:.4f}".format(self.__lscale_per_beam_epsilon_b2_y_err))) + ")"
        self.__plt_plots[plot_name_y] = plotting.plot_scatter_from_dict(self.__epsilon_y_in_plotting_structure_dict,
                                                                        xlabel="Scans for Y",
                                                                        ylabel=self.name + "/Nominal slope residual",
                                                                        new_xticks=scan_names_y,
                                                                        title=top_right_text,
                                                                        title_loc='right',
                                                                        ymin=slope_min,
                                                                        ymax=slope_max,
                                                                        summary_text=summary_text_y,
                                                                        h_line_1=self.__lscale_per_beam_epsilon_b1_y,
                                                                        h_line_1_err=self.__lscale_per_beam_epsilon_b1_y_err,
                                                                        h_line_2=self.__lscale_per_beam_epsilon_b2_y,
                                                                        h_line_2_err=self.__lscale_per_beam_epsilon_b2_y_err
                                                                        )


    def plot_orbit_drift_result(self):
        data_frame = self.__data_in_zero_beam_position
        x_data_label = BPM.__col_time_min
        # y_data_label = [BPM.__col_H_diff, BPM.__col_V_diff]
        legend_labels = ["H", "V"]
        color_for_v = 'blue'
        color_for_h = 'green'
        # h_v_colors = [color_for_h, color_for_v]
        fig_size_shape = "nsq"
        title = self.name + " Orbit Drift (beam2 - beam1) " + " in Fill" + str(self.__fill)
        plot_file_name = "OrbitDrift_" + self.name
        ylabel = "Orbit Drift [" + BPM.__distance_unit + "]"
        ymin = self.__y_range_orbit_drift[0]
        ymax = self.__y_range_orbit_drift[1]
        xlabel = "time [min]"
        marker_size = 2.0

        fig_size = def_setts.fig_sizes[fig_size_shape]
        fig, ax = plt.subplots(figsize=fig_size)

        data_frame.plot(kind='scatter', x=x_data_label, y=BPM.__col_H_diff, s=marker_size, color=color_for_h, ax=ax)
        data_frame.plot(kind='scatter', x=x_data_label, y=BPM.__col_V_diff, s=marker_size, color=color_for_v, ax=ax)

        ncol_legend = len(legend_labels)
        leg_marker_sc = def_setts.leg_vs_plots_marker_scale
        legend_position = def_setts.leg_vs_plots_pos
        leg_text_s = def_setts.leg_vs_plots_text_s
        ax.legend(legend_labels, ncol=ncol_legend, markerscale=leg_marker_sc, fontsize=leg_text_s, loc=legend_position)

        ax.set_title(title)
        ax.set_ylabel(ylabel, labelpad=def_setts.axis_labelpad_y, weight=def_setts.axis_weight,
                      size=def_setts.axis_case_size[fig_size_shape])
        ax.set_xlabel(xlabel, labelpad=def_setts.axis_labelpad_x, weight=def_setts.axis_weight,
                      size=def_setts.axis_case_size[fig_size_shape])

        plt.xticks(fontsize=def_setts.axis_thicks_case_size[fig_size_shape])
        plt.yticks(fontsize=def_setts.axis_thicks_case_size[fig_size_shape])

        if ymin and ymax:
            plt.ylim(ymin, ymax)

        scans_names = list(self.__orbit_drifts_per_scan)
        y_pos_scans = ymin + abs(ymin * setts.delta_pos_for_scan_labels)
        draw_vertical_line_pos = []
        draw_labels_pos_dict = {}
        for scan in scans_names:
            i_x = self.__orbit_drifts_per_scan[scan]["TimeWindows_min"]
            i_y_h = self.__orbit_drifts_per_scan[scan]["OrbitDrifts_X"]
            i_y_v = self.__orbit_drifts_per_scan[scan]["OrbitDrifts_Y"]

            draw_vertical_line_pos.extend([i_x[0], i_x[2]])
            draw_labels_pos_dict[scan] = [i_x[1], y_pos_scans]

            # print(i_x, i_y_h, i_y_v)

            plt.plot(i_x, i_y_h, c=color_for_h)
            plt.plot(i_x, i_y_v, c=color_for_v)

        # Draw vertical lines
        if draw_vertical_line_pos is not None:
            for x_pos in draw_vertical_line_pos:
                plt.axvline(x=x_pos, linestyle='dashed', alpha=0.3)

        # Draw text from draw_labels_pos_dict -> {"text": [x, y], ... }
        if draw_labels_pos_dict is not None:
            labels_text = list(draw_labels_pos_dict)
            for i_label in labels_text:
                plt.text(draw_labels_pos_dict[i_label][0], draw_labels_pos_dict[i_label][1], i_label,
                         fontweight='bold', alpha=0.5, horizontalalignment='center', verticalalignment='center'
                         )

        if fig_size == (12, 4):
            plt.subplots_adjust(left=0.1, right=0.97, top=0.9, bottom=0.2)

        self.__plt_plots[plot_file_name] = fig

    def check_for_boolean_setting(self, setting_label_name):
        settled = False
        if setting_label_name in self.__settings_list:
            settled = self.__settings[setting_label_name]
        return settled

    def save_plots(self):
        ltools.color_print('\n\n Saving plots:', "green")
        plotting.save_plots(self.__plt_plots, self.__output_dir, save_pickle=setts.save_figures_as_pickle)
        plotting.save_plots(self.__sns_plots, self.__output_dir, save_pickle=setts.save_figures_as_pickle)

    @property
    def name(self):
        return self.__name

    @name.setter
    def name(self, name):
        self.__name = name

    @property
    def data(self):
        return self.__in_data_def_format

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

    @property
    def orbit_drifts_per_scan(self):
        return self.__orbit_drifts_per_scan

    @property
    def scans_info_dict(self):
        return self.__scans_info_dict

    @property
    def ref_col_names(self):
        return self.__ref_col_names

    @property
    def col_time_min(self):
        return self.__col_time_min

    @property
    def distance_unit(self):
        return self.__distance_unit

    @property
    def col_H_diff(self):
        return self.__col_H_diff

    @property
    def col_V_diff(self):
        return self.__col_V_diff

    @property
    def fill(self):
        return self.__fill

    @property
    def data_in_zero_beam_position(self):
        return self.__data_in_zero_beam_position

    @property
    def y_range_orbit_drift(self):
        return self.__y_range_orbit_drift

    @property
    def settings(self):
        return self.__settings

    @property
    def settings_list(self):
        return self.__settings_list


def get_base_name_from_DOROS_LR(full_name: str):
    if "_L" in full_name:
        return full_name.replace("_L", "")
    elif "_R" in full_name:
        return full_name.replace("_R", "")
    else:
        return full_name


def get_base_name(full_name: str):
    base_names = ["B1_H", "B1_V", "B2_H", "B2_V"]
    output_name = ""
    for base_name in base_names:
        if base_name in full_name:
            output_name = base_name
            break
    return output_name


def get_beam_name(full_name: str):
    base_names = ["B1", "B2"]
    beam_name = ""
    for base_name in base_names:
        if base_name in full_name:
            beam_name = base_name
            break
    return beam_name

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
