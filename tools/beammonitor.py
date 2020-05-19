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

    __cols_to_be_position_scaled = (__col_b1h, __col_b1v, __col_b2h, __col_b2v,
                                    __col_b1hR, __col_b1vR, __col_b2hR, __col_b2vR,
                                    __col_b1hL, __col_b1vL, __col_b2hL, __col_b2vL)

    __ref_col_names = (__col_b1h, __col_b1v, __col_b2h, __col_b2v)
    __ref_col_diff_suffix_label = '_diff'

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

        if self.name == BPM.__nominal_name:
            self.__is_nominal = True
            self.__interpolation_to_nominal_time = None
            self.__ref_col_diff_names = []
        else:
            self.__is_nominal = False
            self.__ref_col_diff_names = [s + BPM.__ref_col_diff_suffix_label for s in BPM.__ref_col_names]
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
        # print(self.__in_data_def_format)

        if not self.__is_nominal:
            self.get_interpolation_to_nominal_time()
            self.get_nominal_diff()
            self.add_time_min_col(self.__interpolation_to_nominal_time)

        # Plotting
        if not get_only_data:
            self.plot_detector_data()
            if setts.conf_label_special_time_intervals in list(setts.config_dict[self.name][self.__fill]):
                for xrange in setts.config_dict[self.name][self.__fill][setts.conf_label_special_time_intervals]:
                    self.plot_detector_data(xrange=xrange)
            if self.name != BPM.__nominal_name:
                self.plot_detector_all_data()
                self.plot_detector_data_diff_with_nominal()
                if setts.conf_label_special_time_intervals in list(setts.config_dict[self.name][self.__fill]):
                    for xrange in setts.config_dict[self.name][self.__fill][setts.conf_label_special_time_intervals]:
                        self.plot_detector_data_diff_with_nominal(xrange=xrange)

            self.save_plots()

        #   Saving data to file
        cols_to_save = [BPM.__col_time, BPM.__col_time_min] + list(BPM.__ref_col_names) + self.__ref_col_diff_names
        if self.__is_nominal:
            data_to_save = self.__in_data_def_format
        else:
            data_to_save = self.__interpolation_to_nominal_time
        ltools.save_columns_from_pandas_to_file(data_to_save, cols_to_save, self.__output_dir +
                                                "plotted_data.csv")

    @property
    def name(self):
        return self.__name

    @name.setter
    def name(self, name):
        self.__name = name

    @property
    def data(self):
        return self.__in_data_def_format

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
                    data_in_unified_format[BPM.__col_time_min] = (data_in_unified_format[BPM.__col_time] - time0) * BPM.__scale_time
                elif i_col in BPM.__cols_to_be_position_scaled:
                    rescale_factor = BPM.__scale_position[self.name]
                    data_in_unified_format[i_col] = data_in_unified_format[i_col] * rescale_factor
                else:
                    print("     *** WARNING -->> column " + i_col + " not scaled!")

        if self.name == BPM.__doros_name:
            offsets_time = setts.config_dict[self.name][self.__fill][setts.conf_label_offset_time]
            offsets_values = setts.config_dict[self.name][self.__fill][setts.conf_label_offset_values]
            n_offsets_intervals = len(offsets_time)

            # Check if arrays are ordered in time
            if not sorted(offsets_time) == offsets_time:
                raise AssertionError("Offsets values most be ordered in time")

            col_b1h = []
            col_b1v = []
            col_b2h = []
            col_b2v = []

            # read offsets:
            pos_in_offsets_time_array = 0
            print("Applying offsets for time: " + str(offsets_time[pos_in_offsets_time_array]))
            off_b1h = offsets_values[0][setts.col_pos_offset_array[BPM.__col_b1h]]
            off_b1v = offsets_values[0][setts.col_pos_offset_array[BPM.__col_b1v]]
            off_b2h = offsets_values[0][setts.col_pos_offset_array[BPM.__col_b2h]]
            off_b2v = offsets_values[0][setts.col_pos_offset_array[BPM.__col_b2v]]

            for i_iter in range(0, len(data_in_unified_format)):

                if pos_in_offsets_time_array < n_offsets_intervals - 1 and \
                        data_in_unified_format[BPM.__col_time][i_iter] >= offsets_time[pos_in_offsets_time_array+1]:
                    pos_in_offsets_time_array += 1
                    print("Applying offsets for time: " + str(offsets_time[pos_in_offsets_time_array]))
                    off_b1h = offsets_values[pos_in_offsets_time_array][setts.col_pos_offset_array[BPM.__col_b1h]]
                    off_b1v = offsets_values[pos_in_offsets_time_array][setts.col_pos_offset_array[BPM.__col_b1v]]
                    off_b2h = offsets_values[pos_in_offsets_time_array][setts.col_pos_offset_array[BPM.__col_b2h]]
                    off_b2v = offsets_values[pos_in_offsets_time_array][setts.col_pos_offset_array[BPM.__col_b2v]]

                mean_value_LR_b1h = (data_in_unified_format[BPM.__col_b1hL][i_iter] +
                                     data_in_unified_format[BPM.__col_b1hR][i_iter]) / 2
                mean_value_LR_b1v = (data_in_unified_format[BPM.__col_b1vL][i_iter] +
                                     data_in_unified_format[BPM.__col_b1vR][i_iter]) / 2
                mean_value_LR_b2h = (data_in_unified_format[BPM.__col_b2hL][i_iter] +
                                     data_in_unified_format[BPM.__col_b2hR][i_iter]) / 2
                mean_value_LR_b2v = (data_in_unified_format[BPM.__col_b2vL][i_iter] +
                                     data_in_unified_format[BPM.__col_b2vR][i_iter]) / 2

                col_b1h.append(mean_value_LR_b1h - off_b1h)
                col_b1v.append(mean_value_LR_b1v - off_b1v)
                col_b2h.append(mean_value_LR_b2h - off_b2h)
                col_b2v.append(mean_value_LR_b2v - off_b2v)

            data_in_unified_format[BPM.__col_b1h] = np.array(col_b1h)
            data_in_unified_format[BPM.__col_b1v] = np.array(col_b1v)
            data_in_unified_format[BPM.__col_b2h] = np.array(col_b2h)
            data_in_unified_format[BPM.__col_b2v] = np.array(col_b2v)

        return data_in_unified_format

    def get_nominal_diff(self):
        in_data = self.__interpolation_to_nominal_time
        for col_name in BPM.__ref_col_names:
            in_data[col_name + BPM.__ref_col_diff_suffix_label] = in_data[col_name] - self.__nominal_data[col_name]

    def get_interpolation_to_nominal_time(self):
        in_data = self.__in_data_def_format
        self.__interpolation_to_nominal_time = pd.DataFrame()
        x = self.__nominal_data[BPM.__col_time]
        self.__interpolation_to_nominal_time[BPM.__col_time] = x
        for col_name in BPM.__ref_col_names:
            tck = interpolate.splrep(in_data[BPM.__col_time], in_data[col_name], s=0)
            y = interpolate.splev(x, tck, der=0)
            self.__interpolation_to_nominal_time[col_name] = np.array(y)

    def add_time_min_col(self, in_df):
        if BPM.__col_time in in_df.columns:
            time0 = self.__zero_time
            in_df[BPM.__col_time_min] = (in_df[BPM.__col_time] - time0) * BPM.__scale_time
        else:
            raise AssertionError("No time column available")

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

    def plot_detector_data(self, xrange=None, save_data_to_file=True):
        plot_name = self.name + "_data"
        x_min = None
        x_max = None
        if xrange:
            x_min = xrange[0]
            x_max = xrange[1]
            plot_name += "_" + str(x_min) + "_" + str(x_max)
        print("plotting " + plot_name + "...")

        cols_to_plot = BPM.__ref_col_names
        self.__plt_plots[plot_name] = plotting.scatter_plot_from_pandas_frame(self.__in_data_def_format,
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

    def plot_detector_data_diff_with_nominal(self, xrange=None, save_data_to_file=True):
        plot_name = self.name + "_data_diff_with_nominal"
        x_min = None
        x_max = None
        if xrange:
            x_min = xrange[0]
            x_max = xrange[1]
            plot_name += "_" + str(x_min) + "_" + str(x_max)
        print("plotting " + plot_name + "...")
        cols_to_plot = self.__ref_col_diff_names

        self.__plt_plots[plot_name] = plotting.scatter_plot_from_pandas_frame(self.__interpolation_to_nominal_time,
                                                                              x_data_label=BPM.__col_time_min,
                                                                              y_data_label=list(cols_to_plot),
                                                                              ymin=self.__y_diff_range[0],
                                                                              ymax=self.__y_diff_range[1],
                                                                              xmin=x_min, xmax=x_max,
                                                                              xlabel="time [min]",
                                                                              ylabel= self.name + " - Nominal diff." + " [" +
                                                                                     BPM.__distance_unit + "]",
                                                                              ncol_legend=4,
                                                                              plot_style='o',
                                                                              label_cms_status=False
                                                                              ).get_figure()

    def save_plots(self):
        print('\n\n Saving plots:')
        plotting.save_plots(self.__plt_plots, self.__output_dir, save_pickle=setts.save_figures_as_pickle)
        plotting.save_plots(self.__sns_plots, self.__output_dir, save_pickle=setts.save_figures_as_pickle)