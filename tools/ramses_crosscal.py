import pandas as pd
from tools.luminometer import Luminometer as L
from scipy import interpolate
import numpy as np
import tools.plotting_tools as plotting
import tools.lumi_tools as ltools
import settings_ramses as setts_ram

# Cross calibration to reference detector luminosity units

class RamsesCrossCal:
    __csv_channel_file_col_names = ('time', 'dose')
    __cols_format = ("run:fill", "ls", "time", "beamstatus", "E(GeV)",
                     "delivered", "recorded", "avgpu", "source")
    __cols_to_fill_with_calibration = ("delivered", "recorded", "source")
    __ramses_name = 'RAMSES'
    __output_file_name = setts_ram.ramses_output_file_name
    __variation_studies_nbins = setts_ram.variation_studies_nbins
    __to_ref_studies_nbins = setts_ram.to_ref_studies_nbins
    __soft_ratio_lim = setts_ram.soft_ratio_lim
    __medium_ratio_lim = setts_ram.medium_ratio_lim
    __ratio_to_ref_min_rec = setts_ram.ratio_to_ref_min_rec
    __ratio_to_ref_max_rec = setts_ram.ratio_to_ref_max_rec
    __ratio_to_ref_min_del = setts_ram.ratio_to_ref_min_del
    __ratio_to_ref_max_del = setts_ram.ratio_to_ref_max_del
    __min_ratio = setts_ram.ramses_cal_min_ratio
    __max_ratio = setts_ram.ramses_cal_max_ratio

    # columns for uploading format and link to column name in data frame
    __upload_result_format = {
        'run': 'run',
        'ls': 'ls',
        'result_avglumi': 'normalized_to_ref_del'
        # 'result_bxlumi': 'delivered'
    }

    __comment_line = '#Data tag : cross_calibration , Norm tag: None \n' \
                     '#run:fill,ls,time,beamstatus,E(GeV),delivered(hz/ub),recorded(hz/ub),avgpu,source \n'

    def __init__(self, ramses_channels_plus_ref: list, remove_outliers_in_raw_channels: bool = False) -> None:
        if len(ramses_channels_plus_ref) >= 2:
            self.__output_dir = RamsesCrossCal.__output_file_name
            ltools.check_and_create_folder(self.__output_dir)
            self.__remove_outliers_in_raw_channels = remove_outliers_in_raw_channels
            self.__plt_plots = {}
            self.__var_plot_id = 0
            self.__ch_raw_plot_id = 0
            ramses_raw_files_path = []
            self.__channels_data = []

            ref_detector_file_path = ramses_channels_plus_ref[-1]
            ref_name = ref_detector_file_path.split("/")[-1].replace(".csv", "")
            self.__ref_det = L(ref_name, ref_detector_file_path, remove_extra_cols=False, split_run_fill_ls_cols=False)
            self.__label_lumi_ref_rec = self.__ref_det.lumi_rec_label_original_units
            self.__label_lumi_ref_del = self.__ref_det.lumi_del_label_original_units

            self.__label_normalized_to_ref_rec = 'normalized_to_ref_rec'
            self.__label_normalized_to_ref_del = 'normalized_to_ref_del'
            self.__label_ratio_to_ref_rec = 'ratio_to_ref_rec'
            self.__label_ratio_to_ref_del = 'ratio_to_ref_del'

            self.__all_cols_info_calibration = None
            # Only time, dose, __label_normalized_to_ref, __label_ratio_to_ref
            self.__calibrated_data = None

            for det_id in range(0, len(ramses_channels_plus_ref) - 1):
                ramses_raw_files_path.append(ramses_channels_plus_ref[det_id])
                self.__channels_data.append(self.read_channel_from_file(ramses_channels_plus_ref[det_id]))
            print("RAMSES channels files: " + str(ramses_raw_files_path))
            print("Taken " + str(ref_detector_file_path) + " as the reference detector")

            self.__number_of_channels = len(self.__channels_data)

            if self.__number_of_channels == 2:
                mean_ratio_raw_ch, std_ratio_raw_ch = self.get_ratio_between_channels(self.__channels_data[0],
                                                                                      self.__channels_data[1])

            # ramses_raw_files = []
            # for file_path in ramses_raw_files_path:
            #     try:
            #         ramses_raw_files.append(open(file_path))
            #     except IOError:
            #         print("Problem opening file: " + str(file_path))
            #         raise

            self.__channels_data_interpolated = []
            for channel_data in self.__channels_data:
                self.__channels_data_interpolated.append(
                    self.interpolate_to_ref_time(channel_df=channel_data, ref_df=self.__ref_det.data))

            if self.__number_of_channels == 2:
                mean_ratio_interp_ch, std_ratio_interp_ch = self.get_ratio_between_channels(self.__channels_data_interpolated[0],
                                                                                            self.__channels_data_interpolated[1],
                                                                                            extra_label="_interpolated")
            self.get_calibrated_data(mode=0)

            self.save_calibration_to_upload_format()

        else:
            raise AssertionError(
                "RAMSES cross calibration needs only 2 inputs: ramses_raw_file_path and ref_detector_file_path.")

    def merge_ramses_data(self, channels_list: list):
        print("Merging data from: " + str(channels_list))

    @staticmethod
    def interpolate_to_ref_time(channel_df, ref_df):
        print("Interpolating to all reference detector times.")
        tck = interpolate.splrep(channel_df['time'], channel_df['dose'], s=0)
        x = ref_df['time']
        y = interpolate.splev(x, tck, der=0)
        dict_interpolated = {'time': x, 'dose': y}

        return pd.DataFrame(dict_interpolated, columns=['time', 'dose'])

    def read_channel_from_file(self, file_name):
        channel_pd = pd.read_csv(file_name, comment='#', index_col=False)
        self.__ch_raw_plot_id += 1
        plot_raw_channel_data = plotting.hist_from_array(channel_pd['dose'],
                                                            nbins=RamsesCrossCal.__variation_studies_nbins,
                                                            xlabel='dose', ylabel='counts',
                                                            fig_size_shape='sq')
        plotting.save_py_fig_to_file(plot_raw_channel_data,
                                     self.__output_file_name + 'raw_channel_data_' + str(self.__ch_raw_plot_id))
        channel_pd = self.clean_channel_df(channel_pd)
        return channel_pd

    def clean_channel_df(self, df: pd.DataFrame):
        df.dropna(inplace=True)
        df = df[df['dose'] >= 0.0]
        df = df.reset_index(drop=True)
        initial_no_nan_size = len(df)
        # compute variation mean value
        if self.__remove_outliers_in_raw_channels:
            self.__var_plot_id += 1
            print("\n Removing outliers in raw channel data:")
            variations_array = []
            for index_data in range(1, len(df)):
                variations_array.append(abs(df['dose'][index_data] - df['dose'][index_data - 1]))
            variations_np = np.array(variations_array)
            deviation_mean = variations_np.mean()
            deviation_stdv = variations_np.std()
            dif_th = deviation_mean + 2*deviation_stdv
            print("Removing outliers using mean=" + str(deviation_mean) + " ,with std=" + str(deviation_stdv))
            plot_variation_diff_hist = plotting.hist_from_array(variations_array, xmin=0.0, xmax=deviation_mean + 4*deviation_stdv,
                                                                nbins=RamsesCrossCal.__variation_studies_nbins, xlabel='diffs', ylabel='counts',
                                                                mean=deviation_mean, stdv=deviation_stdv,
                                                                fig_size_shape='sq')
            plotting.save_py_fig_to_file(plot_variation_diff_hist, self.__output_file_name + 'variation_diff_hist_' + str(self.__var_plot_id))
            # excl_size = 0
            # for index_data in range(1, len(df)):
            #     dif = abs(df['dose'][index_data] - df['dose'][index_data - 1])
            #     if dif > dif_th:
            #         excl_size += 1
            #         df.drop(index_data, inplace=True, axis=0)
            # print("Remaining percent of data: " + str((initial_no_nan_size - excl_size)*100./initial_no_nan_size))
        return df

    def get_calibrated_data(self, mode=0):
        if mode == 0:
            print(
                "\n Normalization procedure: mean value between channels -> cross calibrate mean-channels detector to reference \n")
            if self.__number_of_channels == 2:
                raw_data = self.get_mean_channel_detector_df()
            elif self.__number_of_channels == 1:
                raw_data = self.__channels_data_interpolated[0]
            else:
                raise AssertionError("Number of channels not implemented")
            self.normalized_to_ref(raw_data)
        else:
            raise AssertionError("mode not implemented")

    def get_mean_channel_detector_df(self):
        array_means = []
        for df_channel in self.__channels_data_interpolated:
            array_means.append(df_channel['dose'])
        np_array = np.array(array_means)
        row_means = np.mean(np_array, axis=0)

        assert (len(self.__channels_data_interpolated[0]) == len(row_means))

        return pd.DataFrame({'time': self.__channels_data_interpolated[0]['time'], 'dose': row_means})

    def normalized_to_ref(self, data):
        ref_data = self.__ref_det.data

        ratio_data_ref_del = data['dose'] / ref_data[self.__label_lumi_ref_del]
        ratio_data_ref_del.dropna(inplace=True)
        ratio_data_ref_del = ratio_data_ref_del[(ratio_data_ref_del < RamsesCrossCal.__ratio_to_ref_max_del) &
                                                (ratio_data_ref_del > RamsesCrossCal.__ratio_to_ref_min_del)]
        mean_del = np.mean(ratio_data_ref_del)
        std_del = np.std(ratio_data_ref_del)

        ratio_data_ref_rec = data['dose'] / ref_data[self.__label_lumi_ref_rec]
        ratio_data_ref_rec.dropna(inplace=True)
        ratio_data_ref_rec = ratio_data_ref_rec[(ratio_data_ref_rec < RamsesCrossCal.__ratio_to_ref_max_rec) &
                                                (ratio_data_ref_rec > RamsesCrossCal.__ratio_to_ref_min_rec)]
        mean_rec = np.mean(ratio_data_ref_rec)
        std_rec = np.std(ratio_data_ref_rec)

        print("Normalizing to mean (delivered)=" + str(mean_del) + " ,with std=" + str(std_del))
        print("Using reference data: " + str(self.__label_lumi_ref_del))
        print("Normalizing to mean (recorded)=" + str(mean_rec) + " ,with std=" + str(std_rec))
        print("Using reference data: " + str(self.__label_lumi_ref_rec))

        plot_ratio_to_ref_rec = plotting.hist_from_array(np.array(ratio_data_ref_rec),
                                                         nbins=RamsesCrossCal.__to_ref_studies_nbins,
                                                         xlabel='uncalibrated detector/reference detector', ylabel='counts',
                                                         mean=mean_rec/setts_ram.ramses_scale,
                                                         stdv=std_rec/setts_ram.ramses_scale,
                                                         fig_size_shape='sq')
        plotting.save_py_fig_to_file(plot_ratio_to_ref_rec,
                                     self.__output_file_name + 'ratio_to_ref_recorded')

        plot_ratio_to_ref_del = plotting.hist_from_array(np.array(ratio_data_ref_del),
                                                         nbins=RamsesCrossCal.__to_ref_studies_nbins,
                                                         xlabel='uncalibrated detector/reference detector',
                                                         ylabel='counts',
                                                         mean=mean_del / setts_ram.ramses_scale,
                                                         stdv=std_del / setts_ram.ramses_scale,
                                                         fig_size_shape='sq')
        plotting.save_py_fig_to_file(plot_ratio_to_ref_del,
                                     self.__output_file_name + 'ratio_to_ref_delivered')

        data[self.__label_normalized_to_ref_rec] = data['dose']/mean_rec
        data[self.__label_ratio_to_ref_rec] = data[self.__label_normalized_to_ref_rec]/ref_data[self.__label_lumi_ref_rec]

        data[self.__label_normalized_to_ref_del] = data['dose'] / mean_del
        data[self.__label_ratio_to_ref_del] = data[self.__label_normalized_to_ref_del] / ref_data[
            self.__label_lumi_ref_del]

        ratio_norm_array_rec = np.array(data[self.__label_ratio_to_ref_rec].copy())
        ratio_norm_array_del = np.array(data[self.__label_ratio_to_ref_del].copy())
        ratio_norm_array_rec = ratio_norm_array_rec[(ratio_norm_array_rec < RamsesCrossCal.__max_ratio) &
                                                    (ratio_norm_array_rec > RamsesCrossCal.__min_ratio)]
        ratio_norm_array_del = ratio_norm_array_del[(ratio_norm_array_del < RamsesCrossCal.__max_ratio) &
                                                    (ratio_norm_array_del > RamsesCrossCal.__min_ratio)]
        mean_norm_rec = np.mean(ratio_norm_array_rec)
        mean_norm_del = np.mean(ratio_norm_array_del)
        std_norm_rec = np.std(ratio_norm_array_rec)
        std_norm_del = np.std(ratio_norm_array_del)

        plot_ratio_to_ref_norm_rec = plotting.hist_from_array(ratio_norm_array_rec,
                                                              nbins=RamsesCrossCal.__to_ref_studies_nbins,
                                                              xlabel='normalized detector/reference detector', ylabel='counts',
                                                              mean=mean_norm_rec, stdv=std_norm_rec,
                                                              fig_size_shape='sq')
        plotting.save_py_fig_to_file(plot_ratio_to_ref_norm_rec,
                                     self.__output_file_name + 'ratio_to_ref_norm_recorded')

        plot_ratio_to_ref_norm_del = plotting.hist_from_array(ratio_norm_array_del,
                                                              nbins=RamsesCrossCal.__to_ref_studies_nbins,
                                                              xlabel='normalized detector/reference detector',
                                                              ylabel='counts',
                                                              mean=mean_norm_del, stdv=std_norm_del,
                                                              fig_size_shape='sq')
        plotting.save_py_fig_to_file(plot_ratio_to_ref_norm_del,
                                     self.__output_file_name + 'ratio_to_ref_norm_delivered')

        # plot_ratio_to_ref_norm_csv = plotting.hist_from_array(data[self.__label_ratio_to_ref],
        #                                                       nbins=RamsesCrossCal.__to_ref_studies_nbins,
        #                                                       xlabel='normalized detector/reference detector',
        #                                                       ylabel='counts',
        #                                                       xmin=RamsesCrossCal.__min_ratio,
        #                                                       xmax=RamsesCrossCal.__max_ratio,
        #                                                       mean=mean_norm, stdv=std_norm,
        #                                                       fig_size_shape='sq')
        # plotting.save_py_fig_to_file(plot_ratio_to_ref_norm_csv,
        #                              self.__output_file_name + 'ratio_to_ref_norm_csv')

        self.__calibrated_data = data
        ## setting extra columns for calibrated data
        run_fill_cols = self.__ref_det.data["run:fill"].str.split(":", n=1, expand=True)
        self.__calibrated_data["run"] = run_fill_cols[0].astype(str).astype(int)
        self.__calibrated_data["fill"] = run_fill_cols[1].astype(str).astype(int)
        ls_double = self.__ref_det.data["ls"].str.split(":", n=1, expand=True)
        self.__calibrated_data["ls"] = ls_double[0].astype(str).astype(int)


        self.__all_cols_info_calibration = pd.DataFrame()

        for col in RamsesCrossCal.__cols_format:
            if col not in RamsesCrossCal.__cols_to_fill_with_calibration:
                self.__all_cols_info_calibration[col] = self.__ref_det.data[col]
            else:
                if col == 'source':
                    self.__all_cols_info_calibration[col] = RamsesCrossCal.__ramses_name
                elif col == "delivered":
                    self.__all_cols_info_calibration[col] = data[self.__label_normalized_to_ref_del]
                elif col == "recorded":
                    self.__all_cols_info_calibration[col] = data[self.__label_normalized_to_ref_rec]
                else:
                    raise BrokenPipeError("something wrong!")

        self.__all_cols_info_calibration.dropna(inplace=True)
        print('\n Output summary: \n')
        print(self.__all_cols_info_calibration)

        saving_file = self.__output_dir + "ramses_calibrated.csv"
        with open(saving_file, 'w') as f:
            f.write(RamsesCrossCal.__comment_line)
        self.__all_cols_info_calibration.to_csv(saving_file, index=False, header=False, mode='a')

    def get_ratio_between_channels(self, ch1_data, ch2_data, extra_label=""):
        ratio_data_ref = ch1_data['dose'] / ch2_data['dose']
        ratio_data_ref.dropna(inplace=True)
        no_nan_size = len(ratio_data_ref)
        ratio_data_ref_soft = ratio_data_ref[(ratio_data_ref < RamsesCrossCal.__soft_ratio_lim) &
                                             (ratio_data_ref > 0.0)]
        ratio_data_ref_med = ratio_data_ref[(ratio_data_ref < RamsesCrossCal.__medium_ratio_lim) &
                                            (ratio_data_ref > 0.0)]
        after_soft_ratio_lim = len(ratio_data_ref_soft)
        after_med_ratio_lim = len(ratio_data_ref_med)
        # ratio_data_ref[ratio_data_ref == np.inf] = 0.0
        mean_med = np.mean(ratio_data_ref_med)
        std_med = np.std(ratio_data_ref_med)
        mean_soft = np.mean(ratio_data_ref_soft)
        std_soft = np.std(ratio_data_ref_soft)
        print(mean_med, std_med, np.max(ratio_data_ref_med), (no_nan_size-after_med_ratio_lim)*100./no_nan_size)
        plot_ratio_channels_med = plotting.hist_from_array(np.array(ratio_data_ref_med),
                                                           nbins=RamsesCrossCal.__variation_studies_nbins,
                                                           xlabel='ch1/ch2 dose ratio', ylabel='counts',
                                                           mean=mean_med, stdv=std_med,
                                                           fig_size_shape='sq')
        plot_ratio_channels_soft = plotting.hist_from_array(np.array(ratio_data_ref_soft),
                                                            nbins=RamsesCrossCal.__variation_studies_nbins,
                                                            xlabel='ch1/ch2 dose ratio', ylabel='counts',
                                                            mean=mean_soft, stdv=std_soft,
                                                            fig_size_shape='sq')
        plotting.save_py_fig_to_file(plot_ratio_channels_med,
                                     self.__output_file_name + 'ratio_channels_1_2_med' + extra_label)
        plotting.save_py_fig_to_file(plot_ratio_channels_soft,
                                     self.__output_file_name + 'ratio_channels_1_2_soft' + extra_label)
        return mean_med, std_med

    def save_calibration_to_upload_format(self):
        data_to_upload = pd.DataFrame()
        labels_dict = RamsesCrossCal.__upload_result_format

        for up_col_name in list(labels_dict):
            data_to_upload[up_col_name] = self.__calibrated_data[labels_dict[up_col_name]]

        # Save into file
        to_upload_file = self.__output_dir + "data_to_upload.csv"
        data_to_upload.to_csv(to_upload_file, index=False, header=False, mode='w')

