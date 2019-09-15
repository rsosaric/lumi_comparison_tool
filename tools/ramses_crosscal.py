import pandas as pd
from tools.luminometer import Luminometer as L
from scipy import interpolate
import numpy as np
import tools.plotting_tools as plotting
import tools.lumi_tools as ltools
import settings as setts

# Cross calibration to reference detector luminosity units

class RamsesCrossCal:
    __csv_channel_file_col_names = ('time', 'dose')
    __output_file_name = setts.ramses_output_file_name
    __variation_studies_nbins = setts.variation_studies_nbins
    __to_ref_studies_nbins = setts.to_ref_studies_nbins
    __soft_ratio_lim = setts.soft_ratio_lim
    __medium_ratio_lim = setts.medium_ratio_lim
    __ratio_to_ref_min = setts.ratio_to_ref_min
    __ratio_to_ref_max = setts.ratio_to_ref_max
    __min_ratio = setts.ramses_cal_min_ratio
    __max_ratio = setts.ramses_cal_max_ratio

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
            self.__ref_det = L(ref_name, ref_detector_file_path)
            self.__label_normalized_to_ref = 'normalized_to_ref'
            self.__label_ratio_to_ref = 'ratio_to_ref'

            for det_id in range(0, len(ramses_channels_plus_ref) - 1):
                ramses_raw_files_path.append(ramses_channels_plus_ref[det_id])
                self.__channels_data.append(self.read_channel_from_file(ramses_channels_plus_ref[det_id]))
            print("RAMSES channels files: " + str(ramses_raw_files_path))
            print("Taken " + str(ref_detector_file_path) + " as the reference detector")

            if len(self.__channels_data) == 2:
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

            if len(self.__channels_data) == 2:
                mean_ratio_interp_ch, std_ratio_interp_ch = self.get_ratio_between_channels(self.__channels_data_interpolated[0],
                                                                                            self.__channels_data_interpolated[1],
                                                                                            extra_label="_interpolated")

            self.get_calibrated_data(mode=0)
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
            mean_channels_detector_df = self.get_mean_channel_detector_df()
            self.normalized_to_ref(mean_channels_detector_df)
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
        ratio_data_ref = data['dose'] / ref_data[self.__ref_det.lumi_del_label]
        ratio_data_ref.dropna(inplace=True)

        ratio_data_ref = ratio_data_ref[(ratio_data_ref < RamsesCrossCal.__ratio_to_ref_max) &
                                        (ratio_data_ref > RamsesCrossCal.__ratio_to_ref_min)]

        mean = np.mean(ratio_data_ref)
        std = np.std(ratio_data_ref)

        print("Normalizing to mean=" + str(mean) + " ,with std=" + str(std))

        plot_ratio_to_ref = plotting.hist_from_array(np.array(ratio_data_ref),
                                                     nbins=RamsesCrossCal.__to_ref_studies_nbins,
                                                     xlabel='uncalibrated detector/reference detector', ylabel='counts',
                                                     mean=mean, stdv=std,
                                                     fig_size_shape='sq')
        plotting.save_py_fig_to_file(plot_ratio_to_ref,
                                     self.__output_file_name + 'ratio_to_ref')
        data[self.__label_normalized_to_ref] = data['dose']/mean
        data[self.__label_ratio_to_ref] = data[self.__label_normalized_to_ref]/ref_data[self.__ref_det.lumi_del_label]

        ratio_norm_array = np.array(data[self.__label_ratio_to_ref])
        ratio_norm_array = ratio_norm_array[(ratio_norm_array < RamsesCrossCal.__max_ratio) &
                                            (ratio_norm_array > RamsesCrossCal.__min_ratio)]
        mean_norm = np.mean(ratio_norm_array)
        std_norm = np.std(ratio_norm_array)

        plot_ratio_to_ref_norm = plotting.hist_from_array(ratio_norm_array,
                                                          nbins=RamsesCrossCal.__to_ref_studies_nbins,
                                                          xlabel='normalized detector/reference detector', ylabel='counts',
                                                          mean=mean_norm, stdv=std_norm,
                                                          fig_size_shape='sq')
        plotting.save_py_fig_to_file(plot_ratio_to_ref_norm,
                                     self.__output_file_name + 'ratio_to_ref_norm')








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