from tools.luminometer import Luminometer as L
import numpy as np
import pandas as pd
import settings as setts
from tools import lumi_tools as ltools
import tools.plotting_tools as plotting


class DetectorsRatio(L):

    def __init__(self, det1: L, det2: L, year: str = None, energy: str = None,
                 fill_nls_data: bool = True, fill_stats=True, load_all_data = False, c_years = False,
                 nls = None) -> None:

        self.__detcs = (det1, det2)

        if fill_stats:
            fill_nls_data = True

        if det1.all_data_analysis_included and det1.all_data_analysis_included:
            self.__all_data_analysis_included = True
        else:
            self.__all_data_analysis_included = False

        if year is None:
            if (det1.year == det2.year) and (det1.energy == det2.energy):
                self.year = det1.year
            else:
                raise ValueError("Comparing detectors from different years or energies makes no sense :(")
        else:
            if c_years:
                list_years = year.split(',')
                self.year = list_years[0]
            else:
                self.year = year

        if energy is None:
            if (det1.year == det2.year) and (det1.energy == det2.energy):
                self.energy = det1.energy
            else:
                raise ValueError("Comparing detectors from different years or energies makes no sense :(")
        else:
            self.energy = energy

        if nls:
            self.__nls = nls
        else:
            try:
                self.__nls = setts.nls_year[int(self.year), int(self.energy)]
            except Warning:
                print("** Warning: no Nls configuration found for " + str(self.year) + " [" + str(self.energy) + "tev]")
                print("Using default value from settings: " + str(setts.nls_default) + "nls")
                self.__nls = setts.nls_default

        # column labels for dataframes
        if c_years:
            self.output_dir = setts.default_output_dir + year.replace(',', '-') + '/' + det1.name + '-' + det2.name + '/'
            self.__year_energy_label = str(self.energy) + 'TeV(' + year + ')'
        else:
            self.output_dir = setts.default_output_dir + str(self.year) + '/' + det1.name + '-' + det2.name + '/'
            self.__year_energy_label = str(self.energy) + 'TeV(' + str(self.year) + ')'

        self.__label_ratio = det1.name + '/' + det2.name
        self.__by_nls_label_ratio = 'by_' + str(self.__nls) + 'nls_' + self.__label_ratio
        self.__lumi_unit = det2.lumi_unit
        self.__det1 = det1
        self.__det2 = det2
        self.__by_nls_lumi_rec_label1 = det1.lumi_rec_label + '_by_' + str(self.__nls) + 'nls'
        self.__by_nls_lumi_rec_label2 = det2.lumi_rec_label + '_by_' + str(self.__nls) + 'nls'
        self.__by_nls_label_ratio_err = self.__by_nls_label_ratio + '_err'
        self.__by_nls_lumi_label = 'by_' + str(self.__nls) + 'nls_lumi_' + self.__label_ratio
        self.__accumulated_rec_lumi1_label = det1.lumi_rec_label + '_accumulated'
        self.__accumulated_rec_lumi2_label = det2.lumi_rec_label + '_accumulated'
        self.__by_nls_accumulated_rec_lumi1_label = self.__by_nls_lumi_rec_label1 + '_accumulated'
        self.__by_nls_accumulated_rec_lumi2_label = self.__by_nls_lumi_rec_label2 + '_accumulated'
        self.__rec_label = det2.lumi_rec_label
        self.__ratio_excluded_label = 'Exclusion info (' + self.__label_ratio + ')'

        keys_for_merging = ['ls', 'time', 'run', 'fill']
        print ("Merging " + det1.name + " and " + det2.name + " ...")
        common_data_in = pd.merge(det1.data, det2.data, on=keys_for_merging, how='outer')

        # Computing single ls ratio
        common_data_in[self.label_ratio] = common_data_in[det1.lumi_rec_label] / common_data_in[det2.lumi_rec_label]

        # Adding cumulative luminosity column
        common_data_in[self.accumulated_rec_lumi2_label] = common_data_in[det2.lumi_rec_label].cumsum()
        common_data_in[self.accumulated_rec_lumi1_label] = common_data_in[det1.lumi_rec_label].cumsum()
        ltools.add_date_column(common_data_in)

        # Applying filters to ratios and store in common_data_filtered
        temp_mean = common_data_in[self.label_ratio].mean()
        temp_stdv = common_data_in[self.label_ratio].std()

        allowed_min_ratio = temp_mean - setts.allowed_ratio_stdv_factor * temp_stdv
        allowed_max_ratio = temp_mean + setts.allowed_ratio_stdv_factor * temp_stdv

        if np.isnan(allowed_min_ratio) or np.isnan(allowed_max_ratio):
            print("WARNING ----- Mean and Std Deviation in " + self.label_ratio + " are undefined")
            allowed_min_ratio = 0
            allowed_max_ratio = 2

        if allowed_min_ratio > setts.ratio_min or np.isnan(allowed_min_ratio):
            allowed_min_ratio = setts.ratio_min
        if allowed_max_ratio < setts.ratio_max or np.isnan(allowed_max_ratio):
            allowed_max_ratio = setts.ratio_max

        common_data_filtered_in = common_data_in[(common_data_in[self.label_ratio] >= allowed_min_ratio) &
                                                 (common_data_in[self.label_ratio] <= allowed_max_ratio)]

        # Resetting indexes
        common_data_filtered_in = common_data_filtered_in.reset_index(drop=True)
        self.__common_data_filtered = common_data_filtered_in
        self.__common_data = common_data_in

        self.__data_exclusion_percent = {"After null and negative luminosity exclusion (excluded data %)": (1.0 - len(
            self.common_data) / len(common_data_in)) * 100}

        self.__fills = np.unique(self.__common_data_filtered['fill'])
        self.__runs = np.unique(self.__common_data_filtered['run'])

        # compute and fill nls data in common_data_filtered
        if fill_nls_data:
            self.fill_nls_data()

        self.__common_data_filtered_no_nan = self.__common_data_filtered.dropna()

        self.__data_all = None
        self.__data_excluded = None

        if self.all_data_analysis_included:
            print ("Merging " + det1.name + " and " + det2.name + " (all data) ... \n")
            common_data_in_all = pd.merge(det1.all_data, det2.all_data, on=keys_for_merging, how='outer')

            common_data_in_all[det1.excluded_label] = common_data_in_all[det1.excluded_label].fillna(det1.label_for_excluded)
            common_data_in_all[det2.excluded_label] = common_data_in_all[det2.excluded_label].fillna(det2.label_for_excluded)

            common_data_in_all[self.label_ratio] = common_data_in_all[det1.lumi_rec_label] / common_data_in_all[det2.lumi_rec_label]

            common_data_in_all[self.__ratio_excluded_label] = np.where(
                (common_data_in_all[det1.excluded_label] == det1.label_for_excluded) | (common_data_in_all[det2.excluded_label] == det2.label_for_excluded),
                det2.label_for_excluded, det2.label_for_included)

            common_data_in_all[self.accumulated_rec_lumi2_label] = common_data_in_all[det2.lumi_rec_label].cumsum()
            common_data_in_all[self.accumulated_rec_lumi1_label] = common_data_in_all[det1.lumi_rec_label].cumsum()

            self.__data_all = common_data_in_all


        # -> Initializing stats vars
        # single ls ratios
        self.__ratios_mean = None
        self.__ratios_stdv = None
        self.__ratios_unfiltered_mean = None
        self.__ratios_unfiltered_stdv = None
        # lumi weighted
        self.__ratios_lw_mean = None
        self.__ratios_lw_stdv = None
        self.__ratios_lw_stdv_dof_corr = None
        # by nls ratios
        self.__nls_ratios_mean = None
        self.__nls_ratios_stdv = None
        # lumi weighted
        self.__nls_ratios_lw_mean = None
        self.__nls_ratios_lw_stdv = None
        self.__nls_ratios_lw_stdv_dof_corr = None

        if fill_stats:
            self.fill_stats()

        # Initializing other variables:
        self.__plt_plots = {}
        self.__sns_plots = {}
        self.__nbx_data = None

    # Other functions
    # TODO: exclude_outlayers function to exclude wrong points from data
    def exclude_outlayers(self):
        return "do something here :)"

    def fill_nbx(self, constant_period=True):
        if constant_period:
            nbx_file_name = 'NBX_files/NBX_' + str(self.year) + '_' + str(self.energy) + 'tev.csv'
            try:
                nbx_by_fill_data = pd.read_csv(nbx_file_name, index_col=False, names=('fill', self.det2.nbx_label))
            except IOError as io_err:
                print("error message : ", io_err)
                print("Problem reading " + nbx_file_name)
                raise
            print("NBX info extracted from " + nbx_file_name)
            nbx_by_fill_data = nbx_by_fill_data.sort_values(by="fill", ascending=True)
            self.nbx_data = nbx_by_fill_data

            self.__common_data_filtered = pd.merge(self.common_data_filtered, nbx_by_fill_data,
                                                   on='fill', how='left', sort=True)

    def fill_nls_data(self):
        # compute by Nls ratios
        nls = self.__nls
        inls = 0
        sum_lumi1 = 0.0
        sum_lumi2 = 0.0
        nls_ratios_array = []

        all_range_nls_array = []
        all_range_nls_err_array = []
        all_range_lumi_in_nls_array = []

        data_len = len(self.common_data_filtered)
        new_fill_in_next_ls = False

        for index_data in range(0, len(self.common_data_filtered)):
            inls += 1
            sum_lumi1 += self.common_data_filtered[self.__det1.lumi_rec_label][index_data]
            sum_lumi2 += self.common_data_filtered[self.__det2.lumi_rec_label][index_data]
            nls_ratios_array.append(self.common_data_filtered[self.label_ratio][index_data])

            all_range_nls = np.nan
            all_range_lumi_in_nls = np.nan
            all_range_nls_err = np.nan

            if index_data != data_len - 1:
                if self.common_data_filtered['fill'][index_data] != self.common_data_filtered['fill'][index_data + 1]:
                    new_fill_in_next_ls = True
                else:
                    new_fill_in_next_ls = False

            if inls == nls or index_data == data_len - 1 or new_fill_in_next_ls:
                by_nls_ratio = sum_lumi1 / sum_lumi2
                by_nls_ratios_stdv = np.array(nls_ratios_array).std()
                by_nls_ratios_mean_err = by_nls_ratios_stdv / np.sqrt(inls)

                all_range_nls = by_nls_ratio
                all_range_lumi_in_nls = sum_lumi2
                all_range_nls_err = by_nls_ratios_mean_err

                inls = 0
                sum_lumi1 = 0.0
                sum_lumi2 = 0.0
                nls_ratios_array = []

            all_range_nls_array.append(all_range_nls)
            all_range_lumi_in_nls_array.append(all_range_lumi_in_nls)
            all_range_nls_err_array.append(all_range_nls_err)

        self.common_data_filtered[self.by_nls_label_ratio] = np.array(all_range_nls_array)
        self.common_data_filtered[self.__by_nls_label_ratio_err] = np.array(all_range_nls_err_array)
        self.common_data_filtered[self.__by_nls_lumi_label] = np.array(all_range_lumi_in_nls_array)

    def fill_stats(self):
        data = self.common_data_filtered

        # non_weighted_stats
        self.__ratios_unfiltered_mean = self.__common_data[self.label_ratio].mean()
        self.__ratios_mean = self.__common_data_filtered[self.label_ratio].mean()
        self.__ratios_unfiltered_stdv = self.__common_data[self.label_ratio].std()
        self.__ratios_stdv = self.__common_data_filtered[self.label_ratio].std()

        self.__nls_ratios_mean = self.__common_data_filtered[self.by_nls_label_ratio].mean()
        self.__nls_ratios_stdv = self.__common_data_filtered[self.by_nls_label_ratio].std()

        lw_stats = ltools.get_w_stats(data[self.label_ratio], data[self.det2.lumi_rec_label])
        nls_lw_stats = ltools.get_w_stats(data[self.by_nls_label_ratio], data[self.by_nls_lumi_label])

        self.__ratios_lw_mean = lw_stats.mean
        self.__ratios_lw_stdv = lw_stats.std_mean
        self.__ratios_lw_stdv_dof_corr = lw_stats.std

        self.__nls_ratios_lw_mean = nls_lw_stats.mean
        self.__nls_ratios_lw_stdv = nls_lw_stats.std_mean
        self.__nls_ratios_lw_stdv_dof_corr = nls_lw_stats.std

    # TODO: implement def fill_equal_lumi_data(self):

    def plot_ratio_hist(self):
        ratio_hist = plotting.hist_from_pandas_frame(data_frame=self.__common_data_filtered, col_label=self.label_ratio,
                                                     nbins=setts.nbins,
                                                     xlabel=self.__label_ratio + " ratio", ylabel='Counts',
                                                     # title='Detectors Ratios Histogram',
                                                     xmin=setts.ratio_min, xmax=setts.ratio_max,
                                                     mean=self.__ratios_mean, stdv=self.__ratios_stdv,
                                                     energy_year_label=self.__year_energy_label)
        self.__plt_plots['ratio_hist'] = ratio_hist[0][0].get_figure()

    def plot_ratio_hist_weighted(self):
        ratio_hist_lumi2_w = plotting.hist_from_pandas_frame(data_frame=self.__common_data_filtered,
                                                             col_label=self.label_ratio,
                                                             nbins=setts.nbins,
                                                             xlabel=self.__label_ratio + " ratio",
                                                             ylabel="Integrated luminosity [$" +
                                                                    self.lumi_unit + "^{-1}$]",
                                                             # title='Detectors Ratios Histogram (lumi weighted)',
                                                             xmin=setts.ratio_min, xmax=setts.ratio_max,
                                                             mean=self.__ratios_lw_mean,
                                                             stdv=self.__ratios_lw_stdv_dof_corr,
                                                             energy_year_label=self.__year_energy_label,
                                                             weight_label=self.det2.lumi_rec_label)
        self.__plt_plots['ratio_hist_lw'] = ratio_hist_lumi2_w[0][0].get_figure()

    def plot_nls_ratio_hist(self):
        ratio_hist = plotting.hist_from_pandas_frame(data_frame=self.common_data_filtered,
                                                     col_label=self.by_nls_label_ratio,
                                                     nbins=setts.nbins,
                                                     xlabel=self.__label_ratio + " ratios in " + str(
                                                         self.__nls) + ' LS',
                                                     ylabel='Counts',
                                                     # title='by Nls Detectors Ratios Histogram',
                                                     xmin=setts.ratio_min, xmax=setts.ratio_max,
                                                     mean=self.__nls_ratios_mean, stdv=self.__nls_ratios_stdv,
                                                     energy_year_label=self.__year_energy_label)
        self.__plt_plots['ratio_nls_hist'] = ratio_hist[0][0].get_figure()

    def plot_nls_ratio_hist_weighted(self):
        ratio_hist_lumi2_w = plotting.hist_from_pandas_frame(data_frame=self.__common_data_filtered_no_nan,
                                                             col_label=self.by_nls_label_ratio,
                                                             nbins=setts.nbins,
                                                             xlabel=self.__label_ratio + " ratios in " + str(
                                                                 self.__nls) + ' LS',
                                                             ylabel="Integrated luminosity [$" +
                                                                    self.lumi_unit + "^{-1}$]",
                                                             # title='Detectors Ratios Histogram (lumi weighted)',
                                                             xmin=setts.ratio_min, xmax=setts.ratio_max,
                                                             mean=self.__nls_ratios_lw_mean,
                                                             stdv=self.__nls_ratios_lw_stdv_dof_corr,
                                                             energy_year_label=self.__year_energy_label,
                                                             weight_label=self.__by_nls_lumi_label)
        self.__plt_plots['nls_ratio_hist_lw'] = ratio_hist_lumi2_w[0][0].get_figure()

    # TODO: stability plots (double axis config)
    # TODO: binning, nls tests plots

    def plot_ratio_vs_time(self):
        ratio_vs_time = plotting.scatter_plot_from_pandas_frame(data_frame=self.__common_data_filtered,
                                                                y_data_label=self.label_ratio, x_data_label='time',
                                                                xlabel='timestamp[s]',
                                                                ylabel=self.__label_ratio + " ratios",
                                                                ymin=setts.ratio_min, ymax=setts.ratio_max,
                                                                energy_year_label=self.__year_energy_label)
        self.__plt_plots['ratio_vs_time'] = ratio_vs_time.get_figure()

    def plot_ratio_vs_date(self):
        ratio_vs_date = plotting.scatter_plot_from_pandas_frame(data_frame=self.__common_data_filtered,
                                                                y_data_label=self.label_ratio, x_data_label='date',
                                                                # xlabel='date',
                                                                ylabel=self.__label_ratio + " ratios",
                                                                ymin=setts.ratio_min, ymax=setts.ratio_max,
                                                                energy_year_label=self.__year_energy_label)
        self.__plt_plots['ratio_vs_date'] = ratio_vs_date.get_figure()

    def plot_ratio_vs_lumi2(self):
        ratio_vs_time = plotting.scatter_plot_from_pandas_frame(data_frame=self.__common_data_filtered,
                                                                y_data_label=self.label_ratio,
                                                                x_data_label=self.accumulated_rec_lumi2_label,
                                                                xlabel="Integrated luminosity [$" +
                                                                       self.lumi_unit + "^{-1}$]",
                                                                ylabel=self.__label_ratio + " ratios",
                                                                ymin=setts.ratio_min, ymax=setts.ratio_max,
                                                                energy_year_label=self.__year_energy_label)
        self.__plt_plots['ratio_vs_lumi2'] = ratio_vs_time.get_figure()

    def plot_nls_ratio_vs_time(self):
        ratio_vs_time = plotting.scatter_plot_from_pandas_frame(data_frame=self.common_data_filtered,
                                                                y_data_label=self.by_nls_label_ratio,
                                                                x_data_label='time',
                                                                xlabel='timestamp[s]',
                                                                ylabel=self.__label_ratio + " ratios in " +
                                                                       str(self.__nls) + ' LS',
                                                                ymin=setts.ratio_min, ymax=setts.ratio_max,
                                                                energy_year_label=self.__year_energy_label)
        self.__plt_plots['by_nls_ratio_vs_time'] = ratio_vs_time.get_figure()

    def plot_nls_ratio_vs_date(self):
        ratio_vs_date = plotting.scatter_plot_from_pandas_frame(data_frame=self.common_data_filtered,
                                                                y_data_label=self.by_nls_label_ratio,
                                                                x_data_label='date',
                                                                ylabel=self.__label_ratio + " ratios in " +
                                                                       str(self.__nls) + ' LS',
                                                                ymin=setts.ratio_min, ymax=setts.ratio_max,
                                                                energy_year_label=self.__year_energy_label)
        self.__plt_plots['by_nls_ratio_vs_date'] = ratio_vs_date.get_figure()

    def plot_nls_ratio_vs_lumi2(self):
        ratio_vs_lumi = plotting.scatter_plot_from_pandas_frame(data_frame=self.common_data_filtered,
                                                                y_data_label=self.by_nls_label_ratio,
                                                                x_data_label=self.accumulated_rec_lumi2_label,
                                                                xlabel="Integrated luminosity [$" +
                                                                       self.lumi_unit + "^{-1}$]",
                                                                ylabel=self.__label_ratio + " ratios in " +
                                                                       str(self.__nls) + ' LS',
                                                                ymin=setts.ratio_min, ymax=setts.ratio_max,
                                                                energy_year_label=self.__year_energy_label)
        self.__plt_plots['by_nls_ratio_vs_lumi2'] = ratio_vs_lumi.get_figure()

    # Things I added

    def plot_nls_ratio_vs_run(self):
        ratio_vs_run = plotting.scatter_plot_from_pandas_frame(data_frame=self.common_data_filtered,
                                                               y_data_label=self.by_nls_label_ratio,
                                                               x_data_label='run',
                                                               xlabel="Run Number",
                                                               ylabel=self.__label_ratio + " ratios in " +
                                                                      str(self.__nls) + ' LS',
                                                               ymin=setts.ratio_min, ymax=setts.ratio_max,
                                                               energy_year_label=self.__year_energy_label)
        self.__plt_plots['by_nls_ratio_vs_run'] = ratio_vs_run.get_figure()

    def plot_nls_ratio_vs_fill(self):
        ratio_vs_fill = plotting.scatter_plot_from_pandas_frame(data_frame=self.common_data_filtered,
                                                                y_data_label=self.by_nls_label_ratio,
                                                                x_data_label='fill',
                                                                xlabel="Fill Number",
                                                                ylabel=self.__label_ratio + " ratios in " +
                                                                       str(self.__nls) + ' LS',
                                                                ymin=setts.ratio_min, ymax=setts.ratio_max,
                                                                energy_year_label=self.__year_energy_label)
        self.__plt_plots['by_nls_ratio_vs_fill'] = ratio_vs_fill.get_figure()

    def plot_ratio_vs_run(self):
        ratio_vs_run = plotting.scatter_plot_from_pandas_frame(data_frame=self.common_data_filtered,
                                                               y_data_label=self.label_ratio,
                                                               x_data_label='run',
                                                               xlabel="Run Number",
                                                               ylabel=self.__label_ratio + " ratios",
                                                               ymin=setts.ratio_min, ymax=setts.ratio_max,
                                                               energy_year_label=self.__year_energy_label)
        self.__plt_plots['by_ratio_vs_run'] = ratio_vs_run.get_figure()

    def plot_ratio_vs_fill(self):
        ratio_vs_fill = plotting.scatter_plot_from_pandas_frame(data_frame=self.common_data_filtered,
                                                                y_data_label=self.label_ratio,
                                                                x_data_label='fill',
                                                                xlabel="Fill Number",
                                                                ylabel=self.__label_ratio + " ratios",
                                                                ymin=setts.ratio_min, ymax=setts.ratio_max,
                                                                energy_year_label=self.__year_energy_label)
        self.__plt_plots['by_ratio_vs_fill'] = ratio_vs_fill.get_figure()

    # End of things I added

    # Fill test
    def plot_bad_fills(self):
        ratio_vs_fill = plotting.plot_bad_fill_info(data_frame=self.common_data_filtered_no_nan, x_data_label='fill',
                                                    y_data_label=self.label_ratio, z_data_label=self.__rec_label,
                                                    xlabel= 'Fill Number', ylabel=self.__label_ratio + " ratios",
                                                    mean=self.__nls_ratios_lw_mean,
                                                    stdv=self.__nls_ratios_lw_stdv_dof_corr,
                                                    ratio_acceptance=setts.bad_ratio_to_plot_stdv_factor,
                                                    filePath=self.output_dir + 'txt/',
                                                    ymin=self.__nls_ratios_lw_mean-3*setts.bad_ratio_to_plot_stdv_factor
                                                         * self.__nls_ratios_lw_stdv_dof_corr,
                                                    ymax=self.__nls_ratios_lw_mean+3*setts.bad_ratio_to_plot_stdv_factor
                                                         * self.__nls_ratios_lw_stdv_dof_corr,
                                                    energy_year_label=self.__year_energy_label,
                                                    txtfileName='Bad_fills')
        self.__plt_plots['by_ratio_vs_fill_bad_fills'] = ratio_vs_fill

    def plot_bad_runs(self):
        ratio_vs_run = plotting.plot_bad_fill_info(data_frame=self.common_data_filtered_no_nan, x_data_label='run',
                                                    y_data_label=self.label_ratio, z_data_label=self.__rec_label,
                                                    xlabel= 'Run Number', ylabel=self.__label_ratio + " ratios",
                                                    mean=self.__nls_ratios_lw_mean,
                                                    stdv=self.__nls_ratios_lw_stdv_dof_corr,
                                                    ratio_acceptance=setts.bad_ratio_to_plot_stdv_factor,
                                                    filePath=self.output_dir + 'txt/',
                                                    ymin=self.__nls_ratios_lw_mean-3*setts.bad_ratio_to_plot_stdv_factor
                                                         *self.__nls_ratios_lw_stdv_dof_corr,
                                                    ymax=self.__nls_ratios_lw_mean+3*setts.bad_ratio_to_plot_stdv_factor
                                                         *self.__nls_ratios_lw_stdv_dof_corr,
                                                    energy_year_label=self.__year_energy_label,
                                                    txtfileName='Bad_runs')
        self.__plt_plots['by_ratio_vs_run_bad_runs'] = ratio_vs_run

    # All/excluded data analysis plots
    def __plot_all_and_excluded(self, det_index):
        if type(det_index)==int and det_index<=1:
            det = self.__detcs[det_index]
            suffix_plot_name = det.excluded_label.replace(" ","_").replace(":","")
            scatter_all_and_excluded_lumi2 = plotting.snsplot_detector_all_and_excluded(self.__data_all,
                                                                                  x_data_label=self.accumulated_rec_lumi2_label,
                                                                                  y_data_label=self.label_ratio,
                                                                                  conditional_label=det.excluded_label,
                                                                                  xlabel="Integrated luminosity [$" +
                                                                                         self.lumi_unit + "^{-1}$]",
                                                                                  ylabel=self.__label_ratio + " ratios",
                                                                                  ymin=setts.ratio_min, ymax=setts.ratio_max,
                                                                                  energy_year_label=self.__year_energy_label)
            self.__sns_plots['scatter_all_and_excluded_vs_lum2_' + suffix_plot_name] = scatter_all_and_excluded_lumi2

            scatter_all_and_excluded_run = plotting.snsplot_detector_all_and_excluded(self.__data_all,
                                                                                  x_data_label="run",
                                                                                  y_data_label=self.label_ratio,
                                                                                  conditional_label=det.excluded_label,
                                                                                  xlabel="Run",
                                                                                  ylabel=self.__label_ratio + " ratios",
                                                                                  ymin=setts.ratio_min,
                                                                                  ymax=setts.ratio_max,
                                                                                  energy_year_label=self.__year_energy_label)
            self.__sns_plots['scatter_all_and_excluded_vs_run_' + suffix_plot_name] = scatter_all_and_excluded_run
        else:
            raise AssertionError("Index must be 0 or 1")

    def plot_all_and_excluded_by_detc(self):
        self.__plot_all_and_excluded(0)
        self.__plot_all_and_excluded(1)

    def save_plots(self):
        print('\n\n Saving plots:')
        plotting.save_plots(self.plt_plots, self.output_dir)
        plotting.save_plots(self.sns_plots, self.output_dir)

    @property
    def output_dir(self):
        return self.__output_dir

    @property
    def label_ratio(self):
        return self.__label_ratio

    @property
    def all_data_analysis_included(self):
        return self.__all_data_analysis_included

    @property
    def by_nls_label_ratio(self):
        return self.__by_nls_label_ratio

    @property
    def by_nls_label_ratio_err(self):
        return self.__by_nls_label_ratio_err

    @property
    def by_nls_lumi1(self):
        return self.__by_nls_lumi_rec_label1

    @property
    def by_nls_lumi2(self):
        return self.__by_nls_lumi_rec_label2

    @property
    def by_nls_lumi_label(self):
        return self.__by_nls_lumi_label

    @property
    def det1(self):
        return self.__det1

    @property
    def det2(self):
        return self.__det2

    @property
    def by_nls_accumulated_rec_lumi1_label(self):
        return self.__by_nls_accumulated_rec_lumi1_label

    @property
    def by_nls_accumulated_rec_lumi2_label(self):
        return self.__by_nls_accumulated_rec_lumi2_label

    @property
    def accumulated_rec_lumi1_label(self):
        return self.__accumulated_rec_lumi1_label

    @property
    def accumulated_rec_lumi2_label(self):
        return self.__accumulated_rec_lumi2_label

    @property
    def year_energy_label(self):
        return self.__year_energy_label

    @property
    def ratio_excluded_label(self):
        return self.__ratio_excluded_label

    @property
    def common_data(self):
        return self.__common_data

    @property
    def common_data_filtered(self):
        return self.__common_data_filtered

    @property
    def common_data_filtered_no_nan(self):
        return self.__common_data_filtered_no_nan

    @property
    def nls(self):
        return self.__nls

    @property
    def nls_data(self):
        if self.__nls_data is None:
            raise AssertionError("please run first the create_nls_data() method")
        return self.__nls_data

    @property
    def data_exclusion_percent(self):
        return self.__data_exclusion_percent

    @property
    def nbx_data(self):
        return self.__nbx_data

    @property
    def ratios_mean(self):
        return self.__ratios_mean

    @property
    def nls_ratios_mean(self):
        if self.__nls_ratios_mean is None:
            raise AssertionError("please run first the create_nls_data() method")
        return self.__nls_ratios_mean

    @property
    def ratios_stdv(self):
        return self.__ratios_stdv

    @property
    def fills(self):
        return self.__fills

    @property
    def runs(self):
        return self.__runs

    @property
    def nls_ratios_stdv(self):
        if self.__nls_ratios_stdv is None:
            raise AssertionError("please run first the create_nls_data() method")
        return self.__nls_ratios_stdv

    @property
    def plt_plots(self):
        return self.__plt_plots

    @property
    def sns_plots(self):
        return self.__sns_plots

    # setters
    @common_data.setter
    def common_data(self, data):
        self.__common_data = data

    @nls_data.setter
    def nls_data(self, data):
        self.__nls_data = data

    @nbx_data.setter
    def nbx_data(self, nbx):
        self.__nbx_data = nbx

    @output_dir.setter
    def output_dir(self, dir_name):
        self.__output_dir = dir_name


# TODO: make numb of detectors automatic
class MultipleDetectorsRatio:
    def __init__(self, det1: L, det2: L, det3: L, year: str = None, energy: str = None, c_years = False) -> None:

        self.__dets = [det1, det2, det3]
        number_of_dets = len(self.__dets)
        self.__lumi_unit = det3.lumi_unit

        if year is None:
            if (det1.year == det2.year) and (det1.energy == det2.energy) \
                    and (det1.year == det3.year) and (det1.energy == det3.energy):
                self.year = det1.year
            else:
                raise ValueError("Comparing detectors from different years or energies makes no sense :(")
        else:
            if c_years:
                list_years = year.split(',')
                self.year = list_years[0]
            else:
                self.year = year

        if energy is None:
            if (det1.year == det2.year) and (det1.energy == det2.energy) \
                    and (det1.year == det3.year) and (det1.energy == det3.energy):
                self.energy = det1.energy
            else:
                raise ValueError("Comparing detectors from different years or energies makes no sense :(")
        else:
            self.energy = energy

        assert det1.lumi_unit == det3.lumi_unit and det2.lumi_unit == det3.lumi_unit
        self.__lumi_unit = det3.lumi_unit

        if c_years:
            self.__output_dir = setts.default_output_dir + year.replace(',', '-') + '/' + det1.name + '-' + det2.name \
                                + '-' + det3.name + '/'
            self.__year_energy_label = str(self.energy) + 'TeV(' + year + ')'
        else:
            self.__output_dir = setts.default_output_dir + str(self.year) + '/' + det1.name + '-' + det2.name \
                                + '-' + det3.name + '/'
            self.__year_energy_label = str(self.energy) + 'TeV(' + str(self.year) + ')'

        #Checking all_data_analysis_included for all detectors
        self.__all_data_analysis_included = True
        for det in self.__dets:
            if not det.all_data_analysis_included:
                self.__all_data_analysis_included = False

        self.__ratios = []
        self.__ratios_all = []

        # TODO: implement full_data comparison
        for i in range(0, number_of_dets):
            for j in range(i + 1, number_of_dets):
                print('Setting ' + str(self.__dets[i].name) + '/' + str(self.__dets[j].name) + ' ...')
                self.__ratios.append(DetectorsRatio(self.__dets[i], self.__dets[j]))

        self.__ref_ratio = self.__ratios[len(self.__ratios) - 1]
        self.__nls = self.__ref_ratio.nls

        self.__nls_ratio_col_names = []
        self.__ratio_plotting_labels = []
        for i_ratio in self.__ratios:
            i_ratio.fill_nls_data()
            self.__nls_ratio_col_names.append(i_ratio.by_nls_label_ratio)
            self.__ratio_plotting_labels.append(i_ratio.label_ratio)
        self.__lumi3_col_name = self.__ref_ratio.accumulated_rec_lumi2_label
        self.__ratio_col_names = self.__ratio_plotting_labels

        keys_for_merging = ['fill', 'run', 'ls', 'time', 'date']

        # TODO: Make this automatic for the number of detectors
        if len(self.__dets) == 3:
            merge_tmp = pd.merge(self.__ratios[0].common_data_filtered, self.__ratios[1].common_data_filtered,
                                 on=keys_for_merging,
                                 how='outer').merge(self.__ratios[2].common_data_filtered,
                                                    on=keys_for_merging, how='outer')
            # if full_data:
            #     merge_tmp_all = pd.merge(self.__ratios_all[0].common_data_filtered,
            #                              self.__ratios_all[1].common_data_filtered,
            #                              on=keys_for_merging,
            #                              how='outer').merge(self.__ratios_all[2].common_data_filtered,
            #                                                 on=keys_for_merging, how='outer')
            #     merge_tmp_all = merge_tmp_all.reset_index(drop=True)
            #     ltools.check_and_clean_after_merging(merge_tmp_all)
            #
            #     merge_tmp_excluded = pd.concat([merge_tmp_all, merge_tmp, merge_tmp], sort=True).drop_duplicates(keep=False)
            #     merge_tmp_excluded = merge_tmp_excluded.reset_index(drop=True)
            #     ltools.check_and_clean_after_merging(merge_tmp_excluded)
            #
            #     self.__combined_data_all = merge_tmp_all
            #     self.__combined_data_excluded = merge_tmp_excluded
            #
            #     print(len(merge_tmp_excluded), len(merge_tmp), len(merge_tmp_all))
            #     print (" \n Data percent excluded in normtags (%): " + str((len(merge_tmp_all) - len(merge_tmp))*100/len(merge_tmp_all)) + "\n")

        else:
            raise ValueError('Number of detectors not implemented yet :(')

        merge_tmp = merge_tmp.reset_index(drop=True)
        ltools.check_and_clean_after_merging(merge_tmp)

        self.__combined_data = merge_tmp
        #print(list(merge_tmp))

        self.__plt_plots = {}
        self.__sns_plots = []

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

    # Things I added

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

    # End of things I added

    # TODO: create combined hists plot

    @property
    def combined_data(self):
        return self.__combined_data

    def save_plots(self):
        print('\n\n Saving plots:')
        plotting.save_plots(self.__plt_plots, self.__output_dir)
