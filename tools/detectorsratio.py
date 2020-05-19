from tools.luminometer import Luminometer as L
import numpy as np
import pandas as pd
import settings as setts
from tools import lumi_tools as ltools
import tools.plotting_tools as plotting
import warnings


class DetectorsRatio(L):

    def __init__(self, det1: L, det2: L, year: str = None, energy: str = None,
                 fill_nls_data: bool = True, fill_stats=True, c_years=False,
                 nls=None, only_stats: bool = False, lumi_type: str = 'rec',
                 compute_by_run_by_fill: bool = False, fill_norm_ratios: bool = False) -> None:

        self.__detcs = (det1, det2)
        self.__only_stats = only_stats
        self.__compute_by_run_by_fill = compute_by_run_by_fill
        self.__stats_DF = None
        self.__normalized_data_filled = False

        if lumi_type == 'rec':
            self.__label_det1_lumi_to_use = det1.lumi_rec_label
            self.__label_det2_lumi_to_use = det2.lumi_rec_label
        elif lumi_type == 'del':
            self.__label_det1_lumi_to_use = det1.lumi_del_label
            self.__label_det2_lumi_to_use = det2.lumi_del_label
        else:
            raise AssertionError("lumitype not implemented!: " + lumi_type)

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
            self.output_dir = setts.default_output_dir + year.replace(',',
                                                                      '-') + '/' + det1.name + '-' + det2.name + "_" \
                              + str(self.energy) + 'TeV/'
            if (year == "2016,2017,2018"):
                year_label = "full RunII"
            else:
                year_label = year
            self.__year_energy_label = str(self.energy) + 'TeV(' + year_label + ')'
        else:
            self.output_dir = setts.default_output_dir + str(self.year) + "_" + str(self.energy) + 'TeV/' + det1.name \
                              + '-' + det2.name + '/'
            self.__year_energy_label = str(self.energy) + 'TeV(' + str(self.year) + ')'

        self.__label_ratio = det1.name + '/' + det2.name
        self.__label_ratio_norm = det1.name + '/' + det2.name + "_norm"
        self.__by_nls_label_ratio = 'by_' + str(self.__nls) + 'nls_' + self.__label_ratio
        self.__by_nls_label_ratio_norm = 'by_' + str(self.__nls) + 'nls_' + self.__label_ratio + "_norm"
        self.__lumi_unit = det2.lumi_unit
        self.__det1 = det1
        self.__det2 = det2
        self.__by_nls_lumi_label1 = self.__label_det1_lumi_to_use + '_by_' + str(self.__nls) + 'nls'
        self.__by_nls_lumi_label2 = self.__label_det2_lumi_to_use + '_by_' + str(self.__nls) + 'nls'
        self.__by_nls_label_ratio_err = self.__by_nls_label_ratio + '_err'
        self.__by_nls_lumi_label = 'by_' + str(self.__nls) + 'nls_lumi_' + self.__label_ratio
        self.__by_nls_lumi_label_norm = 'by_' + str(self.__nls) + 'nls_lumi_' + self.__label_ratio + "_norm"

        if compute_by_run_by_fill:
            self.__by_run_lumi_label1 = self.__label_det1_lumi_to_use + '_by_run'
            self.__by_run_lumi_label2 = self.__label_det2_lumi_to_use + '_by_run'
            self.__by_fill_lumi_label1 = self.__label_det1_lumi_to_use + '_by_fill'
            self.__by_fill_lumi_label2 = self.__label_det2_lumi_to_use + '_by_fill'
            self.__by_run_label_ratio = 'by_run_' + self.__label_ratio
            self.__by_fill_label_ratio = 'by_fill_' + self.__label_ratio
            self.__by_run_label_ratio_err = self.__by_run_label_ratio + '_err'
            self.__by_fill_label_ratio_err = self.__by_fill_label_ratio + '_err'
            self.__by_run_lumi_label = 'by_run_lumi_' + self.__label_ratio
            self.__by_fill_lumi_label = 'by_fill_lumi_' + self.__label_ratio

        self.__accumulated_lumi1_label = self.__label_det1_lumi_to_use + '_accumulated'
        self.__accumulated_lumi2_label = self.__label_det2_lumi_to_use + '_accumulated'
        self.__by_nls_accumulated_lumi1_label = self.__by_nls_lumi_label1 + '_accumulated'
        self.__by_nls_accumulated_lumi2_label = self.__by_nls_lumi_label2 + '_accumulated'
        self.__ratio_excluded_label = 'Exclusion info (' + self.__label_ratio + ')'

        # contains only: fully included, partially excluded, partially included
        self.__by_nls_exclusion_info_label = 'Exclusion info (' + self.__label_ratio + ')(' + str(self.__nls) + ' nls)'

        # contains only: included, excluded
        self.__by_nls_binary_exclusion_info_label = 'by nls exclusion (' + self.__label_ratio + ')'

        # contains percent of exclusion
        self.__by_nls_exclusion_percent_label = 'Exclusion percent (' + self.__label_ratio + ')(' + str(
            self.__nls) + ' nls)'

        keys_for_merging = ['ls', 'time', 'run', 'fill']
        print("Merging " + det1.name + " and " + det2.name + " ...")
        common_data_in = pd.merge(det1.data, det2.data, on=keys_for_merging, how='outer')

        # Computing single ls ratio
        common_data_in[self.label_ratio] = common_data_in[self.__label_det1_lumi_to_use] / common_data_in[
            self.__label_det2_lumi_to_use]

        # Adding cumulative luminosity column
        common_data_in[self.accumulated_lumi2_label] = common_data_in[self.__label_det2_lumi_to_use].cumsum()
        common_data_in[self.accumulated_lumi1_label] = common_data_in[self.__label_det1_lumi_to_use].cumsum()
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

        self.__data_all = None
        self.__data_excluded = None

        # Excluded data analysis
        if self.all_data_analysis_included:
            print("Merging " + det1.name + " and " + det2.name + " (all data) ... \n")
            common_data_in_all = pd.merge(det1.all_data, det2.all_data, on=keys_for_merging, how='outer')

            common_data_in_all[det1.excluded_label] = common_data_in_all[det1.excluded_label].fillna(
                det1.label_for_excluded)
            common_data_in_all[det2.excluded_label] = common_data_in_all[det2.excluded_label].fillna(
                det2.label_for_excluded)

            common_data_in_all[self.label_ratio] = common_data_in_all[self.__label_det1_lumi_to_use] / \
                                                   common_data_in_all[self.__label_det2_lumi_to_use]

            ratio_excluded_label_list = []
            for index, row in common_data_in_all.iterrows():
                if row[det1.excluded_label] != det1.label_for_excluded and row[
                    det2.excluded_label] != det2.label_for_excluded:
                    ratio_excluded_label_list.append("included")
                elif row[det1.excluded_label] == det1.label_for_excluded and row[
                    det2.excluded_label] == det2.label_for_excluded:
                    ratio_excluded_label_list.append("both excluded")
                elif row[det1.excluded_label] == det1.label_for_excluded and row[
                    det2.excluded_label] != det2.label_for_excluded:
                    ratio_excluded_label_list.append("only " + det1.name + " excluded")
                elif row[det1.excluded_label] != det1.label_for_excluded and row[
                    det2.excluded_label] == det2.label_for_excluded:
                    ratio_excluded_label_list.append("only " + det2.name + " excluded")
                else:
                    raise ValueError(
                        "None of the required [if] conditions fulfilled! row info: " + row[det1.excluded_label] + ", " +
                        row[det2.excluded_label])

            common_data_in_all[self.__ratio_excluded_label] = np.array(ratio_excluded_label_list)

            common_data_in_all[self.accumulated_lumi2_label] = common_data_in_all[
                self.__label_det2_lumi_to_use].cumsum()
            common_data_in_all[self.accumulated_lumi1_label] = common_data_in_all[
                self.__label_det1_lumi_to_use].cumsum()

            ltools.add_date_column(common_data_in_all)

            self.__data_all = common_data_in_all

        # compute and fill nls data in common_data_filtered
        if fill_nls_data:
            self.fill_nls_data()
            if self.all_data_analysis_included:
                self.fill_nls_all_data()
        if compute_by_run_by_fill:
            self.fill_by_period_data('run')
            self.fill_by_period_data('fill')

        # -> Initializing stats vars
        # single ls ratios
        self.__ratios_mean = None
        self.__ratios_stdv = None
        self.__ratios_unfiltered_mean = None
        self.__ratios_unfiltered_stdv = None
        # lumi weighted
        self.__ratios_lw_mean = None
        self.__ratios_lw_mean_error = None
        self.__ratios_lw_stdv = None
        # by nls ratios
        self.__nls_ratios_mean = None
        self.__nls_ratios_stdv = None
        # lumi weighted
        self.__nls_ratios_lw_mean = None
        self.__nls_ratios_lw_mean_error = None
        self.__nls_ratios_lw_stdv = None
        # by run, by fill
        self.__by_run_ratios_mean = None
        self.__by_run_ratios_stdv = None
        self.__by_fill_ratios_mean = None
        self.__by_fill_ratios_stdv = None
        self.__by_run_lw_mean = None
        self.__by_run_lw_mean_error = None
        self.__by_run_lw_stdv = None
        self.__by_fill_lw_mean = None
        self.__by_fill_lw_mean_error = None
        self.__by_fill_lw_stdv = None

        if fill_stats:
            self.fill_stats()

        if fill_norm_ratios:
            self.fill_normalized_detector_ratios()

        self.__common_data_filtered_no_nan = self.__common_data_filtered.dropna()

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
        print('Computing ' + self.label_ratio + ' by nls ratios ...')

        nls = self.__nls
        inls = 0
        sum_lumi1 = 0.0
        sum_lumi2 = 0.0
        nls_ratios_array = []

        all_range_nls_array = []
        all_range_nls_err_array = []
        all_range_lumi_in_nls_array = []

        data_to_use = self.common_data_filtered

        data_len = len(data_to_use)
        new_fill_in_next_ls = False

        for index_data in range(0, len(data_to_use)):
            inls += 1
            sum_lumi1 += data_to_use[self.__label_det1_lumi_to_use][index_data]
            sum_lumi2 += data_to_use[self.__label_det2_lumi_to_use][index_data]
            nls_ratios_array.append(data_to_use[self.label_ratio][index_data])

            all_range_nls = np.nan
            all_range_lumi_in_nls = np.nan
            all_range_nls_err = np.nan

            if index_data != data_len - 1:
                if data_to_use['fill'][index_data] != data_to_use['fill'][index_data + 1]:
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
                number_of_excluded_points_in_nls = 0
                nls_ratios_array = []

            all_range_nls_array.append(all_range_nls)
            all_range_lumi_in_nls_array.append(all_range_lumi_in_nls)
            all_range_nls_err_array.append(all_range_nls_err)

        data_to_use[self.by_nls_label_ratio] = np.array(all_range_nls_array)
        data_to_use[self.__by_nls_label_ratio_err] = np.array(all_range_nls_err_array)
        data_to_use[self.__by_nls_lumi_label] = np.array(all_range_lumi_in_nls_array)

    def fill_by_period_data(self, period: str):
        # Currently period can be only run or fill
        # compute by x ratios
        print('Computing ' + self.label_ratio + ' by ' + period + ' ratios ...')

        sum_lumi1 = 0.0
        sum_lumi2 = 0.0
        by_period_ratios_array = []

        all_range_by_period_array = []
        all_range_by_period_err_array = []
        all_range_lumi_in_period_array = []

        data_to_use = self.common_data_filtered
        data_len = len(data_to_use)

        new_period_flag = False
        old_period = data_to_use[period][0]
        inls = 0

        all_range_by_period = np.nan
        all_range_lumi_by_period = np.nan
        all_range_by_period_err = np.nan

        for index_data in range(0, len(data_to_use)):

            this_period = data_to_use[period][index_data]
            eof_reached = (index_data == data_len - 1)

            if old_period != this_period:
                new_period_flag = True

            if new_period_flag or eof_reached:
                if eof_reached and not new_period_flag:
                    sum_lumi1 += data_to_use[self.__label_det1_lumi_to_use][index_data]
                    sum_lumi2 += data_to_use[self.__label_det2_lumi_to_use][index_data]
                    by_period_ratios_array.append(data_to_use[self.label_ratio][index_data])
                # print("------->>> " + str(old_period))
                # print("sum_lumi1: " + str(sum_lumi1))
                # print("sum_lumi2: " + str(sum_lumi2))
                # print("ratios: ", by_period_ratios_array)
                by_period_ratio = sum_lumi1 / sum_lumi2
                by_period_ratios_stdv = np.array(by_period_ratios_array).std()
                by_period_ratios_mean_err = by_period_ratios_stdv / np.sqrt(inls)

                all_range_by_period = by_period_ratio
                all_range_lumi_by_period = sum_lumi2
                all_range_by_period_err = by_period_ratios_mean_err

            all_range_by_period_array.append(all_range_by_period)
            all_range_by_period_err_array.append(all_range_by_period_err)
            all_range_lumi_in_period_array.append(all_range_lumi_by_period)

            # resetting new period flag
            if new_period_flag:
                new_period_flag = False
                by_period_ratios_array = []
                inls = 0

            sum_lumi1 += data_to_use[self.__label_det1_lumi_to_use][index_data]
            sum_lumi2 += data_to_use[self.__label_det2_lumi_to_use][index_data]
            by_period_ratios_array.append(data_to_use[self.label_ratio][index_data])
            old_period = this_period
            inls += 1

        if period == "run":
            by_period_label_ratio = self.__by_run_label_ratio
            by_period_label_ratio_err = self.__by_run_label_ratio_err
            by_period_lumi_label = self.__by_run_lumi_label
        elif period == "fill":
            by_period_label_ratio = self.__by_fill_label_ratio
            by_period_label_ratio_err = self.__by_fill_label_ratio_err
            by_period_lumi_label = self.__by_fill_lumi_label
        else:
            raise AssertionError("period option not available!!")

        data_to_use[by_period_label_ratio] = np.array(all_range_by_period_array)
        data_to_use[by_period_label_ratio_err] = np.array(all_range_by_period_err_array)
        data_to_use[by_period_lumi_label] = np.array(all_range_lumi_in_period_array)

    def fill_nls_all_data(self):
        # compute by Nls ratios
        if not self.all_data_analysis_included:
            raise AssertionError("It is impossible running this function without loading _all.csv data")

        print('Computing ' + self.label_ratio + '(all data) by nls ratios ...')

        nls = self.__nls
        inls = 0
        sum_lumi1 = 0.0
        sum_lumi2 = 0.0
        nls_ratios_array = []

        all_range_nls_array = []
        all_range_nls_err_array = []
        all_range_lumi_in_nls_array = []

        data_to_use = self.__data_all

        data_len = len(data_to_use)
        new_fill_in_next_ls = False

        number_of_excluded_points_in_nls = 0
        by_nls_exclusion_info_list = []
        by_nls_exclusion_binary_info_list = []
        by_nls_exclusion_percent_list = []

        for index_data in range(0, len(data_to_use)):
            inls += 1
            sum_lumi1 += data_to_use[self.__label_det1_lumi_to_use][index_data]
            sum_lumi2 += data_to_use[self.__label_det2_lumi_to_use][index_data]
            nls_ratios_array.append(data_to_use[self.label_ratio][index_data])
            if data_to_use[self.__ratio_excluded_label][index_data] != "included":
                number_of_excluded_points_in_nls += 1

            all_range_nls = np.nan
            all_range_lumi_in_nls = np.nan
            all_range_nls_err = np.nan
            exclusion_percent = 0.0
            exclusion_decision_label = "fully included"
            binary_exclusion_decision_label = "included"

            if index_data != data_len - 1:
                if data_to_use['fill'][index_data] != data_to_use['fill'][index_data + 1]:
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

                exclusion_percent = number_of_excluded_points_in_nls * 100. / inls
                if exclusion_percent == 0:
                    exclusion_decision_label = "fully included"
                    binary_exclusion_decision_label = "included"
                elif exclusion_percent > setts.by_nls_exclusion_threshold:
                    exclusion_decision_label = "partially excluded"
                    binary_exclusion_decision_label = "excluded"
                else:
                    exclusion_decision_label = "partially included"
                    binary_exclusion_decision_label = "included"

                inls = 0
                sum_lumi1 = 0.0
                sum_lumi2 = 0.0
                number_of_excluded_points_in_nls = 0
                nls_ratios_array = []

            all_range_nls_array.append(all_range_nls)
            all_range_lumi_in_nls_array.append(all_range_lumi_in_nls)
            all_range_nls_err_array.append(all_range_nls_err)

            by_nls_exclusion_percent_list.append(exclusion_percent)
            by_nls_exclusion_info_list.append(exclusion_decision_label)
            by_nls_exclusion_binary_info_list.append(binary_exclusion_decision_label)

        data_to_use[self.by_nls_label_ratio] = np.array(all_range_nls_array)
        data_to_use[self.__by_nls_label_ratio_err] = np.array(all_range_nls_err_array)
        data_to_use[self.__by_nls_lumi_label] = np.array(all_range_lumi_in_nls_array)
        data_to_use[self.__by_nls_exclusion_info_label] = np.array(by_nls_exclusion_info_list)
        data_to_use[self.__by_nls_binary_exclusion_info_label] = np.array(by_nls_exclusion_binary_info_list)
        data_to_use[self.__by_nls_exclusion_percent_label] = np.array(by_nls_exclusion_percent_list)

    def fill_stats(self):
        print("Computing " + self.label_ratio + " statistics ...")
        data_to_use = self.common_data_filtered
        stats_names = []
        stats_values = []
        # non_weighted_stats
        self.__ratios_unfiltered_mean = data_to_use[self.label_ratio].mean()
        self.__ratios_mean = data_to_use[self.label_ratio].mean()
        self.__ratios_unfiltered_stdv = data_to_use[self.label_ratio].std()
        self.__ratios_stdv = data_to_use[self.label_ratio].std()
        self.__nls_ratios_mean = data_to_use[self.by_nls_label_ratio].mean()
        self.__nls_ratios_stdv = data_to_use[self.by_nls_label_ratio].std()
        lw_stats = ltools.get_w_stats(data_to_use[self.label_ratio], data_to_use[self.__label_det2_lumi_to_use],
                                      min_val=setts.ratio_min, max_val=setts.ratio_max)
        nls_lw_stats = ltools.get_w_stats(data_to_use[self.by_nls_label_ratio], data_to_use[self.by_nls_lumi_label],
                                          min_val=setts.ratio_min, max_val=setts.ratio_max)
        if self.__compute_by_run_by_fill:
            self.__by_run_ratios_mean = data_to_use[self.__by_run_label_ratio].mean()
            self.__by_run_ratios_stdv = data_to_use[self.__by_run_label_ratio].std()
            self.__by_fill_ratios_mean = data_to_use[self.__by_fill_label_ratio].mean()
            self.__by_fill_ratios_stdv = data_to_use[self.__by_fill_label_ratio].std()
            by_run_lw_stats = ltools.get_w_stats(data_to_use[self.__by_run_label_ratio],
                                                 data_to_use[self.__by_run_lumi_label],
                                                 min_val=setts.ratio_min, max_val=setts.ratio_max)
            by_fill_lw_stats = ltools.get_w_stats(data_to_use[self.__by_fill_label_ratio],
                                                  data_to_use[self.__by_fill_lumi_label],
                                                  min_val=setts.ratio_min, max_val=setts.ratio_max)
            self.__by_run_lw_mean = by_run_lw_stats.mean
            self.__by_run_lw_mean_error = by_run_lw_stats.std_mean
            self.__by_run_lw_stdv = by_run_lw_stats.std
            self.__by_fill_lw_mean = by_fill_lw_stats.mean
            self.__by_fill_lw_mean_error = by_fill_lw_stats.std_mean
            self.__by_fill_lw_stdv = by_fill_lw_stats.std

        self.__ratios_lw_mean = lw_stats.mean
        self.__ratios_lw_mean_error = lw_stats.std_mean
        self.__ratios_lw_stdv = lw_stats.std
        self.__nls_ratios_lw_mean = nls_lw_stats.mean
        self.__nls_ratios_lw_mean_error = nls_lw_stats.std_mean
        self.__nls_ratios_lw_stdv = nls_lw_stats.std

        # Filling stats DF arrays:
        stats_names.append("ratios_unfiltered_mean")
        stats_values.append(self.__ratios_unfiltered_mean)
        stats_names.append("ratios_mean")
        stats_values.append(self.__ratios_mean)
        stats_names.append("ratios_unfiltered_stdv")
        stats_values.append(self.__ratios_unfiltered_stdv)
        stats_names.append("ratios_stdv")
        stats_values.append(self.__ratios_stdv)
        stats_names.append("nls_ratios_mean")
        stats_values.append(self.__nls_ratios_mean)
        stats_names.append("nls_ratios_stdv")
        stats_values.append(self.__nls_ratios_stdv)
        stats_names.append("ratios_lw_mean")
        stats_values.append(self.__ratios_lw_mean)
        stats_names.append("ratios_lw_mean_error")
        stats_values.append(self.__ratios_lw_mean_error)
        stats_names.append("ratios_lw_stdv")
        stats_values.append(self.__ratios_lw_stdv)
        stats_names.append("nls_ratios_lw_mean")
        stats_values.append(self.__nls_ratios_lw_mean)
        stats_names.append("nls_ratios_lw_mean_error")
        stats_values.append(self.__nls_ratios_lw_mean_error)
        stats_names.append("nls_ratios_lw_stdv")
        stats_values.append(self.__nls_ratios_lw_stdv)

        self.__stats_DF = pd.DataFrame()
        for col_id in range(0, len(stats_names)):
            self.__stats_DF[stats_names[col_id]] = [stats_values[col_id]]
        self.save_stats_to_file()

    def fill_normalized_detector_ratios(self):
        print("Adding " + self.label_ratio + " normalized data ...")
        # pd.options.mode.chained_assignment = None
        data_to_use = self.__common_data_filtered
        ratio_mean_factor = self.__ratios_mean
        by_nls_ratio_mean_factor = self.__nls_ratios_mean

        ratios = []
        nls_ratios = []

        for index_data in range(0, len(data_to_use)):
            ratios.append(data_to_use[self.__label_ratio][index_data]/ratio_mean_factor)
            nls_ratios.append(data_to_use[self.__by_nls_label_ratio][index_data]/by_nls_ratio_mean_factor)

        data_to_use[self.__label_ratio_norm] = np.array(ratios)
        data_to_use[self.__by_nls_label_ratio_norm] = np.array(nls_ratios)

        # data_to_use[self.__label_ratio_norm] = data_to_use[self.__label_ratio] / ratio_mean_factor
        # data_to_use[self.__by_nls_label_ratio_norm] = data_to_use[self.__by_nls_label_ratio] / by_nls_ratio_mean_factor

        self.__normalized_data_filled = True

    # TODO: implement def fill_equal_lumi_data(self):

    def save_stats_to_file(self):
        # stats_file = open(self.output_dir + "stats.csv")
        ltools.check_and_create_folder(self.output_dir, creation_info=False)
        self.__stats_DF.to_csv(index=False, path_or_buf=self.output_dir + "stats.csv")
        # stats_file.close()

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
                                                             stdv=self.__ratios_lw_stdv,
                                                             energy_year_label=self.__year_energy_label,
                                                             weight_label=self.__label_det2_lumi_to_use)
        self.__plt_plots['ratio_hist_lw'] = ratio_hist_lumi2_w[0][0].get_figure()

    def plot_nls_ratio_hist(self):
        ratio_hist = plotting.hist_from_pandas_frame(data_frame=self.common_data_filtered,
                                                     col_label=self.by_nls_label_ratio,
                                                     nbins=setts.nbins,
                                                     xlabel=self.__label_ratio + " ratio (" + str(
                                                         self.__nls) + ' LS)',
                                                     ylabel='Counts',
                                                     # title='by Nls Detectors Ratios Histogram',
                                                     xmin=setts.ratio_min, xmax=setts.ratio_max,
                                                     mean=self.__nls_ratios_mean, stdv=self.__nls_ratios_stdv,
                                                     # err_mean=self.__nls_ratios_lw_mean_error,
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
                                                             stdv=self.__nls_ratios_lw_stdv,
                                                             # err_mean=self.__nls_ratios_lw_mean_error,
                                                             energy_year_label=self.__year_energy_label,
                                                             weight_label=self.__by_nls_lumi_label)
        self.__plt_plots['nls_ratio_hist_lw'] = ratio_hist_lumi2_w[0][0].get_figure()

    def plot_nls_ratio_hist_weighted_norm(self):
        if self.__normalized_data_filled:
            ratio_hist_lumi2_w = plotting.hist_from_pandas_frame(data_frame=self.__common_data_filtered_no_nan,
                                                                 col_label=self.__by_nls_label_ratio_norm,
                                                                 nbins=setts.nbins,
                                                                 xlabel=self.__label_ratio + " ratios in " + str(
                                                                     self.__nls) + ' LS',
                                                                 ylabel="Integrated luminosity [$" +
                                                                        self.lumi_unit + "^{-1}$]",
                                                                 # title='Detectors Ratios Histogram (lumi weighted)',
                                                                 xmin=setts.ratio_min, xmax=setts.ratio_max,
                                                                 stdv=self.__nls_ratios_lw_stdv,
                                                                 # err_mean=self.__nls_ratios_lw_mean_error,
                                                                 energy_year_label=self.__year_energy_label,
                                                                 weight_label=self.__by_nls_lumi_label)
            self.__plt_plots['normalized_nls_ratio_hist_lw'] = ratio_hist_lumi2_w[0][0].get_figure()
        else:
            warnings.warn("No normalized data has been filled. Plot normalized_nls_ratio_hist_lw won't be produced")

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
                                                                x_data_label=self.accumulated_lumi2_label,
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
                                                                x_data_label=self.accumulated_lumi2_label,
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

    def plot_nls_ratio_vs_lumi2_norm(self):
        if self.__normalized_data_filled:
            ratio_vs_lumi = plotting.scatter_plot_from_pandas_frame(data_frame=self.common_data_filtered_no_nan,
                                                                    y_data_label=self.__by_nls_label_ratio_norm,
                                                                    x_data_label=self.accumulated_lumi2_label,
                                                                    xlabel="Integrated luminosity [$" +
                                                                           self.lumi_unit + "^{-1}$]",
                                                                    ylabel=self.__label_ratio + " ratios in " +
                                                                           str(self.__nls) + ' LS',
                                                                    ymin=setts.ratio_min, ymax=setts.ratio_max,
                                                                    energy_year_label=self.__year_energy_label)
            self.__plt_plots['normalized_by_nls_ratio_vs_lumi2'] = ratio_vs_lumi.get_figure()
        else:
            warnings.warn("No normalized data has been filled. Plot normalized_by_nls_ratio_vs_lumi2 won't be produced")

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

    def plot_by_fill_ratio_vs_fill(self):
        by_fill_ratio_vs_fill = plotting.scatter_plot_from_pandas_frame(data_frame=self.common_data_filtered,
                                                                        y_data_label=self.__by_fill_label_ratio,
                                                                        x_data_label='fill',
                                                                        xlabel="Fill Number",
                                                                        ylabel=self.__label_ratio + " ratios per Fill",
                                                                        ymin=setts.ratio_min_reduced,
                                                                        ymax=setts.ratio_max_reduced,
                                                                        energy_year_label=self.__year_energy_label,
                                                                        marker_size=3.0)
        self.__plt_plots['by_fill_ratio_vs_fill'] = by_fill_ratio_vs_fill.get_figure()

    def plot_by_run_ratio_vs_run(self):
        by_run_ratio_vs_run = plotting.scatter_plot_from_pandas_frame(data_frame=self.common_data_filtered,
                                                                      y_data_label=self.__by_run_label_ratio,
                                                                      x_data_label='run',
                                                                      xlabel="Run Number",
                                                                      ylabel=self.__label_ratio + " ratios per Run",
                                                                      ymin=setts.ratio_min_reduced,
                                                                      ymax=setts.ratio_max_reduced,
                                                                      energy_year_label=self.__year_energy_label,
                                                                      marker_size=3.0)
        self.__plt_plots['by_run_ratio_vs_run'] = by_run_ratio_vs_run.get_figure()

    def plot_by_run_ratio_vs_run_with_errors(self):
        by_run_ratio_vs_run_with_errors = plotting.plot_scatter_and_errors(data_frame=self.common_data_filtered,
                                                                           y_data_label=self.__by_run_label_ratio,
                                                                           x_data_label='run',
                                                                           err_y_data_label=self.__by_run_label_ratio_err,
                                                                           xlabel="Run Number",
                                                                           ylabel=self.__label_ratio + " ratios per Run",
                                                                           ymin=setts.ratio_min_reduced,
                                                                           ymax=setts.ratio_max_reduced,
                                                                           energy_year_label=self.__year_energy_label)
        self.__plt_plots['by_run_ratio_vs_run_with_errors'] = by_run_ratio_vs_run_with_errors

    def plot_by_run_hist_weighted(self):
        ratio_hist_lumi2_w = plotting.hist_from_pandas_frame(data_frame=self.__common_data_filtered_no_nan,
                                                             col_label=self.__by_run_label_ratio,
                                                             nbins=setts.nbins_by_run,
                                                             xlabel=self.__label_ratio + "by run ratio",
                                                             ylabel="Integrated luminosity [$" +
                                                                    self.lumi_unit + "^{-1}$]",
                                                             # title='Detectors Ratios Histogram (lumi weighted)',
                                                             xmin=setts.ratio_min_reduced,
                                                             xmax=setts.ratio_max_reduced,
                                                             mean=self.__by_run_lw_mean,
                                                             stdv=self.__by_run_lw_stdv,
                                                             energy_year_label=self.__year_energy_label,
                                                             weight_label=self.__by_run_lumi_label)
        self.__plt_plots['by_run_hist_lw'] = ratio_hist_lumi2_w[0][0].get_figure()

    def plot_by_run_hist(self):
        ratio_hist_lumi2_w = plotting.hist_from_pandas_frame(data_frame=self.__common_data_filtered_no_nan,
                                                             col_label=self.__by_run_label_ratio,
                                                             nbins=setts.nbins_by_run,
                                                             xlabel=self.__label_ratio + "by run ratio",
                                                             ylabel="Counts",
                                                             # title='Detectors Ratios Histogram (lumi weighted)',
                                                             xmin=setts.ratio_min_reduced,
                                                             xmax=setts.ratio_max_reduced,
                                                             mean=self.__by_run_ratios_mean,
                                                             stdv=self.__by_run_ratios_stdv,
                                                             energy_year_label=self.__year_energy_label)
        self.__plt_plots['by_run_hist'] = ratio_hist_lumi2_w[0][0].get_figure()

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

    # Fill test
    def plot_bad_fills(self):
        ratio_vs_fill = plotting.plot_bad_fill_info(data_frame=self.common_data_filtered_no_nan, x_data_label='fill',
                                                    y_data_label=self.label_ratio,
                                                    z_data_label=self.__label_det2_lumi_to_use,
                                                    xlabel='Fill Number', ylabel=self.__label_ratio + " ratios",
                                                    mean=self.__nls_ratios_lw_mean,
                                                    stdv=self.__nls_ratios_lw_stdv,
                                                    ratio_acceptance=setts.bad_ratio_to_plot_stdv_factor,
                                                    filePath=self.output_dir + 'txt/',
                                                    ymin=self.__nls_ratios_lw_mean - 3 * setts.bad_ratio_to_plot_stdv_factor
                                                         * self.__nls_ratios_lw_stdv,
                                                    ymax=self.__nls_ratios_lw_mean + 3 * setts.bad_ratio_to_plot_stdv_factor
                                                         * self.__nls_ratios_lw_stdv,
                                                    energy_year_label=self.__year_energy_label,
                                                    txtfileName='Bad_fills')
        self.__plt_plots['by_ratio_vs_fill_bad_fills'] = ratio_vs_fill

    def plot_bad_runs(self):
        ratio_vs_run = plotting.plot_bad_fill_info(data_frame=self.common_data_filtered_no_nan, x_data_label='run',
                                                   y_data_label=self.label_ratio,
                                                   z_data_label=self.__label_det2_lumi_to_use,
                                                   xlabel='Run Number', ylabel=self.__label_ratio + " ratios",
                                                   mean=self.__nls_ratios_lw_mean,
                                                   stdv=self.__nls_ratios_lw_stdv,
                                                   ratio_acceptance=setts.bad_ratio_to_plot_stdv_factor,
                                                   filePath=self.output_dir + 'txt/',
                                                   ymin=self.__nls_ratios_lw_mean - 3 * setts.bad_ratio_to_plot_stdv_factor * self.__nls_ratios_lw_stdv,
                                                   ymax=self.__nls_ratios_lw_mean + 3 * setts.bad_ratio_to_plot_stdv_factor * self.__nls_ratios_lw_stdv,
                                                   energy_year_label=self.__year_energy_label,
                                                   txtfileName='Bad_runs')
        self.__plt_plots['by_ratio_vs_run_bad_runs'] = ratio_vs_run

    # All/excluded data analysis plots
    def __plot_all_and_excluded(self, det_index):
        if type(det_index) == int and det_index <= 1:
            det = self.__detcs[det_index]
            suffix_plot_name = det.excluded_label.replace(" ", "_").replace(":", "")
            scatter_all_and_excluded_lumi2 = plotting.snsplot_detector_all_and_excluded(self.__data_all,
                                                                                        x_data_label=self.accumulated_lumi2_label,
                                                                                        y_data_label=self.label_ratio,
                                                                                        conditional_label=det.excluded_label,
                                                                                        xlabel="Integrated luminosity [$" +
                                                                                               self.lumi_unit + "^{-1}$]",
                                                                                        ylabel=self.__label_ratio + " ratios",
                                                                                        ymin=setts.ratio_min,
                                                                                        ymax=setts.ratio_max,
                                                                                        energy_year_label=self.__year_energy_label)
            self.__sns_plots['all_and_excluded_vs_lum2_' + suffix_plot_name] = scatter_all_and_excluded_lumi2

            scatter_all_and_excluded_run = plotting.snsplot_detector_all_and_excluded(self.__data_all,
                                                                                      x_data_label="run",
                                                                                      y_data_label=self.label_ratio,
                                                                                      conditional_label=det.excluded_label,
                                                                                      xlabel="Run",
                                                                                      ylabel=self.__label_ratio + " ratios",
                                                                                      ymin=setts.ratio_min,
                                                                                      ymax=setts.ratio_max,
                                                                                      energy_year_label=self.__year_energy_label)
            self.__sns_plots['all_and_excluded_vs_run_' + suffix_plot_name] = scatter_all_and_excluded_run
        else:
            raise AssertionError("Index must be 0 or 1")

    def plot_all_and_excluded_by_detc(self, single_detector_info_plots=False, all_ls_plots=True):
        if single_detector_info_plots:
            self.__plot_all_and_excluded(0)
            self.__plot_all_and_excluded(1)

        if all_ls_plots:
            scatter_all_and_excluded_run_comb = plotting.snsplot_detector_all_and_excluded(self.__data_all,
                                                                                           x_data_label="run",
                                                                                           y_data_label=self.label_ratio,
                                                                                           conditional_label=self.__ratio_excluded_label,
                                                                                           # conditional_label_extra=self.det2.excluded_label,
                                                                                           xlabel="Run",
                                                                                           ylabel=self.__label_ratio + " ratios",
                                                                                           ymin=setts.ratio_min,
                                                                                           ymax=setts.ratio_max,
                                                                                           energy_year_label=self.__year_energy_label)
            self.__sns_plots['all_and_excluded_vs_run_combined'] = scatter_all_and_excluded_run_comb

            scatter_all_and_excluded_vs_lumi2_comb = plotting.snsplot_detector_all_and_excluded(self.__data_all,
                                                                                                x_data_label=self.accumulated_lumi2_label,
                                                                                                y_data_label=self.label_ratio,
                                                                                                conditional_label=self.__ratio_excluded_label,
                                                                                                # conditional_label_extra=self.det2.excluded_label,
                                                                                                xlabel="Integrated luminosity [$" +
                                                                                                       self.lumi_unit + "^{-1}$]",
                                                                                                ylabel=self.__label_ratio + " ratios",
                                                                                                ymin=setts.ratio_min,
                                                                                                ymax=setts.ratio_max,
                                                                                                energy_year_label=self.__year_energy_label)
            self.__sns_plots['all_and_excluded_vs_lumi2_combined'] = scatter_all_and_excluded_vs_lumi2_comb

            scatter_all_and_excluded_vs_date_comb = plotting.snsplot_detector_all_and_excluded(self.__data_all,
                                                                                               x_data_label='date',
                                                                                               y_data_label=self.label_ratio,
                                                                                               conditional_label=self.__ratio_excluded_label,
                                                                                               # conditional_label_extra=self.det2.excluded_label,
                                                                                               xlabel="Date",
                                                                                               ylabel=self.__label_ratio + " ratios",
                                                                                               ymin=setts.ratio_min,
                                                                                               ymax=setts.ratio_max,
                                                                                               xmin=
                                                                                               self.__data_all['date'][
                                                                                                   0],
                                                                                               xmax=
                                                                                               self.__data_all['date'][
                                                                                                   len(
                                                                                                       self.__data_all) - 1],
                                                                                               xlabel_rotation=30,
                                                                                               energy_year_label=self.__year_energy_label)
            self.__sns_plots['all_and_excluded_vs_date_combined'] = scatter_all_and_excluded_vs_date_comb

        scatter_all_and_excluded_vs_date_comb_by_nls = plotting.snsplot_detector_all_and_excluded(self.__data_all,
                                                                                                  x_data_label='date',
                                                                                                  y_data_label=self.by_nls_label_ratio,
                                                                                                  conditional_label=self.by_nls_exclusion_info_label,
                                                                                                  # conditional_label_extra=self.det2.excluded_label,
                                                                                                  xlabel="Date",
                                                                                                  ylabel=self.__label_ratio + " ratios in "
                                                                                                         + str(
                                                                                                      self.__nls) + ' LS',
                                                                                                  ymin=setts.ratio_min,
                                                                                                  ymax=setts.ratio_max,
                                                                                                  xmin=self.__data_all[
                                                                                                      'date'][0],
                                                                                                  xmax=self.__data_all[
                                                                                                      'date'][len(
                                                                                                      self.__data_all) - 1],
                                                                                                  xlabel_rotation=30,
                                                                                                  energy_year_label=self.__year_energy_label)
        self.__sns_plots['all_and_excluded_by_nls_vs_date_combined'] = scatter_all_and_excluded_vs_date_comb_by_nls

        scatter_all_and_excluded_vs_lumi2_comb_by_nls = plotting.snsplot_detector_all_and_excluded(self.__data_all,
                                                                                                   x_data_label=self.accumulated_lumi2_label,
                                                                                                   y_data_label=self.by_nls_label_ratio,
                                                                                                   conditional_label=self.by_nls_exclusion_info_label,
                                                                                                   # conditional_label_extra=self.det2.excluded_label,
                                                                                                   xlabel="Integrated luminosity [$" +
                                                                                                          self.lumi_unit + "^{-1}$]",
                                                                                                   ylabel=self.__label_ratio + " ratios in "
                                                                                                          + str(
                                                                                                       self.__nls) + ' LS',
                                                                                                   ymin=setts.ratio_min,
                                                                                                   ymax=setts.ratio_max,
                                                                                                   energy_year_label=self.__year_energy_label)
        self.__sns_plots['all_and_excluded_by_nls_vs_lumi2_combined'] = scatter_all_and_excluded_vs_lumi2_comb_by_nls

        scatter_all_and_excluded_percent_vs_lumi2_comb_by_nls = plotting.snsplot_detector_all_and_excluded(
            self.__data_all,
            x_data_label=self.accumulated_lumi2_label,
            y_data_label=self.by_nls_label_ratio,
            conditional_label=self.by_nls_exclusion_percent_label,
            # conditional_label_extra=self.det2.excluded_label,
            xlabel="Integrated luminosity [$" +
                   self.lumi_unit + "^{-1}$]",
            ylabel=self.__label_ratio + " ratios in "
                   + str(
                self.__nls) + ' LS',
            ymin=setts.ratio_min,
            ymax=setts.ratio_max,
            energy_year_label=self.__year_energy_label)
        self.__sns_plots[
            'all_and_excluded_percent_by_nls_vs_lumi2_combined'] = scatter_all_and_excluded_percent_vs_lumi2_comb_by_nls

        scatter_all_and_excluded_vs_run_comb_by_nls = plotting.snsplot_detector_all_and_excluded(self.__data_all,
                                                                                                 x_data_label='run',
                                                                                                 y_data_label=self.by_nls_label_ratio,
                                                                                                 conditional_label=self.by_nls_exclusion_info_label,
                                                                                                 # conditional_label_extra=self.det2.excluded_label,
                                                                                                 xlabel="Run",
                                                                                                 ylabel=self.__label_ratio + " ratios in "
                                                                                                        + str(
                                                                                                     self.__nls) + ' LS',
                                                                                                 ymin=setts.ratio_min,
                                                                                                 ymax=setts.ratio_max,
                                                                                                 energy_year_label=self.__year_energy_label)
        self.__sns_plots['all_and_excluded_by_nls_vs_run_combined'] = scatter_all_and_excluded_vs_run_comb_by_nls

        # all/excluded histograms
        scatter_all_and_excluded_hist = plotting.snsplot_hist_all_and_excluded(self.__data_all,
                                                                               x_data_label=self.label_ratio,
                                                                               conditional_label=self.ratio_excluded_label,
                                                                               # conditional_label_extra=self.det2.excluded_label,
                                                                               xlabel=self.__label_ratio + " ratios",
                                                                               ylabel="Counts",
                                                                               xmin=setts.ratio_min,
                                                                               xmax=setts.ratio_max,
                                                                               bins=setts.nbins,
                                                                               energy_year_label=self.__year_energy_label)
        self.__sns_plots['all_and_excluded_hist'] = scatter_all_and_excluded_hist

        scatter_all_and_excluded_hist_by_nls = plotting.snsplot_hist_all_and_excluded(self.__data_all,
                                                                                      x_data_label=self.by_nls_label_ratio,
                                                                                      conditional_label=self.by_nls_exclusion_info_label,
                                                                                      # conditional_label_extra=self.det2.excluded_label,
                                                                                      xlabel=self.__label_ratio + " ratios in " + str(
                                                                                          self.__nls) + ' LS',
                                                                                      ylabel="Counts",
                                                                                      xmin=setts.ratio_min,
                                                                                      xmax=setts.ratio_max,
                                                                                      bins=setts.nbins,
                                                                                      energy_year_label=self.__year_energy_label)
        self.__sns_plots['all_and_excluded_by_nls_hist'] = scatter_all_and_excluded_hist_by_nls

    def save_plots(self):
        if self.__only_stats:
            print("Only stats mode. Plots won't be produced")
        else:
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
        return self.__by_nls_lumi_label1

    @property
    def by_nls_lumi2(self):
        return self.__by_nls_lumi_label2

    @property
    def by_nls_lumi_label(self):
        return self.__by_nls_lumi_label

    @property
    def by_nls_lumi_label(self):
        return self.__by_nls_lumi_label

    @property
    def by_nls_label_ratio_norm(self):
        return self.__by_nls_label_ratio_norm

    @property
    def det1(self):
        return self.__det1

    @property
    def det2(self):
        return self.__det2

    @property
    def by_nls_accumulated_lumi1_label(self):
        return self.__by_nls_accumulated_lumi1_label

    @property
    def by_nls_accumulated_lumi2_label(self):
        return self.__by_nls_accumulated_lumi2_label

    @property
    def accumulated_lumi1_label(self):
        return self.__accumulated_lumi1_label

    @property
    def accumulated_lumi2_label(self):
        return self.__accumulated_lumi2_label

    @property
    def lumi2_label(self):
        return self.__label_det2_lumi_to_use

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
    def data_all(self):
        if self.__data_all is None:
            raise AssertionError("please run first the create_nls_data() method")
        return self.__data_all

    @property
    def data_exclusion_percent(self):
        return self.__data_exclusion_percent

    @property
    def by_nls_exclusion_info_label(self):
        return self.__by_nls_exclusion_info_label

    @property
    def by_nls_binary_exclusion_info_label(self):
        return self.__by_nls_binary_exclusion_info_label

    @property
    def by_nls_exclusion_percent_label(self):
        return self.__by_nls_exclusion_percent_label

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
    def nls_ratios_lw_mean(self):
        if self.__nls_ratios_lw_mean is None:
            raise AssertionError("please run first the create_nls_data() method")
        return self.__nls_ratios_lw_mean

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

    @property
    def year_energy_label(self):
        return self.__year_energy_label

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