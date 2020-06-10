import numpy as np
import pandas as pd
from tools import lumi_tools as ltools
import settings as setts
from datetime import datetime


class Luminometer:
    __allowed_detectors = ["PXL", "HFOC", "HFET", "DT", "BCM1F", "PLTZERO", "RAMSES", "BCM1FSI", "ZMonitoring",
                           "ZCounts", "ZCountsBB", "ZCountsBB_MC", "ZCountsEE_MC", "ZCounts_MC", "ZCountsEE", "zcounts"]
    __allowed_detectors_file_labels = ["pcc", "hfoc", "hfet", "dt", "bcm1f", "plt", "ram", "bcm1f_si", "bcm1f_d"]

    __lumi_unit = 'fb'
    __nbx_label = 'nbx'

    # _lumi_unit_label = "r'$"+ __lumi_unit +"^{-1}$'"

    def __init__(self, name: str, data_file_name: str, full_data_analysis: bool = False,
                 mixed_data=False, remove_extra_cols: bool = True, split_run_fill_ls_cols: bool = True) -> None:
        # name contains detector name+possible calibration tag. Example:pcc18v3
        self.name = name.upper()
        self.__all_data_analysis_included = full_data_analysis
        self.__mixed_data = mixed_data
        self.__lumi_unit = Luminometer.__lumi_unit
        self.__lumi_rec_label = 'lumi_recorded_' + self.name
        self.__lumi_del_label = 'lumi_delivered_' + self.name
        self.__lumi_rec_label_original_units = 'lumi_recorded_ou_' + self.name
        self.__lumi_del_label_original_units = 'lumi_delivered_ou_' + self.name
        self.__detector_name_label = 'detector_name_' + self.name
        self.__excluded_label = self.name + ' exclusion info:'
        self.__label_for_excluded = 'excluded'
        self.__label_for_included = 'included'

        __csv_file_col_names = ['run:fill', 'ls', 'time', 'beamstatus', 'E(GeV)', self.__lumi_del_label_original_units,
                                self.__lumi_rec_label_original_units, 'avgpu', self.detector_name_label]
        __not_needed_columns = ['run:fill', 'beamstatus', 'E(GeV)', 'avgpu']

        # read data from csv file
        # data_file_pd = pd.read_csv(data_file_name, comment='#', index_col=False, names=__csv_file_col_names)

        # Get number of commented lines and find header commented line
        rows_to_skip_in_csv = setts.rows_to_skip_in_csv
        rows_to_skip_at_the_end = setts.rows_to_skip_at_the_end

        self.__rec_info_in_file = False
        self.__del_info_in_file = False

        data_file_pd = pd.read_csv(data_file_name, index_col=False, skiprows=rows_to_skip_in_csv,
                                   skipfooter=rows_to_skip_at_the_end, engine='python')

        if 'delivered(hz/ub)' in data_file_pd.columns and 'recorded(hz/ub)' in data_file_pd.columns:
            self.__rec_info_in_file = True
            self.__del_info_in_file = True
            data_file_pd.rename(columns={'#run:fill': 'run:fill',
                                         'delivered(hz/ub)': self.__lumi_del_label_original_units,
                                         'recorded(hz/ub)': self.__lumi_rec_label_original_units,
                                         'source': self.detector_name_label}, inplace=True)
        elif 'recorded(hz/ub)' in data_file_pd.columns:
            print("**** WARNING: Reading only recorded. delivered lumi not found in the input file!")
            self.__rec_info_in_file = True
            data_file_pd.rename(columns={'#run:fill': 'run:fill',
                                         'recorded(hz/ub)': self.__lumi_rec_label_original_units,
                                         'source': self.detector_name_label}, inplace=True)
        elif 'delivered(hz/ub)' in data_file_pd.columns:
            print("**** WARNING: Reading only recorded. recorded lumi not found in the input file!")
            self.__del_info_in_file = True
            data_file_pd.rename(columns={'#run:fill': 'run:fill',
                                         'delivered(hz/ub)': self.__lumi_del_label_original_units,
                                         'source': self.detector_name_label}, inplace=True)
        else:
            raise AssertionError("*** ERROR: No luminosity values found. At least a column with delivered(hz/ub) "
                                 "or recorded(hz/ub) must be included in the input file!")

        run_fill_cols = data_file_pd["run:fill"].str.split(":", n=1, expand=True)
        ls_double = data_file_pd["ls"].str.split(":", n=1, expand=True)
        if split_run_fill_ls_cols:
            data_file_pd["run"] = run_fill_cols[0].astype(str).astype(int)
            data_file_pd["fill"] = run_fill_cols[1].astype(str).astype(int)
            data_file_pd["ls"] = ls_double[0].astype(str).astype(int)
        if remove_extra_cols:
            try:
                data_file_pd.drop(columns=__not_needed_columns, inplace=True)
            except:
                print("** Warning! Extra columns not removed!")

        year, energy = ltools.get_year_and_energy(run_fill_cols[1][0])
        self.year = year
        self.energy = energy

        # Convert lumi to [1/fb], including multiplication by 23.3s
        # self.__lumi_csv_unit = ltools.get_lumi_unit_from_csv(data_file_name)
        self.__lumi_csv_unit = 'hz/ub'
        lumi_convertion_factor = ltools.get_lumi_unit_convertion_factor_to_inv_fb(self.__lumi_csv_unit)

        columns_names_to_remove = []
        if self.__del_info_in_file:
            data_file_pd[self.__lumi_del_label] = data_file_pd[
                                                      self.__lumi_del_label_original_units] * lumi_convertion_factor
            columns_names_to_remove.append(self.__lumi_del_label_original_units)
        if self.__rec_info_in_file:
            data_file_pd[self.__lumi_rec_label] = data_file_pd[
                                                      self.__lumi_rec_label_original_units] * lumi_convertion_factor
            columns_names_to_remove.append(self.__lumi_rec_label_original_units)

        if remove_extra_cols:
            data_file_pd.drop(columns=columns_names_to_remove, inplace=True)

        # TODO: Find if timestamp is for UTC or some CERN time
        # ltools.add_date_column(data_file_pd)

        self.__data = data_file_pd
        self.__check_detector_type()

        self.__all_data = None
        self.__excluded_data = None

        print('Initialized detector from file: ' + data_file_name)

        if full_data_analysis:
            # data_file_pd_all = pd.read_csv(data_file_name.replace(".csv", "_all.csv"), comment='#', index_col=False, names=__csv_file_col_names)

            data_file_pd_all = pd.read_csv(data_file_name.replace(".csv", "_all.csv"), index_col=False,
                                           skiprows=rows_to_skip_in_csv, skipfooter=rows_to_skip_at_the_end,
                                           engine='python')
            if 'delivered(hz/ub)' in data_file_pd_all.columns and 'recorded(hz/ub)' in data_file_pd_all.columns:
                data_file_pd_all.rename(columns={'#run:fill': 'run:fill',
                                                 'delivered(hz/ub)': self.__lumi_del_label_original_units,
                                                 'recorded(hz/ub)': self.__lumi_rec_label_original_units,
                                                 'source': self.detector_name_label}, inplace=True)
            elif 'recorded(hz/ub)' in data_file_pd_all.columns:
                print("**** WARNING: Reading only recorded. delivered lumi not found in the input file!")
                data_file_pd_all.rename(columns={'#run:fill': 'run:fill',
                                                 'delivered(hz/ub)': self.__lumi_del_label_original_units,
                                                 'recorded(hz/ub)': self.__lumi_rec_label_original_units,
                                                 'source': self.detector_name_label}, inplace=True)
            else:
                raise AssertionError("*** ERROR: No luminosity values found. At least a column with delivered(hz/ub) "
                                     "or recorded(hz/ub) must be included in the input file!")

            run_fill_cols_all = data_file_pd_all["run:fill"].str.split(":", n=1, expand=True)
            ls_double_all = data_file_pd_all["ls"].str.split(":", n=1, expand=True)
            data_file_pd_all["run"] = run_fill_cols_all[0].astype(str).astype(int)
            data_file_pd_all["fill"] = run_fill_cols_all[1].astype(str).astype(int)
            data_file_pd_all["ls"] = ls_double_all[0].astype(str).astype(int)
            data_file_pd_all.drop(columns=__not_needed_columns, inplace=True)

            data_file_pd_all[self.__lumi_del_label] = data_file_pd_all[
                                                          self.__lumi_del_label_original_units] * lumi_convertion_factor
            data_file_pd_all[self.__lumi_rec_label] = data_file_pd_all[
                                                          self.__lumi_rec_label_original_units] * lumi_convertion_factor

            self.__all_data = data_file_pd_all
            self.__all_data = self.__all_data.sort_values(by="time", ascending=True)

            print('Initialized detector from file: ' + data_file_name.replace(".csv", "_all.csv"))

            merge_tmp_excluded = data_file_pd_all[~data_file_pd_all['time'].isin(data_file_pd['time'])]
            self.__all_data[self.__excluded_label]= ~data_file_pd_all['time'].isin(data_file_pd['time'])

            # Changing name from True -> excluded; False -> included
            self.__all_data[self.__excluded_label].replace([True, False],
                                                           [self.__label_for_excluded, self.__label_for_included],
                                                           inplace=True)

            merge_tmp_excluded = merge_tmp_excluded.reset_index(drop=True)

            self.__excluded_data = merge_tmp_excluded

            if abs(len(data_file_pd_all) - len(data_file_pd)) != len(merge_tmp_excluded):
                raise AssertionError("Problem during extraction of excluded data")

            # print(len(merge_tmp_excluded), len(data_file_pd), len(data_file_pd_all))
            print("    Data percent excluded in normtag (%): " +
                  str((len(data_file_pd_all) - len(data_file_pd))*100/len(data_file_pd_all)))

            self.__excluded_data = self.__excluded_data.sort_values(by="time", ascending=True)

        self.__data = self.__data.sort_values(by="time", ascending=True)

    # Check detector type consistency between file_name, csv file label and allowed detectors
    def __check_detector_type(self):
        # check detector consistency:
        if not self.__mixed_data:
            label_from_csv_list = np.unique(self.__data[self.detector_name_label])
            if len(label_from_csv_list) > 1:
                raise AssertionError('More than one detector in .csv file. Detectors in csv: ' + str(label_from_csv_list))
            elif len(label_from_csv_list) == 1:
                label_from_csv = label_from_csv_list[0]
                if label_from_csv not in Luminometer.__allowed_detectors:
                    raise AssertionError('detector not allowed. Allowed detectors: ' +
                                         str(Luminometer.__allowed_detectors))
    @property
    def name(self):
        return self.__name

    @property
    def all_data_analysis_included(self):
        return self.__all_data_analysis_included

    @property
    def data(self):
        return self.__data

    @property
    def all_data(self):
        return self.__all_data

    @property
    def excluded_data(self):
        return self.__excluded_data

    @property
    def lumi_unit(self):
        return self.__lumi_unit

    @property
    def year(self):
        return self.__year

    @property
    def energy(self):
        return self.__energy

    @property
    def lumi_rec_label(self):
        return self.__lumi_rec_label

    @property
    def lumi_del_label(self):
        return self.__lumi_del_label

    @property
    def lumi_rec_label_original_units(self):
        return self.__lumi_rec_label_original_units

    @property
    def lumi_del_label_original_units(self):
        return self.__lumi_del_label_original_units

    @property
    def detector_name_label(self):
        return self.__detector_name_label

    @property
    def nbx_label(self):
        return Luminometer.__nbx_label

    @property
    def excluded_label(self):
        return self.__excluded_label

    @property
    def label_for_included(self):
        return self.__label_for_included

    @property
    def label_for_excluded(self):
        return self.__label_for_excluded

    @name.setter
    def name(self, name):
        self.__name = name

    @data.setter
    def data(self, data):
        self.__data = data

    @year.setter
    def year(self, year):
        valid_years = np.array(list(setts.rangeFillYearly))[:, 0]
        if int(year) in valid_years:
            self.__year = year
        else:
            raise ValueError("Year not implemented")

    @energy.setter
    def energy(self, energy):
        valid_energies = np.array(list(setts.rangeFillYearly))[:, 1]
        if int(energy) in valid_energies:
            self.__energy = energy
        else:
            raise ValueError("Energy not implemented")
