import numpy as np
import pandas as pd
from tools import lumi_tools as ltools
import settings as setts
from datetime import datetime


class Luminometer:
    __allowed_detectors = ["PXL", "HFOC", "HFET", "DT", "BCM1F", "PLTZERO", "RAMSES"]
    __allowed_detectors_file_labels = ["pcc", "hfoc", "hfet", "dt", "bcm1f", "plt", "ram"]

    __lumi_unit = 'fb'
    __nbx_label = 'nbx'

    # _lumi_unit_label = "r'$"+ __lumi_unit +"^{-1}$'"

    def __init__(self, name: str, data_file_name: str, all_data: bool = False) -> None:
        self.name = name
        if all_data:
            self.__name = name
        self.__lumi_unit = Luminometer.__lumi_unit
        self.__lumi_rec_label = 'lumi_recorded_' + self.name
        self.__lumi_del_label = 'lumi_delivered_' + self.name
        self.__detector_name_label = 'detector_name_' + self.name

        __csv_file_col_names = ['run_fill', 'ls', 'time', 'beam_status', 'energy', self.lumi_del_label,
                                self.lumi_rec_label, 'avg_PU_' + self.name, self.detector_name_label]
        __not_needed_columns = ['run_fill', 'beam_status', 'energy', 'avg_PU_' + self.name]

        # read data from csv file
        data_file_pd = pd.read_csv(data_file_name, comment='#', index_col=False, names=__csv_file_col_names)
        run_fill_cols = data_file_pd["run_fill"].str.split(":", n=1, expand=True)
        ls_double = data_file_pd["ls"].str.split(":", n=1, expand=True)
        data_file_pd["run"] = run_fill_cols[0].astype(str).astype(int)
        data_file_pd["fill"] = run_fill_cols[1].astype(str).astype(int)
        data_file_pd["ls"] = ls_double[0].astype(str).astype(int)
        data_file_pd.drop(columns=__not_needed_columns, inplace=True)

        year, energy = ltools.get_year_and_energy(run_fill_cols[1][0])
        self.year = year
        self.energy = energy

        # Convert lumi to [1/fb], including multiplication by 23.3s
        self.__lumi_csv_unit = ltools.get_lumi_unit_from_csv(data_file_name)
        lumi_convertion_factor = ltools.get_lumi_unit_convertion_factor_to_inv_fb(self.__lumi_csv_unit)
        data_file_pd['lumi_delivered_' + self.name] = data_file_pd[
                                                          'lumi_delivered_' + self.name] * lumi_convertion_factor
        data_file_pd['lumi_recorded_' + self.name] = data_file_pd['lumi_recorded_' + self.name] * lumi_convertion_factor

        # TODO: Find if timestamp is for UTC or some CERN time
        # ltools.add_date_column(data_file_pd)

        self.__data = data_file_pd

        print('Initialized detector from file: ' + data_file_name)

    @property
    def name(self):
        return self.__name

    @property
    def data(self):
        return self.__data

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
    def detector_name_label(self):
        return self.__detector_name_label

    @property
    def nbx_label(self):
        return Luminometer.__nbx_label

    @name.setter
    def name(self, name):
        if name in Luminometer.__allowed_detectors_file_labels:
            self.__name = name
        else:
            raise ValueError("Luminometer not implemented")

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
