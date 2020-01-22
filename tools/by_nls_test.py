from typing import Dict, Any

from tools.detectorsratio import DetectorsRatio as Ratios
import tools.plotting_tools as plotting
from tools.luminometer import Luminometer as L
from tools.linearityanalysis import Stats as StatNls
import tools.lumi_tools as ltools
import settings as setts
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.api as sm


class ByNlsTest:
    __hist_name = 'nls_ratio_hist_lw'

    def __init__(self, det1: L, det2: L, array_by_step=None) -> None:
        self.__year_energy_label = None
        self.__plt_plots = {}
        self.__data = []
        self.__by_nls_lumi_label = None
        self.__label_ratio = None
        self.__by_nls_label_ratio = None

        self.__all_nls_data = pd.DataFrame()
        self.__all_nls_weights = pd.DataFrame()

        if array_by_step is not None:
            self.__nls = [1]
            for i in range(1, 150 // array_by_step):
                self.__nls.append(i * array_by_step)
            nls_array = self.__nls
        else:
            nls_array = setts.nls_list_to_test

        self.__nls = nls_array

        for nls_iterator in nls_array:
            ratios12 = Ratios(det1, det2, nls=nls_iterator)

            self.__all_nls_data[str(nls_iterator)] = ratios12.common_data_filtered[ratios12.by_nls_label_ratio]
            self.__all_nls_weights[str(nls_iterator)] = ratios12.common_data_filtered[ratios12.by_nls_lumi_label]

            # Fill needed methods
            ratios12.fill_stats()

            self.__data.append([nls_iterator, ratios12.nls_ratios_stdv])
            self.__year_energy_label = ratios12.year_energy_label
            self.__output_dir = ratios12.output_dir
            self.__by_nls_lumi_label = ratios12.by_nls_lumi_label
            self.__label_ratio = ratios12.label_ratio
            self.__by_nls_label_ratio = ratios12.by_nls_label_ratio

        self.__data_file_pd = pd.DataFrame(self.__data, columns=['nls', 'stddev'])
        self.__output_dir = self.__output_dir + 'stdtest/'

    def plot_stddev_vs_nls(self):
        stddev_vs_nls = plotting.scatter_plot_from_pandas_frame(data_frame=self.__data_file_pd,
                                                                y_data_label='stddev',
                                                                x_data_label='nls',
                                                                xlabel='NLS',
                                                                ylabel='Ratios Standard Deviation',
                                                                energy_year_label=self.__year_energy_label,
                                                                plot_style='c',
                                                                marker_size=1)
        self.__plt_plots['by_nls_stddev_vs_nls'] = stddev_vs_nls.get_figure()

    def plot_nls_histograms(self):
        print('creating nls histograms ...')
        hist_list = plotting.hist_list_from_pandas_frame(self.__all_nls_data, 200,
                                                         xlabel=self.__label_ratio + " ratios in " + ' LS',
                                                         ylabel="Integrated luminosity [$" + "fb" + "^{-1}$]",
                                                         xmin=setts.ratio_min, xmax=setts.ratio_max,
                                                         energy_year_label=self.__year_energy_label,
                                                         all_nls_weights=self.__all_nls_weights,
                                                         legend_labels=self.__nls, legend_position='upper right',
                                                         fill=False)
        self.__plt_plots['Several_nls_histograms'] = hist_list

    def save_plots(self):
        print('\n\n Saving test plots:')
        plotting.save_plots(self.__plt_plots, self.__output_dir)
