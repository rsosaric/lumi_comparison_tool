from tools.detectorsratio import DetectorsRatio as Ratios
import tools.plotting_tools as plotting
import tools.lumi_tools as ltools
import settings as setts
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import statsmodels.api as sm


def lin_func(x, a, b):
    return a * x + b


class LinearityAnalysis:
    __conversion_fb_to_ub = 1000000000.0 / setts.xLS
    __linearity_lumi_conversion = __conversion_fb_to_ub

    # TODO: by fill fitting
    # TODO: by fill slopes histogram (mean, std)
    def __init__(self, ratios: Ratios) -> None:
        print('\n Initializing ' + ratios.label_ratio + ' Linearity Analysis')
        self.__plt_plots = []
        self.__output_dir = ratios.output_dir + 'linearity_analysis/'
        self.__by_fill_fits_dir = self.__output_dir + '/by_fill_fits/'
        # Check and load NBX information
        if ratios.nbx_data is None:
            ratios.fill_nbx()

        ini_data = ratios.common_data_filtered
        self.__ratios = ratios

        # compute SBIL
        self.__sbil_label = 'sbil'
        self.__nls_sbil_label = 'nls_sbil'
        self.__nbx_label = ratios.det2.nbx_label
        self.__nls_ratio_label = ratios.by_nls_label_ratio
        self.__by_nls_label_ratio_err = ratios.by_nls_label_ratio_err
        self.__nls = setts.nls_for_lin_year[ratios.year, ratios.energy]
        ini_data[self.__sbil_label] = \
            ini_data[ratios.det2.lumi_rec_label] \
            * LinearityAnalysis.__linearity_lumi_conversion / ini_data[self.__nbx_label]
        ini_data[self.__nls_sbil_label] = \
            ini_data[ratios.by_nls_lumi_label] \
            * LinearityAnalysis.__linearity_lumi_conversion / (ini_data[self.__nbx_label] * ratios.nls)

        self.__sbil_mean = ini_data[self.__sbil_label].mean()
        self.__sbil_stdv = ini_data[self.__sbil_label].std()

        self.__lin_data = ini_data

        self.__year_energy_label = ratios.year_energy_label
        self.__label_ratio = ratios.label_ratio

        self.__by_fill_data = {}
        self.__by_fill_stats = {}
        self.__by_fill_slopes = []
        self.__by_fill_slopes_err = []
        self.__by_fill_nls_slopes = []
        self.__by_fill_nls_slopes_err = []
        self.__by_fill_slopes_mean = None
        self.__by_fill_slopes_stdv = None
        self.__by_fill_nls_slopes_mean = None
        self.__by_fill_nls_slopes_stdv = None
        self.__excluded_fills_in_running = []

        self.create_by_fill_data()

    def create_by_fill_data(self):
        print('Creating by fill data ...')
        ltools.check_and_create_folder(self.__by_fill_fits_dir, creation_info=False)
        fills = self.__ratios.fills
        fills.sort()
        # organize data by single fill
        for fill in fills:
            exclude_fill = False
            self.__by_fill_data[fill] = self.__lin_data[self.__lin_data['fill'] == fill]
            # get stats for each fill
            self.__by_fill_stats[fill] = Stats(self.__by_fill_data[fill])
            # remove outliers
            if setts.remove_outliers_in_lin_by_fill:
                mean_ratio = self.__by_fill_stats[fill].mean[self.__label_ratio]
                stdv_ratio = self.__by_fill_stats[fill].stdv[self.__label_ratio]
                min_ratio_val = mean_ratio - setts.allowed_ratio_stdv_factor*stdv_ratio
                max_ratio_val = mean_ratio + setts.allowed_ratio_stdv_factor*stdv_ratio

                mean_sbil = self.__by_fill_stats[fill].mean[self.__sbil_label]
                stdv_sbil = self.__by_fill_stats[fill].stdv[self.__sbil_label]
                min_sbil_val = mean_sbil - setts.allowed_ratio_stdv_factor*stdv_sbil
                max_sbil_val = mean_sbil + setts.allowed_ratio_stdv_factor*stdv_sbil

                self.__by_fill_data[fill] = self.__by_fill_data[fill][(self.__by_fill_data[fill][self.__label_ratio] >= min_ratio_val) &
                                                                      (self.__by_fill_data[fill][self.__label_ratio] <= max_ratio_val) &
                                                                      (self.__by_fill_data[fill][self.__sbil_label] >= min_sbil_val) &
                                                                      (self.__by_fill_data[fill][self.__sbil_label] <= max_sbil_val)]
            if len(self.__by_fill_data[fill][self.__label_ratio]) < setts.min_number_of_points_req or \
                    len(self.__by_fill_data[fill][self.__nls_ratio_label].dropna()) < setts.min_number_of_points_req:
                self.__excluded_fills_in_running.append(fill)
                exclude_fill = True


            # Fitting
            if exclude_fill != True:
                self.fill_fitting(fill)

        self.__by_fill_slopes = np.array(self.__by_fill_slopes)
        self.__by_fill_slopes_mean = self.__by_fill_slopes.mean()
        self.__by_fill_slopes_stdv = self.__by_fill_slopes.std()

        self.__by_fill_nls_slopes = np.array(self.__by_fill_nls_slopes)
        self.__by_fill_nls_slopes_mean = np.array(self.__by_fill_nls_slopes).mean()
        self.__by_fill_nls_slopes_stdv = np.array(self.__by_fill_nls_slopes).std()

        slope_hist = plotting.hist_from_array(self.__by_fill_slopes,
                                              nbins=setts.nbins,
                                              xlabel="slope [hz/" + r'$\mu$' + "b]",
                                              ylabel='Counts',
                                              mean=self.__by_fill_slopes_mean,
                                              stdv=self.__by_fill_slopes_stdv
                                              )
        # TODO: increase point size in plots
        slope_hist_nls = plotting.hist_from_array(self.__by_fill_nls_slopes,
                                                  nbins=setts.nbins,
                                                  xlabel="slope [hz/" + r'$\mu$' + "b]",
                                                  ylabel='Counts',
                                                  mean=self.__by_fill_nls_slopes_mean,
                                                  stdv=self.__by_fill_nls_slopes_stdv
                                                  )

        self.__plt_plots.append(('slope_hist', slope_hist))
        self.__plt_plots.append(('slope_hist_nls', slope_hist_nls))

        # TODO: plot slopes vs Fill
        # TODO: plot slopes vs Lumi

    def fill_fitting(self, fill: int) -> None:
        fill_data = self.__by_fill_data[fill]
        # TODO: check dropna effect
        # TODO: restrict sbil range for fitting using the histogram
        # TODO: get chi2 values or R value
        # TODO: slope histos, mean, stdv, err (weigthed)
        x_nls = fill_data[self.__nls_sbil_label].dropna()
        y_nls = fill_data[self.__nls_ratio_label].dropna()
        ey_nls = fill_data[self.__by_nls_label_ratio_err].dropna()
        x = fill_data[self.__sbil_label].dropna()
        y = fill_data[self.__label_ratio].dropna()
        print(fill)
        popt, pcov = curve_fit(lin_func, x, y)
        slope = popt[0]
        slope_err = np.sqrt(pcov[0, 0])
        intercept = popt[1]
        # intercept_err = np.sqrt(pcov[1, 1])

        nls_popt, nls_pcov = curve_fit(lin_func, x_nls, y_nls, sigma=ey_nls)
        nls_slope = nls_popt[0]
        nls_slope_err = np.sqrt(nls_pcov[0, 0])
        nls_intercept = nls_popt[1]
        # nls_intercept_err = np.sqrt(nls_pcov[1, 1])

        fig = plotting.plot_from_fit(x, y, fitted_slope=slope, fitted_intercept=intercept,
                                     fitted_slope_err=slope_err,
                                     fitted_f=lin_func, xlabel="SBIL [hz/" + r'$\mu$' + "b]",
                                     ylabel=self.__label_ratio + " ratio",
                                     energy_year_label=self.__year_energy_label)
        nls_fig = plotting.plot_from_fit(x_nls, y_nls, y_err=ey_nls,
                                         fitted_slope=nls_slope,
                                         fitted_slope_err=nls_slope_err,
                                         fitted_intercept=nls_intercept,
                                         fitted_f=lin_func,
                                         xlabel="SBIL [hz/" + r'$\mu$' + "b]",
                                         ylabel=self.__label_ratio + " ratios in " +
                                                str(self.__nls) + ' LS',
                                         energy_year_label=self.__year_energy_label)

        fig.savefig(self.__by_fill_fits_dir + str(fill))
        nls_fig.savefig(self.__by_fill_fits_dir + str(fill) + "_nls")
        plt.close(fig)
        plt.close(nls_fig)

        self.__by_fill_slopes.append(slope)
        self.__by_fill_slopes_err.append(slope_err)
        self.__by_fill_nls_slopes.append(nls_slope)
        self.__by_fill_slopes_err.append(nls_slope_err)

    def plot_hist_sbil(self):
        sbil_hist = plotting.hist_from_pandas_frame(data_frame=self.__lin_data,
                                                    col_label=self.__sbil_label,
                                                    nbins=setts.nbins,
                                                    xlabel=self.__ratios.det2.name +
                                                           " SBIL [hz/" + r'$\mu$' + "b]",
                                                    ylabel='Counts',
                                                    # title='by Nls Detectors Ratios Histogram',
                                                    # xmin=setts.ratio_min, xmax=setts.ratio_max,
                                                    mean=self.__sbil_mean, stdv=self.__sbil_stdv,
                                                    energy_year_label=self.__year_energy_label)
        self.__plt_plots.append(['sbil_hist', sbil_hist[0][0].get_figure()])

    # TODO: create same for nls organized and rescaled data
    def plot_ratio_vs_all_data_sbil(self):
        ratio_vs_all_data_sbil = plotting.scatter_plot_from_pandas_frame(data_frame=self.__lin_data,
                                                                         y_data_label=self.__label_ratio,
                                                                         x_data_label=self.__sbil_label,
                                                                         xlabel="SBIL [hz/" + r'$\mu$' + "b]",
                                                                         ylabel=self.__label_ratio + " ratios",
                                                                         ymin=setts.ratio_min, ymax=setts.ratio_max,
                                                                         energy_year_label=self.__year_energy_label,
                                                                         fig_size_shape='sq')
        self.__plt_plots.append(['ratio_vs_all_data_sbil', ratio_vs_all_data_sbil.get_figure()])

    # def plot_linear_fit_by_fill(self, fill, model_fit):
    #     fill_fit_plot = plotting.plot_from_fit(self.__by_fill_data[fill], model_fit=model_fit,
    #                                            x_data_label=self.__nls_sbil_label,
    #                                            y_data_label=self.__nls_ratio_label)
    #     self.__plt_plots.append([self.__by_fill_fits_dir_name + 'fill_' + str(fill), fill_fit_plot.get_figure()])

    def save_plots(self):
        print('\n\n Saving linearity plots for ' + self.__label_ratio + ':')
        plotting.save_plots(self.__plt_plots, self.__output_dir)

    @property
    def sbil_mean(self):
        return self.__sbil_mean

    @property
    def sbil_stdv(self):
        return self.__sbil_stdv


class Stats:
    def __init__(self, array: pd.DataFrame) -> None:
        self.__mean = array.mean()
        self.__stdv = array.std()

    @property
    def mean(self):
        return self.__mean

    @property
    def stdv(self):
        return self.__stdv
