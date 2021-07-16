from tools.detectorsratio import DetectorsRatio as Ratios
from tools.luminometer import Luminometer as L
from scipy.stats.distributions import chi2 as func_chi2
import tools.plotting_tools as plotting
import tools.lumi_tools as ltools
import settings as setts
import numpy as np
import pandas as pd
import math
import json
from scipy.optimize import curve_fit
from lmfit import Model
import matplotlib.pyplot as plt


def lin_func(x, a, b):
    return a * x + b


lin_model = Model(lin_func)


class LinearityAnalysis:
    __conversion_fb_to_ub = 1000000000.0 / setts.xLS
    __linearity_lumi_conversion = __conversion_fb_to_ub

    def __init__(self, dets_file_labels: list, input_dir: str, lumi_type: str = 'rec') -> None:

        self.__label_col_fill = "fill"
        self.__label_col_slopes = "slopes"
        self.__label_col_slopes_err = "err on the slope"
        self.__label_col_slopes_chi2 = "chi2 on the slope"
        self.__label_col_nls_slopes = "by nls slopes"
        self.__label_col_nls_slopes_err = "by nls err on the slope"
        self.__label_col_nls_slopes_rel_err = "by nls relative err on the slope"
        self.__label_col_nls_slopes_redchi2 = "by nls reduced chi2 on the slope"
        self.__label_col_nls_slopes_pvalue = "by nls chi2 pvalue on the slope"
        self.__label_col_fill_lumi = "lumi in fill"
        self.__label_col_accumulated_lumi = "accumulated lumi up to fill"
        self.__label_col_total_w = "total weight"

        n_files = len(dets_file_labels)
        detcs = []
        for det in dets_file_labels:
            try:
                detcs.append(L(det, input_dir + det + ".csv"))

            except IOError as errIO:
                print(errIO)
                print('Please check if default input folder is correctly created: ' + setts.csv_input_base_dir)
                print('Also check that your .csv file is in the correct year folder: ' + input_dir)
                raise

        self.__year = detcs[0].year
        self.__energy = detcs[0].energy

        # Trying to read nls from config file
        try:
            self.__nls = setts.nls_for_lin_year[int(self.__year), int(self.__energy)]
            print("By nls ratios will be computed using nls=" + str(self.__nls))
        except AssertionError:
            print("no Nls configuration found for " + str(self.__year) + " [" + str(self.__energy) + "tev]")

        # Trying to read (min,max) slopes from config file
        try:
            self.__min_allowed_slope = setts.limits_linearity_slopes_year[int(self.__year), int(self.__energy)][0]
            self.__max_allowed_slope = setts.limits_linearity_slopes_year[int(self.__year), int(self.__energy)][1]
            print("Slopes range have been set to: (" + str(self.__min_allowed_slope) + "," +
                  str(self.__max_allowed_slope) + ")")
        except:
            print("no (min,max) slopes configuration found for " + str(self.__year) + " [" + str(
                self.__energy) + "tev] in the settings.py")
            print("Setting no limits for the slopes!")
            self.__min_allowed_slope = None
            self.__max_allowed_slope = None

        ratios = Ratios(detcs[0], detcs[1], lumi_type=lumi_type, nls=self.__nls)

        assert (self.__nls == ratios.nls)

        print('\n Initializing ' + ratios.label_ratio + ' Linearity Analysis')
        self.__plt_plots = {}
        self.__output_dir = ratios.output_dir + setts.linearity_analysis_output_dir
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
        self.__lumi_unit = ratios.lumi_unit
        self.__by_nls_label_ratio_err = ratios.by_nls_label_ratio_err
        ini_data[self.__sbil_label] = ini_data[ratios.det2.lumi_rec_label] \
                                      * LinearityAnalysis.__linearity_lumi_conversion / ini_data[self.__nbx_label]
        ini_data[self.__nls_sbil_label] = \
            ini_data[ratios.by_nls_lumi_label] \
            * LinearityAnalysis.__linearity_lumi_conversion / (ini_data[self.__nbx_label] * self.__nls)

        self.__sbil_mean = ini_data[self.__sbil_label].mean()
        self.__sbil_stdv = ini_data[self.__sbil_label].std()

        self.__lin_data = ini_data

        self.__year_energy_label = ratios.year_energy_label
        self.__label_ratio = ratios.label_ratio
        self.__label_norm_by_fill_ratio = self.__label_ratio + "_by_fill_norm"
        self.__label_norm_by_nls_by_fill_ratio = self.__label_ratio + "by_nls_by_fill_norm"
        self.__label_lumi = ratios.lumi2_label
        self.__accumulated_lumi2_label = ratios.accumulated_lumi2_label

        # Contains the data for each fill
        self.__by_fill_data = {}
        # Contains the data properties (mean, slope, etc) for each fill
        self.__by_fill_df = pd.DataFrame()
        # mean and by nls mean for each fill
        self.__by_fill_means = {}
        self.__label_mean_by_fill = self.__label_ratio + "_mean_by_fill"
        self.__label_by_nls_mean_by_fill = self.__label_ratio + "by_nls_mean_by_fill"
        self.__by_fill_df_all = pd.DataFrame()
        self.__by_fill_lumi = []
        self.__by_fill_acumulative_lumi = []
        self.__by_fill_stats = {}
        self.__by_fill_slopes = []
        self.__by_fill_slopes_err = []
        self.__by_fill_nls_slopes = []
        self.__by_fill_nls_slopes_err = []
        self.__by_fill_nls_slopes_redchi2 = []
        self.__by_fill_nls_slopes_pvalue = []
        self.__by_fill_slopes_mean = None
        self.__by_fill_slopes_stdv = None
        self.__by_fill_nls_slopes_mean = None
        self.__by_fill_nls_slopes_stdv = None
        self.__by_fill_slopes_mean_lw = None
        self.__by_fill_slopes_stdv_lw = None
        self.__by_fill_nls_slopes_mean_lw = None
        self.__by_fill_nls_slopes_stdv_lw = None
        self.__by_fill_nls_slopes_mean_errw = None
        self.__by_fill_nls_slopes_stdv_errw = None
        self.__by_fill_nls_slopes_mean_totalw = None
        self.__by_fill_nls_slopes_stdv_totalw = None
        self.__by_fill_nls_slopes_err_totalw = None
        self.__by_fill_nls_slopes_err_lw = None
        self.__all_data_fitted_slope = None
        self.__all_data_fitted_slope_err = None
        self.__all_data_fitted_slope_chi2 = None
        self.__excluded_fills_in_running = []
        self.__good_fills = []
        self.__average_points_in_summary = setts.points_in_summary_lin_year[self.__year, self.__energy]
        self.__data_avg_summary = pd.DataFrame()

        self.__by_run_data = {}
        self.__by_run_mean = {}
        self.__by_run_std = {}

        self.create_by_run_data()
        self.create_by_fill_data()
        self.get_normalize_by_fill_data()
        # self.by_nls_outlier_removal()
        self.fit_all_data()
        self.get_slope_average_per_intervals()
        # print(self.__by_fill_df)

        self.plot_hist_sbil()
        self.plot_ratio_vs_all_data_sbil()
        self.plot_by_fill_norm_ratio_vs_all_data_sbil()
        self.plot_by_nls_by_fill_norm_ratio_vs_all_data_sbil()
        self.save_plots()
        self.save_summary_csv()

    def create_by_fill_data(self):
        print('Creating by fill data ...')
        ltools.check_and_create_folder(self.__by_fill_fits_dir, creation_info=False)
        fills = self.__ratios.fills
        good_fills = []
        fills.sort()
        # organize data by single fill
        for fill in fills:
            exclude_fill = False
            self.__by_fill_data[fill] = self.__lin_data[self.__lin_data['fill'] == fill]
            # get stats for each fill
            self.__by_fill_stats[fill] = Stats(self.__by_fill_data[fill])
            mean_ratio = self.__by_fill_stats[fill].mean[self.__label_ratio]
            mean_by_nls_ratio = self.__by_fill_stats[fill].mean[self.__nls_ratio_label]
            stdv_ratio = self.__by_fill_stats[fill].stdv[self.__label_ratio]

            self.__by_fill_means[fill] = {self.__label_mean_by_fill: mean_ratio,
                                          self.__label_by_nls_mean_by_fill: mean_by_nls_ratio}
            # Processing fill data previous the fitting  ===============================================================
            # remove outliers
            if setts.remove_outliers_in_lin_by_fill:
                min_ratio_val = mean_ratio - setts.allowed_ratio_stdv_factor * stdv_ratio
                max_ratio_val = mean_ratio + setts.allowed_ratio_stdv_factor * stdv_ratio

                mean_sbil = self.__by_fill_stats[fill].mean[self.__sbil_label]
                stdv_sbil = self.__by_fill_stats[fill].stdv[self.__sbil_label]
                min_sbil_val = mean_sbil - setts.allowed_ratio_stdv_factor * stdv_sbil
                max_sbil_val = mean_sbil + setts.allowed_ratio_stdv_factor * stdv_sbil

                self.__by_fill_data[fill] = self.__by_fill_data[fill][
                    (self.__by_fill_data[fill][self.__label_ratio] >= min_ratio_val) &
                    (self.__by_fill_data[fill][self.__label_ratio] <= max_ratio_val) &
                    (self.__by_fill_data[fill][self.__sbil_label] >= min_sbil_val) &
                    (self.__by_fill_data[fill][self.__sbil_label] <= max_sbil_val)]
            if len(self.__by_fill_data[fill][self.__label_ratio]) < setts.min_number_of_points_req or \
                    len(self.__by_fill_data[fill][self.__nls_ratio_label].dropna()) < setts.min_number_of_points_req:
                self.__excluded_fills_in_running.append(fill)
                exclude_fill = True

            # Fitting  ======================================================================================
            if not exclude_fill:
                fitted_data_ini = len(self.__by_fill_slopes)
                self.fill_fitting(fill)
                fitted_data_changed = len(self.__by_fill_slopes)

                if fitted_data_changed == fitted_data_ini + 1:
                    good_fills.append(fill)
                    self.__by_fill_acumulative_lumi.append(
                        np.array(self.__by_fill_data[fill][self.__accumulated_lumi2_label])[-1])
                    self.__by_fill_lumi.append(np.sum(self.__by_fill_data[fill][self.__label_lumi]))
                else:
                    print("     -->> Fill " + str(fill) + " not included in final results.")
            else:
                print("     -->> Fill " + str(fill) + " not included in final results.")

        # Filling fitting summary dataframe =================================================================
        if len(good_fills) == len(self.__by_fill_lumi) and len(good_fills) == len(self.__by_fill_slopes):
            self.__good_fills = good_fills
            self.__by_fill_df[self.__label_col_fill] = good_fills
            self.__by_fill_df[self.__label_col_fill_lumi] = self.__by_fill_lumi
            self.__by_fill_df[self.__label_col_accumulated_lumi] = self.__by_fill_acumulative_lumi
            self.__by_fill_df[self.__label_col_slopes] = self.__by_fill_slopes
            self.__by_fill_df[self.__label_col_slopes_err] = self.__by_fill_slopes_err
            self.__by_fill_df[self.__label_col_nls_slopes] = self.__by_fill_nls_slopes
            self.__by_fill_df[self.__label_col_nls_slopes_err] = self.__by_fill_nls_slopes_err
            self.__by_fill_df[self.__label_col_nls_slopes_redchi2] = self.__by_fill_nls_slopes_redchi2
            self.__by_fill_df[self.__label_col_nls_slopes_pvalue] = self.__by_fill_nls_slopes_pvalue

            self.__by_fill_df.replace([np.inf, -np.inf], np.nan, inplace=True)
            self.__by_fill_df.dropna(inplace=True)

            self.__by_fill_df[self.__label_col_total_w] = ltools.get_total_weights(
                self.__by_fill_df[self.__label_col_fill_lumi],
                self.__by_fill_df[self.__label_col_nls_slopes_err])
        else:
            raise AssertionError("ERROR in create_by_fill_data(): lengths of good_fills, by_fill_slopes and "
                                 "by_fill_lumi do not match! ")

        # Overall fit quality analysis (by nls fits) =================================================================
        # fit error #
        temp_mean_slope = self.__by_fill_df[self.__label_col_nls_slopes].mean()
        self.__by_fill_df[self.__label_col_nls_slopes_rel_err] \
            = abs(self.__by_fill_df[self.__label_col_nls_slopes_err])/abs(temp_mean_slope)

        mean_of_errors = self.__by_fill_df[self.__label_col_nls_slopes_err].mean()
        stdv_of_errors = self.__by_fill_df[self.__label_col_nls_slopes_err].std()
        mean_of_rel_errors = self.__by_fill_df[self.__label_col_nls_slopes_rel_err].mean()
        stdv_of_rel_errors = self.__by_fill_df[self.__label_col_nls_slopes_rel_err].std()
        min_rel_err = np.min(self.__by_fill_df[self.__label_col_nls_slopes_rel_err])
        max_rel_err = np.max(self.__by_fill_df[self.__label_col_nls_slopes_rel_err])
        ltools.color_print("\n === Fit quality analysis ===", "blue")
        print("*-> Mean of the fitted slope errors (relative errors): " + str(float("{0:.3f}".format(mean_of_errors))) +
              " (" + str(float("{0:.3f}".format(mean_of_rel_errors))) + ")")
        print("*-> Stdv of the fitted slope errors (relative errors): " + str(float("{0:.3f}".format(stdv_of_errors))) +
              " (" + str(float("{0:.3f}".format(stdv_of_rel_errors))) + ")")
        mean_of_redchi2s = self.__by_fill_df[self.__label_col_nls_slopes_redchi2].mean()
        stdv_of_redchi2s = self.__by_fill_df[self.__label_col_nls_slopes_redchi2].std()
        min_redchi2 = np.min(self.__by_fill_df[self.__label_col_nls_slopes_redchi2])
        max_redchi2 = np.max(self.__by_fill_df[self.__label_col_nls_slopes_redchi2])
        print("*-> Mean of the fitted slope chi2s: " + str(float("{0:.3f}".format(mean_of_redchi2s))))
        print("*-> Stdv of the fitted slope chi2s: " + str(float("{0:.3f}".format(stdv_of_redchi2s))))

        # Setting thresholds
        redchi2_good_limit = mean_of_redchi2s + stdv_of_redchi2s
        rel_errors_good_limit = mean_of_rel_errors + stdv_of_rel_errors

        self.__by_fill_df_all = self.__by_fill_df.copy()
        all_slopes_n = len(self.__by_fill_df)

        # drop slopes with bad fit quality
        self.__by_fill_df = self.__by_fill_df[
            (self.__by_fill_df[self.__label_col_nls_slopes_rel_err] <= rel_errors_good_limit) &
            (self.__by_fill_df[self.__label_col_nls_slopes_redchi2] <= redchi2_good_limit)
            ].copy()
        after_excl_slopes_n = len(self.__by_fill_df)
        print(str((all_slopes_n - after_excl_slopes_n) * 100 / all_slopes_n) +
              " % of the fills have been excluded because they had very bad fit quality.")

        # drop slopes outside the range in settings
        temp_slopes_n = len(self.__by_fill_df)
        if self.__min_allowed_slope is not None and self.__max_allowed_slope is not None:
            self.__by_fill_df = self.__by_fill_df[
                (self.__by_fill_df[self.__label_col_slopes] >= self.__min_allowed_slope) &
                (self.__by_fill_df[self.__label_col_slopes] <= self.__max_allowed_slope) &
                (self.__by_fill_df[self.__label_col_nls_slopes] >= self.__min_allowed_slope) &
                (self.__by_fill_df[self.__label_col_nls_slopes] <= self.__max_allowed_slope)]
            after_excl_slopes_n = len(self.__by_fill_df)
            self.__by_fill_df.reset_index(drop=True, inplace=True)
            print(str((temp_slopes_n - after_excl_slopes_n) * 100 / temp_slopes_n) +
                  " % of the fills have been excluded because they had slopes outside the allowed range.")

        ltools.color_print("PLEASE CHECK IF THESE NUMBERS ARE BIG!", "yellow")

        self.__by_fill_slopes_mean = self.__by_fill_df[self.__label_col_slopes].mean()
        self.__by_fill_slopes_stdv = self.__by_fill_df[self.__label_col_slopes].std()

        self.__by_fill_nls_slopes_mean = self.__by_fill_df[self.__label_col_nls_slopes].mean()
        self.__by_fill_nls_slopes_stdv = self.__by_fill_df[self.__label_col_nls_slopes].std()

        # Getting weighted values:
        lw_stats = ltools.get_w_stats(self.__by_fill_df[self.__label_col_slopes],
                                      self.__by_fill_df[self.__label_col_fill_lumi],
                                      min_val=self.__min_allowed_slope, max_val=self.__max_allowed_slope)
        lw_nls_stats = ltools.get_w_stats(self.__by_fill_df[self.__label_col_nls_slopes],
                                          self.__by_fill_df[self.__label_col_fill_lumi],
                                          min_val=self.__min_allowed_slope, max_val=self.__max_allowed_slope)
        errw_nls_stats = ltools.get_w_stats(self.__by_fill_df[self.__label_col_nls_slopes],
                                            self.__by_fill_df[self.__label_col_nls_slopes_err],
                                            min_val=self.__min_allowed_slope, max_val=self.__max_allowed_slope)
        totalw_nls_stats = ltools.get_w_stats(self.__by_fill_df[self.__label_col_nls_slopes],
                                              self.__by_fill_df[self.__label_col_total_w],
                                              min_val=self.__min_allowed_slope, max_val=self.__max_allowed_slope)

        self.__by_fill_slopes_mean_lw = lw_stats.mean
        self.__by_fill_slopes_stdv_lw = lw_stats.std
        self.__by_fill_nls_slopes_mean_lw = lw_nls_stats.mean
        self.__by_fill_nls_slopes_stdv_lw = lw_nls_stats.std
        self.__by_fill_nls_slopes_mean_errw = errw_nls_stats.mean
        self.__by_fill_nls_slopes_stdv_errw = errw_nls_stats.std
        self.__by_fill_nls_slopes_mean_totalw = totalw_nls_stats.mean
        self.__by_fill_nls_slopes_stdv_totalw = totalw_nls_stats.std
        self.__by_fill_nls_slopes_err_totalw = self.__by_fill_nls_slopes_stdv_totalw/math.sqrt(len(self.__by_fill_df))
        self.__by_fill_nls_slopes_err_lw = self.__by_fill_nls_slopes_stdv_lw / math.sqrt(len(self.__by_fill_df))

        slope_hist = plotting.hist_from_array(self.__by_fill_df[self.__label_col_slopes],
                                              nbins=setts.nbins_linearity,
                                              xlabel="slope [hz/" + r'$\mu$' + "b]",
                                              ylabel='Counts',
                                              mean=self.__by_fill_slopes_mean,
                                              stdv=self.__by_fill_slopes_stdv,
                                              xmin=self.__min_allowed_slope,
                                              xmax=self.__max_allowed_slope,
                                              mean_float_format="{0:.4f}"
                                              )

        slope_hist_nls = plotting.hist_from_array(self.__by_fill_df[self.__label_col_nls_slopes],
                                                  nbins=setts.nbins_linearity,
                                                  xlabel="slope [hz/" + r'$\mu$' + "b]",
                                                  ylabel='Counts',
                                                  mean=self.__by_fill_nls_slopes_mean,
                                                  stdv=self.__by_fill_nls_slopes_stdv,
                                                  xmin=self.__min_allowed_slope,
                                                  xmax=self.__max_allowed_slope,
                                                  mean_float_format="{0:.4f}"
                                                  )

        slope_hist_lw = plotting.hist_from_pandas_frame(data_frame=self.__by_fill_df,
                                                        col_label=self.__label_col_slopes,
                                                        nbins=setts.nbins_linearity,
                                                        xlabel="slope [hz/" + r'$\mu$' + "b]",
                                                        ylabel="Integrated luminosity [$" +
                                                               self.__ratios.lumi_unit + "^{-1}$]",
                                                        xmin=self.__min_allowed_slope,
                                                        xmax=self.__max_allowed_slope,
                                                        mean=self.__by_fill_slopes_mean_lw,
                                                        stdv=self.__by_fill_slopes_stdv_lw,
                                                        energy_year_label=self.__year_energy_label,
                                                        weight_label=self.__label_col_fill_lumi,
                                                        mean_float_format="{0:.4f}")

        slope_hist_nls_lw = plotting.hist_from_pandas_frame(data_frame=self.__by_fill_df,
                                                            col_label=self.__label_col_nls_slopes,
                                                            nbins=setts.nbins_linearity,
                                                            xlabel="slope [hz/" + r'$\mu$' + "b]",
                                                            ylabel="Integrated luminosity [$" +
                                                                   self.__ratios.lumi_unit + "^{-1}$]",
                                                            xmin=self.__min_allowed_slope,
                                                            xmax=self.__max_allowed_slope,
                                                            mean=self.__by_fill_nls_slopes_mean_lw,
                                                            stdv=self.__by_fill_nls_slopes_stdv_lw,
                                                            energy_year_label=self.__year_energy_label,
                                                            weight_label=self.__label_col_fill_lumi,
                                                            mean_float_format="{0:.4f}")

        slope_hist_nls_errw = plotting.hist_from_pandas_frame(data_frame=self.__by_fill_df,
                                                              col_label=self.__label_col_nls_slopes,
                                                              nbins=setts.nbins_linearity,
                                                              xlabel="slope [hz/" + r'$\mu$' + "b]",
                                                              ylabel="Integrated luminosity [$" +
                                                                     self.__ratios.lumi_unit + "^{-1}$]",
                                                              xmin=self.__min_allowed_slope,
                                                              xmax=self.__max_allowed_slope,
                                                              mean=self.__by_fill_nls_slopes_mean_errw,
                                                              stdv=self.__by_fill_nls_slopes_stdv_errw,
                                                              energy_year_label=self.__year_energy_label,
                                                              weight_label=self.__label_col_slopes_err,
                                                              mean_float_format="{0:.4f}",
                                                              normlize=True)

        slope_hist_nls_totalw = plotting.hist_from_pandas_frame(data_frame=self.__by_fill_df,
                                                                col_label=self.__label_col_nls_slopes,
                                                                nbins=setts.nbins_linearity,
                                                                xlabel="slope [hz/" + r'$\mu$' + "b]",
                                                                ylabel="events weighted by lumi and fit error",
                                                                xmin=self.__min_allowed_slope,
                                                                xmax=self.__max_allowed_slope,
                                                                mean=self.__by_fill_nls_slopes_mean_totalw,
                                                                stdv=self.__by_fill_nls_slopes_stdv_totalw,
                                                                energy_year_label=self.__year_energy_label,
                                                                weight_label=self.__label_col_total_w,
                                                                mean_float_format="{0:.4f}")

        slope_nls_vs_fill = plotting.plot_scatter_and_errors(data_frame=self.__by_fill_df,
                                                             y_data_label=self.__label_col_nls_slopes,
                                                             x_data_label="fill",
                                                             err_y_data_label=self.__label_col_nls_slopes_err,
                                                             ymin=self.__min_allowed_slope,
                                                             ymax=self.__max_allowed_slope,
                                                             xlabel="fill number",
                                                             ylabel="slope [hz/" + r'$\mu$' + "b]",
                                                             energy_year_label=self.__year_energy_label,
                                                             use_integers_in_x_axis=True,
                                                             fig_size_shape="sq")
        slope_nls_vs_lumi = plotting.plot_scatter_and_errors(data_frame=self.__by_fill_df,
                                                             y_data_label=self.__label_col_nls_slopes,
                                                             x_data_label=self.__label_col_accumulated_lumi,
                                                             err_y_data_label=self.__label_col_nls_slopes_err,
                                                             ymin=self.__min_allowed_slope,
                                                             ymax=self.__max_allowed_slope,
                                                             xlabel="Integrated luminosity [$" +
                                                                    self.__lumi_unit + "^{-1}$]",
                                                             ylabel="slope [hz/" + r'$\mu$' + "b]",
                                                             energy_year_label=self.__year_energy_label,
                                                             fig_size_shape="sq")

        slope_err_hist_nls = plotting.hist_from_pandas_frame(data_frame=self.__by_fill_df_all,
                                                             col_label=self.__label_col_nls_slopes_err,
                                                             nbins=200,
                                                             xlabel="slope error [hz/" + r'$\mu$' + "b]",
                                                             ylabel="Counts",
                                                             mean=mean_of_errors,
                                                             stdv=stdv_of_errors,
                                                             energy_year_label=self.__year_energy_label,
                                                             mean_float_format="{0:.4f}")

        slope_rel_err_hist_nls = plotting.hist_from_pandas_frame(data_frame=self.__by_fill_df_all,
                                                                 col_label=self.__label_col_nls_slopes_rel_err,
                                                                 nbins=200,
                                                                 xlabel="slope error [hz/" + r'$\mu$' + "b]",
                                                                 ylabel="Counts",
                                                                 mean=mean_of_rel_errors,
                                                                 stdv=stdv_of_rel_errors,
                                                                 xmin=min_rel_err, xmax=max_rel_err,
                                                                 energy_year_label=self.__year_energy_label,
                                                                 mean_float_format="{0:.4f}")

        slope_redchi2_hist_nls = plotting.hist_from_pandas_frame(data_frame=self.__by_fill_df_all,
                                                                 col_label=self.__label_col_nls_slopes_redchi2,
                                                                 nbins=200,
                                                                 xlabel="chi2/dof",
                                                                 ylabel="Counts",
                                                                 mean=mean_of_redchi2s,
                                                                 stdv=stdv_of_redchi2s,
                                                                 xmin=min_redchi2, xmax=max_redchi2,
                                                                 energy_year_label=self.__year_energy_label,
                                                                 mean_float_format="{0:.4f}")

        self.__plt_plots['slope_hist'] = slope_hist
        self.__plt_plots['slope_hist_nls'] = slope_hist_nls

        self.__plt_plots['slope_err_hist_nls'] = slope_err_hist_nls[0][0].get_figure()
        self.__plt_plots['slope_rel_err_hist_nls'] = slope_rel_err_hist_nls[0][0].get_figure()
        self.__plt_plots['slope_redchi2_hist_nls'] = slope_redchi2_hist_nls[0][0].get_figure()

        self.__plt_plots['slope_hist_lw'] = slope_hist_lw[0][0].get_figure()
        self.__plt_plots['slope_hist_nls_lw'] = slope_hist_nls_lw[0][0].get_figure()
        self.__plt_plots['slope_hist_nls_errw'] = slope_hist_nls_errw[0][0].get_figure()
        self.__plt_plots['slope_hist_nls_totalw'] = slope_hist_nls_totalw[0][0].get_figure()

        self.__plt_plots['slope_nls_vs_fill'] = slope_nls_vs_fill
        self.__plt_plots['slope_nls_vs_lumi'] = slope_nls_vs_lumi

        # TODO: plot slopes vs Fill
        # TODO: plot slopes vs Lumi

    def create_by_run_data(self):
        print('Creating by run data ...')
        runs = self.__ratios.runs
        runs.sort()

        for run in runs:
            exclude_fill = False
            self.__by_run_data[run] = self.__lin_data[self.__lin_data['run'] == run]
            # get stats for each fill
            stats = Stats(self.__by_run_data[run])
            self.__by_run_mean[run] = stats.mean[self.__label_ratio]
            self.__by_run_std[run] = stats.stdv[self.__label_ratio]

    def get_normalize_by_fill_data(self):
        print('Normalizing all fills data ...')
        norm_ratios = []
        norm_nls_ratios = []

        for index_data in range(0, len(self.__lin_data)):
            fill = self.__lin_data['fill'][index_data]

            mean = self.__by_fill_means[fill][self.__label_mean_by_fill]
            mean_nls = self.__by_fill_means[fill][self.__label_by_nls_mean_by_fill]
            norm_ratios.append(self.__lin_data[self.__label_ratio][index_data]/mean)
            norm_nls_ratios.append(self.__lin_data[self.__nls_ratio_label][index_data]/mean_nls)

        self.__lin_data[self.__label_norm_by_fill_ratio] = np.array(norm_ratios)
        self.__lin_data[self.__label_norm_by_nls_by_fill_ratio] = np.array(norm_nls_ratios)

    # do this in detectors ratios
    # def correct_by_run_instabilities(self):

    def by_nls_outlier_removal(self):
        mean_nls_err = self.__lin_data[self.__by_nls_label_ratio_err].mean()
        std_nls_err = self.__lin_data[self.__by_nls_label_ratio_err].std()

        self.__lin_data = self.__lin_data[
            (self.__lin_data[self.__by_nls_label_ratio_err] >= mean_nls_err - std_nls_err) &
            (self.__lin_data[self.__by_nls_label_ratio_err] <= mean_nls_err + std_nls_err)].copy()

    def fit_all_data(self) -> None:
        if len(self.__lin_data[self.__label_norm_by_nls_by_fill_ratio]) == 0:
            raise AssertionError("norm_by_nls_by_fill_ratio data is needed")

        to_fit_data_nls = pd.DataFrame()

        to_fit_data_nls["x_nls"] = self.__lin_data[self.__nls_sbil_label]
        to_fit_data_nls["y_nls"] = self.__lin_data[self.__label_norm_by_nls_by_fill_ratio]
        to_fit_data_nls["ey_nls"] = self.__lin_data[self.__by_nls_label_ratio_err]

        to_fit_data_nls.dropna(inplace=True)

        # Extra exclusion for by nls
        min_ratio_val = setts.ratio_min
        max_ratio_val = setts.ratio_max
        to_fit_data_nls = to_fit_data_nls[
            (to_fit_data_nls["y_nls"] >= min_ratio_val) &
            (to_fit_data_nls["y_nls"] <= max_ratio_val) &
            (to_fit_data_nls["ey_nls"] > 0.0)]

        x = to_fit_data_nls["x_nls"]
        y = to_fit_data_nls["y_nls"]
        ey = to_fit_data_nls["ey_nls"]

        result_fit = lin_model.fit(y, x=x, a=1, b=1, weights=1.0/ey)

        slope = result_fit.params['a'].value
        slope_err = result_fit.params['a'].stderr
        intercept = result_fit.params['b'].value
        chi2 = result_fit.redchi

        self.__all_data_fitted_slope = slope
        self.__all_data_fitted_slope_err = slope_err
        self.__all_data_fitted_slope_chi2 = chi2

        nls_fig = plotting.plot_from_fit(x, y, y_err=ey,
                                         fitted_slope=slope, chi2=chi2,
                                         fitted_slope_err=slope_err,
                                         fitted_intercept=intercept,
                                         fitted_f=lin_func,
                                         ymin=0.85, ymax=1.15,
                                         xlabel="SBIL [hz/" + r'$\mu$' + "b]",
                                         ylabel=self.__label_ratio + " ratios in " +
                                                str(self.__nls) + ' LS',
                                         energy_year_label=self.__year_energy_label,
                                         add_linearity_special_text=True)

        nls_fig.savefig(self.__output_dir + "all_by_nls_data_fit")
        plt.close(nls_fig)

    def fill_fitting(self, fill: int) -> None:
        print("Processing fill: " + str(fill))
        fill_data = self.__by_fill_data[fill]
        # TODO: get chi2 values or R value

        to_fit_data = pd.DataFrame()
        to_fit_data_nls = pd.DataFrame()

        to_fit_data["x"] = fill_data[self.__sbil_label]
        to_fit_data["y"] = fill_data[self.__label_ratio]
        to_fit_data_nls["x_nls"] = fill_data[self.__nls_sbil_label]
        to_fit_data_nls["y_nls"] = fill_data[self.__nls_ratio_label]
        to_fit_data_nls["ey_nls"] = fill_data[self.__by_nls_label_ratio_err]

        to_fit_data.dropna(inplace=True)
        to_fit_data_nls.dropna(inplace=True)

        # Extra exclusion for by nls
        min_ratio_val = setts.ratio_min
        max_ratio_val = setts.ratio_max
        to_fit_data_nls = to_fit_data_nls[
            (to_fit_data_nls["y_nls"] >= min_ratio_val) &
            (to_fit_data_nls["y_nls"] <= max_ratio_val) &
            (to_fit_data_nls["ey_nls"] > 0.0)]

        x_nls = to_fit_data_nls["x_nls"]
        y_nls = to_fit_data_nls["y_nls"]
        ey_nls = to_fit_data_nls["ey_nls"]
        x = to_fit_data["x"]
        y = to_fit_data["y"]

        # popt, pcov = curve_fit(lin_func, x, y)
        result_fit = lin_model.fit(y, x=x, a=1, b=1)

        slope = result_fit.params['a'].value
        slope_err = result_fit.params['a'].stderr
        intercept = result_fit.params['b'].value
        chi2 = result_fit.redchi

        assert (len(x_nls) == len(y_nls) == len(ey_nls))
        if len(x_nls) >= setts.min_num_data_to_fit:
            # nls_popt, nls_pcov = curve_fit(lin_func, x_nls, y_nls, sigma=ey_nls)
            # nls_slope = nls_popt[0]
            # nls_slope_err = np.sqrt(nls_pcov[0, 0])
            # nls_intercept = nls_popt[1]
            # nls_intercept_err = np.sqrt(nls_pcov[1, 1])

            result_nls_fit = lin_model.fit(y_nls, x=x_nls, a=1, b=1, weights=1.0 / ey_nls)
            # print(result.fit_report())
            nls_slope = result_nls_fit.params['a'].value
            nls_slope_err = result_nls_fit.params['a'].stderr
            nls_intercept = result_nls_fit.params['b'].value
            # reduced chi2 = chi2/dof, more info here: https://lmfit.github.io/lmfit-py/fitting.html
            nls_chi2 = result_nls_fit.redchi
            nls_full_chi2 = result_nls_fit.chisqr
            nls_dof = result_nls_fit.nfree
            p_value = func_chi2(nls_full_chi2, nls_dof)

            fig = plotting.plot_from_fit(x, y, fitted_slope=slope, fitted_intercept=intercept,
                                         fitted_slope_err=slope_err, chi2=chi2,
                                         fitted_f=lin_func, xlabel="SBIL [hz/" + r'$\mu$' + "b]",
                                         ylabel=self.__label_ratio + " ratio",
                                         energy_year_label=self.__year_energy_label,
                                         add_linearity_special_text=True)

            nls_fig = plotting.plot_from_fit(x_nls, y_nls, y_err=ey_nls,
                                             fitted_slope=nls_slope, chi2=nls_chi2,
                                             fitted_slope_err=nls_slope_err,
                                             fitted_intercept=nls_intercept,
                                             fitted_f=lin_func,
                                             xlabel="SBIL [hz/" + r'$\mu$' + "b]",
                                             ylabel=self.__label_ratio + " ratios in " +
                                                    str(self.__nls) + ' LS',
                                             energy_year_label=self.__year_energy_label,
                                             add_linearity_special_text=True)

            fig.savefig(self.__by_fill_fits_dir + str(fill))
            nls_fig.savefig(self.__by_fill_fits_dir + str(fill) + "_nls")
            plt.close(fig)
            plt.close(nls_fig)

            self.__by_fill_slopes.append(slope)
            self.__by_fill_slopes_err.append(slope_err)
            self.__by_fill_nls_slopes.append(nls_slope)
            self.__by_fill_nls_slopes_err.append(nls_slope_err)
            self.__by_fill_nls_slopes_redchi2.append(nls_chi2)
            self.__by_fill_nls_slopes_pvalue.append(p_value)
        else:
            print("** WARNING: Not enough points to fit: fill->" + str(fill) + ", data_size->" + str(len(x_nls)))

    def get_slope_average_per_intervals(self):
        data_to_use = self.__by_fill_df
        npoints = self.__average_points_in_summary
        number_of_fills = len(data_to_use)
        last_fill_index = number_of_fills - 1
        # delta_lumi = data_to_use[self.__label_col_accumulated_lumi][last_fill_index]/npoints
        number_of_fills_per_point = number_of_fills/npoints

        avg_slopes = []
        avg_err_slopes = []
        avg_mean_err_slopes = []
        avg_special_err_slopes = []
        fill_marker = []
        lumi_marker = []

        temp_slope_list = []
        temp_err_slope_list = []
        npoints_counter = 0

        data_to_use.reset_index(drop=True, inplace=True)

        for i in range(0, number_of_fills):
            npoints_counter += 1
            temp_slope_list.append(data_to_use[self.__label_col_nls_slopes][i])
            temp_err_slope_list.append(data_to_use[self.__label_col_nls_slopes_err][i])
            if npoints_counter >= number_of_fills_per_point or i >= number_of_fills - 1:
                stats = ltools.get_w_stats(temp_slope_list, 1/np.array(temp_err_slope_list))
                special_error_from_weighted_mean = ltools.get_w_mean_error(temp_slope_list,
                                                                           1/np.array(temp_err_slope_list),
                                                                           stats.mean)
                # print(special_error_from_weighted_mean)
                # print(stats.mean, stats.std_mean, stats.std, stats.sum_weights)
                avg_slopes.append(stats.mean)
                avg_mean_err_slopes.append(stats.std/math.sqrt(len(temp_slope_list)))
                avg_err_slopes.append(stats.std)
                avg_special_err_slopes.append(special_error_from_weighted_mean)
                fill_marker.append(data_to_use[self.__label_col_fill][i])
                lumi_marker.append(data_to_use[self.__label_col_accumulated_lumi][i])
                npoints_counter = 0
                temp_slope_list = []
                temp_err_slope_list = []

        self.__data_avg_summary[self.__label_col_nls_slopes] = np.array(avg_slopes)
        self.__data_avg_summary[self.__label_col_nls_slopes_err] = np.array(avg_err_slopes)
        self.__data_avg_summary["slope error in lumi interval (mean stdv)"] = np.array(avg_err_slopes)
        self.__data_avg_summary["slope error in lumi interval (mean error)"] = np.array(avg_mean_err_slopes)
        self.__data_avg_summary["slope error in lumi interval (mean weighted error)"] = np.array(avg_special_err_slopes)
        self.__data_avg_summary[self.__label_col_fill] = np.array(fill_marker)
        self.__data_avg_summary[self.__label_col_accumulated_lumi] = np.array(lumi_marker)

        slope_nls_vs_fill = plotting.plot_scatter_and_errors(data_frame=self.__data_avg_summary,
                                                             y_data_label=self.__label_col_nls_slopes,
                                                             x_data_label=self.__label_col_fill,
                                                             err_y_data_label=self.__label_col_nls_slopes_err,
                                                             ymin=self.__min_allowed_slope,
                                                             ymax=self.__max_allowed_slope,
                                                             xlabel="fill number",
                                                             ylabel="slope [hz/" + r'$\mu$' + "b]",
                                                             energy_year_label=self.__year_energy_label,
                                                             use_integers_in_x_axis=True)
        slope_nls_vs_lumi = plotting.plot_scatter_and_errors(data_frame=self.__data_avg_summary,
                                                             y_data_label=self.__label_col_nls_slopes,
                                                             x_data_label=self.__label_col_accumulated_lumi,
                                                             err_y_data_label=self.__label_col_nls_slopes_err,
                                                             ymin=self.__min_allowed_slope,
                                                             ymax=self.__max_allowed_slope,
                                                             xlabel="Integrated luminosity [$" +
                                                                    self.__lumi_unit + "^{-1}$]",
                                                             ylabel="slope [hz/" + r'$\mu$' + "b]",
                                                             energy_year_label=self.__year_energy_label)


        self.__plt_plots['avg_slope_nls_vs_fill'] = slope_nls_vs_fill
        self.__plt_plots['avg_slope_nls_vs_lumi'] = slope_nls_vs_lumi

        self.__data_avg_summary.to_csv(self.__output_dir + "avg_summary.csv", index=False)

    def plot_hist_sbil(self):
        sbil_hist = plotting.hist_from_pandas_frame(data_frame=self.__lin_data,
                                                    col_label=self.__sbil_label,
                                                    nbins=setts.nbins_sbil_histo,
                                                    xlabel=self.__ratios.det2.name +
                                                           " SBIL [hz/" + r'$\mu$' + "b]",
                                                    ylabel='Counts',
                                                    mean=self.__sbil_mean, stdv=self.__sbil_stdv,
                                                    xmin=setts.sbil_min, xmax=setts.sbil_max,
                                                    energy_year_label=self.__year_energy_label)
        self.__plt_plots['sbil_hist'] = sbil_hist[0][0].get_figure()

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
        self.__plt_plots['ratio_vs_all_data_sbil'] = ratio_vs_all_data_sbil.get_figure()

    def plot_by_fill_norm_ratio_vs_all_data_sbil(self):
        ratio_vs_all_data_sbil = plotting.scatter_plot_from_pandas_frame(data_frame=self.__lin_data,
                                                                         y_data_label=self.__label_norm_by_fill_ratio,
                                                                         x_data_label=self.__sbil_label,
                                                                         xlabel="SBIL [hz/" + r'$\mu$' + "b]",
                                                                         ylabel=self.__label_ratio + " ratios",
                                                                         ymin=setts.ratio_min, ymax=setts.ratio_max,
                                                                         energy_year_label=self.__year_energy_label,
                                                                         fig_size_shape='sq')
        self.__plt_plots['by_fill_norm_ratio_vs_all_data_sbil'] = ratio_vs_all_data_sbil.get_figure()

    def plot_by_nls_by_fill_norm_ratio_vs_all_data_sbil(self):
        ratio_vs_all_data_sbil = plotting.scatter_plot_from_pandas_frame(data_frame=self.__lin_data,
                                                                         y_data_label=self.__label_norm_by_nls_by_fill_ratio,
                                                                         x_data_label=self.__sbil_label,
                                                                         xlabel="SBIL [hz/" + r'$\mu$' + "b]",
                                                                         ylabel=self.__label_ratio + " ratios",
                                                                         ymin=setts.ratio_min, ymax=setts.ratio_max,
                                                                         energy_year_label=self.__year_energy_label,
                                                                         fig_size_shape='sq')
        self.__plt_plots['by_nls_by_fill_norm_ratio_vs_all_data_sbil'] = ratio_vs_all_data_sbil.get_figure()

    def save_plots(self):
        print('\n\n Saving linearity plots for ' + self.__label_ratio + ':')
        plotting.save_plots(self.__plt_plots, self.__output_dir)

    def save_summary_csv(self):
        summary_dict = {
            # slopes values from fill-by-fill
            'slope_hist': self.__by_fill_slopes_mean,
            'slope_hist_nls': self.__by_fill_nls_slopes_mean,
            'slope_hist_lw': self.__by_fill_slopes_mean_lw,
            'slope_hist_nls_lw': self.__by_fill_nls_slopes_mean_lw,
            'slope_hist_nls_errw': self.__by_fill_nls_slopes_mean_errw,
            'slope_hist_nls_totalw': self.__by_fill_nls_slopes_mean_totalw,
            'slope_stdv_hist': self.__by_fill_slopes_stdv,
            'slope_stdv_hist_nls': self.__by_fill_nls_slopes_stdv,
            'slope_stdv_hist_lw': self.__by_fill_slopes_stdv_lw,
            'slope_stdv_hist_nls_lw': self.__by_fill_nls_slopes_stdv_lw,
            'slope_stdv_hist_nls_errw': self.__by_fill_nls_slopes_stdv_errw,
            'slope_stdv_hist_nls_totalw': self.__by_fill_nls_slopes_stdv_totalw,
            'slope_err_hist_nls_totalw': self.__by_fill_nls_slopes_err_totalw,
            'slope_err_hist_nls_lw': self.__by_fill_nls_slopes_err_lw,
            # slope from all data fitting
            'all_data_fitted_slope': self.__all_data_fitted_slope,
            'all_data_fitted_slope_err': self.__all_data_fitted_slope_err,
            'all_data_fitted_slope_chi2': self.__all_data_fitted_slope_chi2
        }

        with open(self.__output_dir + "summary_values.json", 'w') as outfile:
            json.dump(summary_dict, outfile, indent=4)

    @property
    def sbil_mean(self):
        return self.__sbil_mean

    @property
    def sbil_stdv(self):
        return self.__sbil_stdv


class Stats:
    def __init__(self, array: pd.DataFrame) -> None:
        self.__mean = array.mean(numeric_only=True)
        self.__stdv = array.std(numeric_only=True)

    @property
    def mean(self):
        return self.__mean

    @property
    def stdv(self):
        return self.__stdv
