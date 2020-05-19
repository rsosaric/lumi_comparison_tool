from tools.detectorsratio import DetectorsRatio as Ratios
from tools.luminometer import Luminometer as L
import tools.plotting_tools as plotting
import tools.lumi_tools as ltools
import settings as setts
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from lmfit import Model
import matplotlib.pyplot as plt


def lin_func(x, a, b):
    return a * x + b


lin_model = Model(lin_func)


class LinearityAnalysis:
    __conversion_fb_to_ub = 1000000000.0 / setts.xLS
    __linearity_lumi_conversion = __conversion_fb_to_ub

    def __init__(self, dets_file_labels: list, input_dir: str, lumi_type: str) -> None:

        self.__label_col_fill = "fill"
        self.__label_col_slopes = "slopes"
        self.__label_col_slopes_err = "err on the slope"
        self.__label_col_slopes_chi2 = "chi2 on the slope"
        self.__label_col_nls_slopes = "by nls slopes"
        self.__label_col_nls_slopes_err = "by nls err on the slope"
        self.__label_col_nls_slopes_chi2 = "by nls chi2 on the slope"
        self.__label_col_fill_lumi = "lumi in fill"
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
        self.__by_fill_stats = {}
        self.__by_fill_slopes = []
        self.__by_fill_slopes_err = []
        self.__by_fill_nls_slopes = []
        self.__by_fill_nls_slopes_err = []
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
        self.__excluded_fills_in_running = []
        self.__good_fills = []

        self.__by_run_data = {}
        self.__by_run_mean = {}
        self.__by_run_std = {}

        self.create_by_run_data()
        self.create_by_fill_data()
        self.get_normalize_by_fill_data()
        self.fit_all_data()

        # print(self.__by_fill_df)

        self.plot_hist_sbil()
        self.plot_ratio_vs_all_data_sbil()
        self.plot_by_fill_norm_ratio_vs_all_data_sbil()
        self.plot_by_nls_by_fill_norm_ratio_vs_all_data_sbil()
        self.save_plots()

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

            # Fitting
            if not exclude_fill:
                fitted_data_ini = len(self.__by_fill_slopes)
                self.fill_fitting(fill)
                fitted_data_changed = len(self.__by_fill_slopes)

                if fitted_data_changed == fitted_data_ini + 1:
                    good_fills.append(fill)
                    self.__by_fill_lumi.append(np.sum(self.__by_fill_data[fill][self.__label_lumi]))
                else:
                    print("     -->> Fill " + str(fill) + " not included in final results.")
            else:
                print("     -->> Fill " + str(fill) + " not included in final results.")

        if len(good_fills) == len(self.__by_fill_lumi) and len(good_fills) == len(self.__by_fill_slopes):
            self.__good_fills = good_fills
            self.__by_fill_df[self.__label_col_fill] = good_fills
            self.__by_fill_df[self.__label_col_fill_lumi] = self.__by_fill_lumi
            self.__by_fill_df[self.__label_col_slopes] = self.__by_fill_slopes
            self.__by_fill_df[self.__label_col_slopes_err] = self.__by_fill_slopes_err
            self.__by_fill_df[self.__label_col_nls_slopes] = self.__by_fill_nls_slopes
            self.__by_fill_df[self.__label_col_nls_slopes_err] = self.__by_fill_nls_slopes_err

            self.__by_fill_df.replace([np.inf, -np.inf], np.nan, inplace=True)
            self.__by_fill_df.dropna(inplace=True)

            self.__by_fill_df[self.__label_col_total_w] = ltools.get_total_weights(
                self.__by_fill_df[self.__label_col_fill_lumi],
                self.__by_fill_df[self.__label_col_nls_slopes_err])
        else:
            raise AssertionError("ERROR in create_by_fill_data(): lengths of good_fills, by_fill_slopes and "
                                 "by_fill_lumi do not match! ")

        self.__by_fill_df_all = self.__by_fill_df

        # drop slopes outside the range in settings
        all_slopes_n = len(self.__by_fill_df)
        if self.__min_allowed_slope is not None and self.__max_allowed_slope is not None:
            self.__by_fill_df = self.__by_fill_df[
                (self.__by_fill_df[self.__label_col_slopes] >= self.__min_allowed_slope) &
                (self.__by_fill_df[self.__label_col_slopes] <= self.__max_allowed_slope) &
                (self.__by_fill_df[self.__label_col_nls_slopes] >= self.__min_allowed_slope) &
                (self.__by_fill_df[self.__label_col_nls_slopes] <= self.__max_allowed_slope)]
            after_excl_slopes_n = len(self.__by_fill_df)
            self.__by_fill_df.reset_index(drop=True)
            print(str((all_slopes_n - after_excl_slopes_n) * 100 / all_slopes_n) +
                  " % of the fills have been excluded because they had slopes ouside the allowed range.")
            print("PLEASE CHECK IF THIS NUMBER IS BIG!")

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
                                                                ylabel="Integrated luminosity [$" +
                                                                       self.__ratios.lumi_unit + "^{-1}$]",
                                                                xmin=self.__min_allowed_slope,
                                                                xmax=self.__max_allowed_slope,
                                                                mean=self.__by_fill_nls_slopes_mean_totalw,
                                                                stdv=self.__by_fill_nls_slopes_stdv_totalw,
                                                                energy_year_label=self.__year_energy_label,
                                                                weight_label=self.__label_col_total_w,
                                                                mean_float_format="{0:.4f}")

        self.__plt_plots['slope_hist'] = slope_hist
        self.__plt_plots['slope_hist_nls'] = slope_hist_nls

        self.__plt_plots['slope_hist_lw'] = slope_hist_lw[0][0].get_figure()
        self.__plt_plots['slope_hist_nls_lw'] = slope_hist_nls_lw[0][0].get_figure()
        self.__plt_plots['slope_hist_nls_errw'] = slope_hist_nls_errw[0][0].get_figure()
        self.__plt_plots['slope_hist_nls_totalw'] = slope_hist_nls_totalw[0][0].get_figure()

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

    def fit_all_data(self) -> None:
        if len(self.__lin_data[self.__label_norm_by_nls_by_fill_ratio]) == 0:
            raise AssertionError("norm_by_nls_by_fill_ratio data is needed")

        to_fit_data_nls = pd.DataFrame()

        to_fit_data_nls["x_nls"] = self.__lin_data[self.__nls_sbil_label]
        to_fit_data_nls["y_nls"] = self.__lin_data[self.__nls_ratio_label]
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

        nls_fig = plotting.plot_from_fit(x, y, y_err=ey,
                                         fitted_slope=slope, chi2=chi2,
                                         fitted_slope_err=slope_err,
                                         fitted_intercept=intercept,
                                         fitted_f=lin_func,
                                         xlabel="SBIL [hz/" + r'$\mu$' + "b]",
                                         ylabel=self.__label_ratio + " ratios in " +
                                                str(self.__nls) + ' LS',
                                         energy_year_label=self.__year_energy_label)

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
            nls_chi2 = result_nls_fit.redchi

            fig = plotting.plot_from_fit(x, y, fitted_slope=slope, fitted_intercept=intercept,
                                         fitted_slope_err=slope_err, chi2=chi2,
                                         fitted_f=lin_func, xlabel="SBIL [hz/" + r'$\mu$' + "b]",
                                         ylabel=self.__label_ratio + " ratio",
                                         energy_year_label=self.__year_energy_label)

            nls_fig = plotting.plot_from_fit(x_nls, y_nls, y_err=ey_nls,
                                             fitted_slope=nls_slope, chi2=nls_chi2,
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
            self.__by_fill_nls_slopes_err.append(nls_slope_err)
        else:
            print("** WARNING: Not enough points to fit: fill->" + str(fill) + ", data_size->" + str(len(x_nls)))

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
