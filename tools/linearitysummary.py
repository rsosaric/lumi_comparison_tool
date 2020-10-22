import settings as setts
import pandas as pd
import numpy as np
import json
import tools.lumi_tools as ltools
from tools.linearityanalysis import LinearityAnalysis
import matplotlib.pyplot as plt
import tools.plotting_tools as plotting
import matplotlib.colors as mcolors


class LinearitySummary:
    def __init__(self, dets_file_labels: list, year: str, energy: str) -> None:
        self.__lumi_unit = 'fb'
        # Results for each loaded detector pair
        self.__results = []
        n_input_dets = len(dets_file_labels)
        det2_name = dets_file_labels[n_input_dets - 1]

        self.__analyis_label = dets_file_labels[0].upper()
        for i in range(1, n_input_dets):
            self.__analyis_label += ("-" + dets_file_labels[i].upper())

        self.__base_dir = setts.default_output_dir + year + "_" + energy + "/"

        for i in range(0, n_input_dets - 1):
            temp_lin_results = LinearityResults(detector_pair=[dets_file_labels[i], det2_name],
                                                base_dir=self.__base_dir, year=year, energy=energy)
            temp_lin_results.load_results()
            self.__results.append(temp_lin_results)

        self.__year = year
        self.__energy = energy

        self.__output_dir = setts.default_output_dir + year + "_" + energy + "/linearity_summary/" \
                            + self.__analyis_label + "/"
        # print(self.__output_dir)

        try:
            temp_energy = self.__energy
            if "TeV" in temp_energy:
                temp_energy = temp_energy.replace("TeV", "")
            self.__min_allowed_slope = setts.limits_linearity_slopes_year[int(self.__year), int(temp_energy)][0]
            self.__max_allowed_slope = setts.limits_linearity_slopes_year[int(self.__year), int(temp_energy)][1]
            print("Slopes range have been set to: (" + str(self.__min_allowed_slope) + "," +
                  str(self.__max_allowed_slope) + ")")
        except:
            print("no (min,max) slopes configuration found for " + str(self.__year) + " [" + str(
                self.__energy) + "TeV] in the settings.py")
            print("Setting no limits for the slopes!")
            self.__min_allowed_slope = None
            self.__max_allowed_slope = None

        self.__energy_year_label = self.__energy + "(" + self.__year + ")"

        self.__plt_plots = {}
        self.plot_avg_for_all_pairs()
        self.plot_avg_for_all_pairs(include_by_fill_slopes_mean=True)
        self.save_plots()

    def plot_avg_for_all_pairs(self, include_overall_fit=False, include_by_fill_slopes_mean=False):
        plot_name = "avg_slopes_for_all_pairs"
        fig_size_shape = "sq"
        fig_size = (8, 8)
        fig, ax = plt.subplots(figsize=fig_size)
        markersize = 5
        # ymin = self.__min_allowed_slope
        # ymax = self.__max_allowed_slope
        ymin = -0.01
        ymax = 0.01
        xlabel = "Integrated luminosity [$" + self.__lumi_unit + "^{-1}$]"
        ylabel = "slope [hz/" + r'$\mu$' + "b]"
        legend_labels = []
        leg_text_s = setts.leg_vs_plots_text_s
        leg_marker_sc = 1.0
        legend_position = 1

        colors = list(mcolors.TABLEAU_COLORS)

        i_color = 0
        for result in self.__results:
            label_ratio = result.label_ratio
            legend_labels.append(label_ratio + " avg. per lumi. interval")
            x = result.lumi_array_for_avg_slopes
            y = result.avg_slopes_array
            y_err = result.avg_slopes_err_array
            plt.errorbar(x, y, yerr=y_err, markersize=markersize, fmt='o', color=colors[i_color],
                         label=label_ratio + " avg. per lumi. interval")
            i_color += 1

        if include_by_fill_slopes_mean:
            i_color = 0
            for result in self.__results:
                plot_name += "_with_by_fill_mean"
                color = colors[i_color]
                label_ratio = result.label_ratio
                legend_labels.append(label_ratio + " per fill fit mean value")
                max_lumi = result.lumi_array_for_avg_slopes[len(result.lumi_array_for_avg_slopes) - 1]
                slope_value = result.slope_nls_totalw
                # slope_err = result.slope_nls_totalw_stdv
                slope_err = result.slope_nls_totalw_err
                n_x_points = 20
                delta_lumi = max_lumi/(n_x_points - 1)
                x = []
                y = []
                for i in range(0, n_x_points):
                    x.append(delta_lumi * i)
                    y.append(slope_value)
                x = np.array(x)
                y = np.array(y)
                plt.fill_between(x, y - slope_err, y + slope_err, label='_nolegend_', color=color, alpha=0.2)
                plt.plot(x, y, '--', label=label_ratio + " per fill fit mean value", color=color, alpha=0.8)
                i_color += 1

        # ax.legend(legend_labels, markerscale=leg_marker_sc, fontsize=leg_text_s, loc=legend_position)
        plt.legend(markerscale=leg_marker_sc, fontsize=leg_text_s, loc=legend_position)

        ax.set_xlabel(xlabel, labelpad=setts.axis_labelpad, weight=setts.axis_weight,
                      size=setts.axis_case_size[fig_size_shape])
        ax.set_ylabel(ylabel, labelpad=setts.axis_labelpad, weight=setts.axis_weight,
                      size=setts.axis_case_size[fig_size_shape])

        plt.xticks(fontsize=setts.axis_thicks_case_size[fig_size_shape])
        plt.yticks(fontsize=setts.axis_thicks_case_size[fig_size_shape])

        if ymin is not None and ymax is not None:
            plt.ylim(ymin, ymax)

        plt.subplots_adjust(left=0.18, right=0.97, top=0.95, bottom=0.1)

        plotting.add_extra_text(ax, fig_size_shape, energy_year_label=self.__energy_year_label,
                                experiment=setts.experiment, work_status=setts.work_status)


        self.__plt_plots[plot_name] = fig

    def save_plots(self):
        print('\n\n Saving linearity plots for ' + self.__analyis_label + ':')
        plotting.save_plots(self.__plt_plots, self.__output_dir)


class LinearityResults:
    def __init__(self, detector_pair: list, base_dir: str, year: str, energy: str) -> None:
        self.__detector_pair = detector_pair
        self.__name_det1 = detector_pair[0].upper()
        self.__name_det2 = detector_pair[1].upper()
        self.__label_ratio = self.__name_det1 + "/" + self.__name_det2
        self.__label = self.__name_det1 + "-" + self.__name_det2
        self.__input_dir = base_dir + self.__label + "/" + setts.linearity_analysis_output_dir
        self.__year = year
        self.__csv_input_dir = setts.csv_input_base_dir + year + '/'
        self.__energy = energy
        self.__path_to_avg_slopes_results_file = self.__input_dir + "avg_summary.csv"
        self.__path_to_results_file = self.__input_dir + "summary_values.json"

        # Variables initialization
        self.__avg_slopes_df = pd.DataFrame()
        self.__slopes_values_dict = {}
        self.__slope_nls_totalw = None
        self.__slope_nls_totalw_err = None
        self.__slope_overall_fit = None
        self.__slope_overall_fit_err = None

    def load_results(self):
        print("Loading linearity results from " + self.__name_det1 + " and " + self.__name_det2 + " from: " +
              self.__input_dir)
        if not ltools.check_folder_existence(self.__path_to_avg_slopes_results_file) \
                or not ltools.check_folder_existence(self.__path_to_results_file):
            print("     -->> Input files not found. Producing linearity results ...")
            LinearityAnalysis(dets_file_labels=self.__detector_pair, input_dir=self.__csv_input_dir)
            if not ltools.check_folder_existence(self.__path_to_avg_slopes_results_file) \
                    or not ltools.check_folder_existence(self.__path_to_results_file):
                ltools.color_print("There was something wrong in the production of the linearity analysis", "red")
                quit()

        # Loading avg summary results
        self.__avg_slopes_df = pd.read_csv(self.__path_to_avg_slopes_results_file, index_col=False, engine='python')
        # print(self.__avg_slopes_df)

        # Loading slope values results
        with open(self.__path_to_results_file) as json_file:
            self.__slopes_values_dict = json.load(json_file)
            # print(self.__slopes_values_dict)

        # Filling important variables
        self.__slope_nls_totalw = self.__slopes_values_dict["slope_hist_nls_totalw"]
        self.__slope_nls_totalw_stdv = self.__slopes_values_dict["slope_stdv_hist_nls_totalw"]
        self.__slope_nls_totalw_err = self.__slopes_values_dict["slope_err_hist_nls_totalw"]
        self.__slope_overall_fit = self.__slopes_values_dict["all_data_fitted_slope"]
        self.__slope_overall_fit_err = self.__slopes_values_dict["all_data_fitted_slope_err"]

        ltools.color_print("    Result summary:", "blue")
        print("         slope_nls_totalw: " + str(self.__slope_nls_totalw))
        print("         slope_overall_fit: " + str(self.__slope_overall_fit))

    @property
    def avg_slopes_df(self):
        return self.__avg_slopes_df

    @property
    def avg_slopes_array(self):
        return self.__avg_slopes_df["by nls slopes"]

    @property
    def avg_slopes_err_array(self):
        return self.__avg_slopes_df["by nls err on the slope"]

    @property
    def fill_array_for_avg_slopes(self):
        return self.__avg_slopes_df["fill"]

    @property
    def lumi_array_for_avg_slopes(self):
        return self.__avg_slopes_df["accumulated lumi up to fill"]

    @property
    def slope_nls_totalw (self):
        return self.__slope_nls_totalw

    @property
    def slope_nls_totalw_stdv(self):
        return self.__slope_nls_totalw_stdv

    @property
    def slope_nls_totalw_err(self):
        return self.__slope_nls_totalw_err

    @property
    def slope_overall_fit(self):
        return self.__slope_overall_fit

    @property
    def slope_overall_fit_err(self):
        return self.__slope_overall_fit_err

    @property
    def label_ratio(self):
        return self.__label_ratio
