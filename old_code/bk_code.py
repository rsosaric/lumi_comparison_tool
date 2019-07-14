def create_all_nls_data_sim_method(self):
    if self.__all_nls_data is not None:
        raise Warning('Overwriting all_nls_data with simultaneous method!')

    for ratio_i in self.__ratios:
        ratio_i.create_nls_data()
        self.__nls_ratio_col_names.append(ratio_i.by_nls_label_ratio)
        self.__nls_ratio_plotting_labels.append(ratio_i.label_ratio)
    self.__ratio_col_names = self.__nls_ratio_plotting_labels

    nls = self.__nls
    print(nls)

    inls = 0
    sum_lumi1 = 0.0
    sum_lumi2 = 0.0
    sum_lumi3 = 0.0
    nls_ratios_array13 = []
    nls_ratios_array23 = []
    nls_ratios_array12 = []
    nls_data_array = []

    df_to_use = self.all_data

    for index_data in range(0, len(df_to_use)):
        inls += 1
        sum_lumi1 += df_to_use[self.__dets[0].lumi_rec_label][index_data]
        sum_lumi2 += df_to_use[self.__dets[1].lumi_rec_label][index_data]
        sum_lumi3 += df_to_use[self.__dets[2].lumi_rec_label][index_data]

        nls_ratios_array13.append(df_to_use[self.__ratios13.label_ratio][index_data])
        nls_ratios_array23.append(df_to_use[self.__ratios23.label_ratio][index_data])
        nls_ratios_array12.append(df_to_use[self.__ratios12.label_ratio][index_data])

        if inls == nls or index_data == len(df_to_use) - 1:

            if sum_lumi2 > 0.0:
                by_nls_ratio12 = sum_lumi1 / sum_lumi2
                by_nls_ratios_stdv12 = np.array(nls_ratios_array12).std()
                by_nls_ratios_mean_err12 = by_nls_ratios_stdv12 / np.sqrt(inls)
            else:
                by_nls_ratio12 = np.nan
                by_nls_ratios_stdv12 = np.nan
                by_nls_ratios_mean_err12 = np.nan

            if sum_lumi3 > 0.0:
                by_nls_ratio13 = sum_lumi1 / sum_lumi3
                by_nls_ratio23 = sum_lumi2 / sum_lumi3
                by_nls_ratios_stdv13 = np.array(nls_ratios_array13).std()
                by_nls_ratios_stdv23 = np.array(nls_ratios_array23).std()
                by_nls_ratios_mean_err13 = by_nls_ratios_stdv13 / np.sqrt(inls)
                by_nls_ratios_mean_err23 = by_nls_ratios_stdv23 / np.sqrt(inls)
            else:
                by_nls_ratio13 = np.nan
                by_nls_ratio23 = np.nan
                by_nls_ratios_stdv13 = np.nan
                by_nls_ratios_stdv23 = np.nan
                by_nls_ratios_mean_err13 = np.nan
                by_nls_ratios_mean_err23 = np.nan

            nls_data_array.append(
                (df_to_use['fill'][index_data], df_to_use['run'][index_data],
                 df_to_use['ls'][index_data], df_to_use['time'][index_data],
                 df_to_use[self.__ratios13.accumulated_rec_lumi1_label][index_data],
                 df_to_use[self.__ratios13.accumulated_rec_lumi2_label][index_data],
                 df_to_use[self.__ratios12.accumulated_rec_lumi2_label][index_data],
                 by_nls_ratio13, by_nls_ratios_mean_err13, by_nls_ratio23, by_nls_ratios_mean_err23,
                 by_nls_ratio12, by_nls_ratios_mean_err12
                 ))

            inls = 0
            sum_lumi1 = 0.0
            sum_lumi2 = 0.0
            sum_lumi3 = 0.0
            nls_ratios_array12 = []
            nls_ratios_array13 = []
            nls_ratios_array23 = []

    nls_data_col_types = (
        [('fill', np.int64), ('run', np.int64), ('ls', np.int64), ('time', np.int64),
         (self.__ratios13.accumulated_rec_lumi1_label, np.float64),
         (self.__ratios13.accumulated_rec_lumi2_label, np.float64),
         (self.__ratios12.accumulated_rec_lumi2_label, np.float64),
         (self.__ratios13.by_nls_label_ratio, np.float64),
         ('mean_err_' + self.__ratios13.by_nls_label_ratio, np.float64),
         (self.__ratios23.by_nls_label_ratio, np.float64),
         ('mean_err_' + self.__ratios23.by_nls_label_ratio, np.float64),
         (self.__ratios12.by_nls_label_ratio, np.float64),
         ('mean_err_' + self.__ratios12.by_nls_label_ratio, np.float64)]
    )

    nls_data_array_np = np.array(nls_data_array, dtype=nls_data_col_types)
    self.__all_nls_data = pd.DataFrame(nls_data_array_np)
    ltools.add_date_column(self.__all_nls_data)

    # self.__nls_ratios_mean = self.__nls_data[self.by_nls_label_ratio].mean()
    # self.__nls_ratios_stdv = self.__nls_data[self.by_nls_label_ratio].std()


# def plot_nls_ratios_vs_lumi3(self):
    #     nls_ratios_vs_lumi3 = plotting.scatter_plot_from_pandas_frame(data_frame=self.all_nls_data,
    #                                                                   y_data_label=self.__nls_ratio_col_names,
    #                                                                   x_data_label=self.__lumi3_col_name,
    #                                                                   ylabel="ratios in " + str(self.__nls) + ' LS',
    #                                                                   ymin=setts.ratio_min, ymax=setts.ratio_max,
    #                                                                   energy_year_label=self.__year_energy_label,
    #                                                                   legend_labels=self.__nls_ratio_plotting_labels)
    #     self.__plt_plots.append(['nls_ratios_vs_lumi3', nls_ratios_vs_lumi3.get_figure()])


# def plot_nls_ratio_hist(self):
    #     ratio_hist = plotting.hist_from_pandas_frame(data_frame=self.nls_data, col_label=self.by_nls_label_ratio,
    #                                                  nbins=setts.nbins,
    #                                                  xlabel=self.__label_ratio + " ratios in " + str(
    #                                                      self.__nls) + ' LS',
    #                                                  ylabel='Counts',
    #                                                  # title='by Nls Detectors Ratios Histogram',
    #                                                  xmin=setts.ratio_min, xmax=setts.ratio_max,
    #                                                  mean=self.__nls_ratios_mean, stdv=self.__nls_ratios_stdv,
    #                                                  energy_year_label=self.__year_energy_label)
    #     self.__plt_plots.append(['ratio_nls_hist', ratio_hist[0][0].get_figure()])