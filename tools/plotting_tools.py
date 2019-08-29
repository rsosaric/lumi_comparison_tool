import matplotlib.pyplot as plt
plt.rcParams.update({'figure.max_open_warning': 0})
import settings as setts
import tools.lumi_tools as ltools
from statsmodels.graphics.api import abline_plot
import seaborn as sns
import numpy as np
import pandas as pd

__cms_label_pos_sq = setts.cms_label_pos_sq
__cms_label_pos_nsq = setts.cms_label_pos_nsq
__delta_y_pos_sq = setts.delta_y_pos_sq
__delta_y_pos_nsq = setts.delta_y_pos_nsq

__year_energy_label_pos_sq = setts.year_energy_label_pos_sq
__year_energy_label_pos_nsq = setts.year_energy_label_pos_nsq


def save_py_fig_to_file(fig, output_name):
    fig.savefig(output_name)
    plt.close(fig)


def save_plots(names_and_plots, output_name):
    ltools.check_and_create_folder(output_name)
    plot_names = list(names_and_plots)
    for plot_name in plot_names:
        plot_object = names_and_plots[plot_name]
        for plot_format in setts.plots_formats:
            plot_full_path = output_name + plot_name + '.' + plot_format
            save_py_fig_to_file(plot_object, plot_full_path)
            print('saved plot: ' + plot_full_path)


# TODO: set histo range for same binning
def hist_from_pandas_frame(data_frame, col_label, nbins, title='', xlabel='', ylabel='', xmin=0.0, xmax=0.0,
                           x_data_range=None,
                           weight_label=None,
                           mean=None, stdv=None, err_mean=None, color='#86bf91',
                           label_cms_status=True, energy_year_label='',
                           fig_size_shape='sq'):
    fig_size = get_fig_size(fig_size_shape)
    if x_data_range is None:
        x_data_range = (xmin, xmax)
    if weight_label:
        ratio_hist = data_frame.hist(bins=nbins, column=col_label, grid=False, weights=data_frame[weight_label],
                                     figsize=fig_size, sharex=True, color=color,
                                     range=x_data_range)
    else:
        ratio_hist = data_frame.hist(bins=nbins, column=col_label, grid=False,
                                     figsize=fig_size, sharex=True, color=color,
                                     range=x_data_range)

    ratio_hist_ax = ratio_hist[0][0]

    ratio_hist_ax.set_title(title)
    ratio_hist_ax.set_xlabel(xlabel, labelpad=setts.axis_labelpad, weight=setts.axis_weight, size=setts.axis_case_size[fig_size_shape])
    ratio_hist_ax.set_ylabel(ylabel, labelpad=setts.axis_labelpad, weight=setts.axis_weight, size=setts.axis_case_size[fig_size_shape])

    plt.xticks(fontsize=setts.axis_thicks_case_size[fig_size_shape])
    plt.yticks(fontsize=setts.axis_thicks_case_size[fig_size_shape])

    if xmin != 0 and xmax != 0:
        plt.xlim(xmin, xmax)

    # histogram info
    if stdv:
        if mean:
            plt.text(setts.pos_leg_1[0], setts.pos_leg_1[1], 'mean: ' + str(float("{0:.3f}".format(mean))), ha='left',
                     fontsize=setts.leg_font_size, fontweight='light', transform=ratio_hist_ax.transAxes)
        plt.text(setts.pos_leg_2[0], setts.pos_leg_2[1], r'$\sigma$: ' + str(float("{0:.4f}".format(stdv))), ha='left',
                 fontsize=setts.leg_font_size, fontweight='light', transform=ratio_hist_ax.transAxes)
        if err_mean:
            plt.text(setts.pos_leg_3[0], setts.pos_leg_3[1], 'err_mean: ' + str(float("{0:.4f}".format(err_mean))),
                     ha='left',
                     fontsize=setts.leg_font_size, fontweight='light', transform=ratio_hist_ax.transAxes)
    elif mean:
        plt.text(setts.pos_leg_1[0], setts.pos_leg_1[1], 'mean: ' + str(float("{0:.3f}".format(mean))), ha='left',
                 fontsize=setts.leg_font_size,
                 fontweight='light', transform=ratio_hist_ax.transAxes)
    elif err_mean:
        plt.text(setts.pos_leg_1[0], setts.pos_leg_1[1], 'err_mean: ' + str(float("{0:.4f}".format(err_mean))),
                 ha='left',
                 fontsize=setts.leg_font_size,
                 fontweight='light', transform=ratio_hist_ax.transAxes)

    if label_cms_status:
        add_extra_text(ratio_hist_ax, fig_size_shape, energy_year_label=energy_year_label,
                       experiment=setts.experiment, work_status=setts.work_status)

    return ratio_hist


def hist_from_array(data, nbins, title='', xlabel='', ylabel='', xmin=0.0, xmax=0.0,
                    weight=None,
                    mean=None, stdv=None, err_mean=None, color='#86bf91',
                    label_cms_status=True, energy_year_label='',
                    fig_size_shape='sq'):
    fig_size = get_fig_size(fig_size_shape)
    fig, ax = plt.subplots(figsize=fig_size)
    if weight:
        plt.hist(data, bins=nbins, weights=weight,
                 color=color, zorder=2, rwidth=0.9)
    else:
        plt.hist(data, bins=nbins,
                 color=color, zorder=2, rwidth=0.9)

    ax.set_title(title)
    ax.set_xlabel(xlabel, labelpad=setts.axis_labelpad, weight=setts.axis_weight, size=setts.axis_case_size[fig_size_shape])
    ax.set_ylabel(ylabel, labelpad=setts.axis_labelpad, weight=setts.axis_weight, size=setts.axis_case_size[fig_size_shape])
    plt.xticks(fontsize=setts.axis_thicks_case_size[fig_size_shape])
    plt.yticks(fontsize=setts.axis_thicks_case_size[fig_size_shape])

    if xmin != 0 and xmax != 0:
        plt.xlim(xmin, xmax)

    # histogram info
    if stdv:
        if mean:
            plt.text(setts.pos_leg_1[0], setts.pos_leg_1[1], 'mean: ' + str(float("{0:.3f}".format(mean))), ha='left',
                     fontsize=setts.leg_font_size, fontweight='light', transform=ax.transAxes)
        plt.text(setts.pos_leg_2[0], setts.pos_leg_2[1], r'$\sigma$: ' + str(float("{0:.4f}".format(stdv))), ha='left',
                 fontsize=setts.leg_font_size, fontweight='light', transform=ax.transAxes)
        if err_mean:
            plt.text(setts.pos_leg_3[0], setts.pos_leg_3[1], 'err_mean: ' + str(float("{0:.4f}".format(err_mean))),
                     ha='left',
                     fontsize=setts.leg_font_size, fontweight='light', transform=ax.transAxes)
    elif mean:
        plt.text(setts.pos_leg_1[0], setts.pos_leg_1[1], 'mean: ' + str(float("{0:.3f}".format(mean))), ha='left',
                 fontsize=setts.leg_font_size,
                 fontweight='light', transform=ax.transAxes)
    elif err_mean:
        plt.text(setts.pos_leg_1[0], setts.pos_leg_1[1], 'err_mean: ' + str(float("{0:.4f}".format(err_mean))),
                 ha='left',
                 fontsize=setts.leg_font_size,
                 fontweight='light', transform=ax.transAxes)

    if label_cms_status:
        add_extra_text(ax, fig_size_shape, energy_year_label=energy_year_label,
                       experiment=setts.experiment, work_status=setts.work_status)

    return fig


def scatter_plot_from_pandas_frame(data_frame, x_data_label, y_data_label, title='', xlabel='', ylabel='',
                                   ymin=0.0, ymax=0.0, color='#86bf91', label_cms_status=True,
                                   energy_year_label='', legend_labels=None, legend_position=setts.leg_vs_plots_pos,
                                   leg_marker_sc=setts.leg_vs_plots_marker_scale,
                                   marker_size=0.7,
                                   leg_text_s=setts.leg_vs_plots_text_s,
                                   plot_style='o',
                                   fig_size_shape='nsq'):
    fig_size = get_fig_size(fig_size_shape)
    if xlabel == '':
        xlabel = x_data_label
    if ylabel == '':
        ylabel = y_data_label
    fig, ax = plt.subplots()
    if type(y_data_label) == list:
        plot = data_frame.plot(x=x_data_label, y=y_data_label, style='o', figsize=fig_size,
                               markersize=marker_size, ax=ax)
        if legend_labels is None:
            legend_labels = y_data_label
        else:
            assert len(legend_labels) == len(y_data_label)
        ax.legend(legend_labels, ncol=len(legend_labels),
                  markerscale=leg_marker_sc, fontsize=leg_text_s, loc=legend_position)
    else:
        plot = data_frame.plot(x=x_data_label, y=y_data_label, style=plot_style, figsize=fig_size,
                               markersize=marker_size, legend=None, ax=ax)
    plot.set_title(title)
    plot.set_ylabel(ylabel, labelpad=setts.axis_labelpad, weight=setts.axis_weight, size=setts.axis_case_size[fig_size_shape])
    plot.set_xlabel(xlabel, labelpad=setts.axis_labelpad, weight=setts.axis_weight, size=setts.axis_case_size[fig_size_shape])

    plt.xticks(fontsize=setts.axis_thicks_case_size[fig_size_shape])
    plt.yticks(fontsize=setts.axis_thicks_case_size[fig_size_shape])

    if label_cms_status:
        add_extra_text(ax, fig_size_shape, energy_year_label=energy_year_label,
                       experiment=setts.experiment, work_status=setts.work_status)

    # margin optimization
    if fig_size == (12, 4):
        plt.subplots_adjust(left=0.1, right=0.97, top=0.9, bottom=0.2)
    # else:
    #     raise Warning('nsq_fig_size not optimized')

    if ymin != 0 and ymax != 0:
        plt.ylim(ymin, ymax)

    return plot


# def plot_from_fit(data_frame, model_fit, x_data_label, y_data_label, fig_size='sq'):
#     fig, ax = plt.subplots()
#     data_plot = data_frame.plot(x=x_data_label, y=y_data_label, style='o', figsize=fig_size,
#                                 markersize=0.5, ax=ax)
#     abline_plot(model_results=model_fit, ax=ax)
#
#     return data_plot
def plot_from_fit(x, y, fitted_slope, fitted_slope_err, fitted_intercept, fitted_f,
                  energy_year_label,
                  xlabel='', ylabel='', y_err=[]):
    fig_size_shape = 'sq'
    fig, ax = plt.subplots(figsize=(8, 8))
    if len(y_err) > 0:
        assert len(y_err) == len(x)
        plt.errorbar(x, y, yerr=y_err, fmt="none")
    else:
        plt.scatter(x, y)
    xfine = np.linspace(np.min(x), np.max(x), 100)
    plt.plot(xfine, fitted_f(xfine, a=fitted_slope, b=fitted_intercept), 'r-')

    ax.set_xlabel(xlabel, labelpad=setts.axis_labelpad, weight=setts.axis_weight, size=setts.axis_case_size[fig_size_shape])
    ax.set_ylabel(ylabel, labelpad=setts.axis_labelpad, weight=setts.axis_weight, size=setts.axis_case_size[fig_size_shape])

    plt.xticks(fontsize=setts.axis_thicks_case_size[fig_size_shape])
    plt.yticks(fontsize=setts.axis_thicks_case_size[fig_size_shape])

    add_extra_text(ax, fig_size_shape, energy_year_label=energy_year_label,
                   experiment=setts.experiment, work_status=setts.work_status)

    # Adding fitting info
    plt.text(setts.pos_by_fill_fit_info[0], setts.pos_by_fill_fit_info[1],
             r"slope[hz/" + r"$\mu$" + "b]: " + str(float("{0:.4f}".format(fitted_slope))) +
             r' $\pm$ ' + str(float("{0:.4f}".format(fitted_slope_err))),
             ha='left',
             fontsize=setts.leg_font_size, fontweight='bold', transform=ax.transAxes)

    return fig


def hist_list_from_pandas_frame(list_nls_data, nbins, title='', xlabel='', ylabel='', xmin=0.0, xmax=0.0,
                                x_data_range=None,
                                all_nls_weights=pd.DataFrame(),
                                label_cms_status=True, energy_year_label='',
                                fig_size_shape='sq',
                                legend_labels=None, legend_position=setts.leg_vs_plots_pos,
                                leg_marker_sc=setts.leg_vs_plots_marker_scale,
                                leg_text_s=setts.leg_vs_plots_text_s, nls_label=True,
                                histtype='step', stacked=True, fill=False
                                ):
    global ratio_hist
    fig_size = get_fig_size(fig_size_shape)
    fig, ax = plt.subplots(figsize=fig_size)
    legend_labels_list = []

    if x_data_range is None:
        x_data_range = (xmin, xmax)

    for hist_name in list(list_nls_data):
        if all_nls_weights.empty:
            list_nls_data[hist_name].hist(bins=nbins, alpha=1, ax=ax, label=hist_name, grid=False, figsize=fig_size,
                                            range=x_data_range, histtype=histtype, stacked=stacked, fill=fill)

        else:
            list_nls_data[hist_name].dropna().hist(bins=nbins, alpha=1, ax=ax, label=hist_name, grid=False, figsize=fig_size,
                                                    range=x_data_range, histtype=histtype, stacked=stacked, fill=fill,
                                                    weights=all_nls_weights[hist_name].dropna())

    if legend_labels is None:
        print("No legend set")
    else:
        if nls_label:
            for nls in legend_labels:
                legend_labels_list.append(str(nls) + " LS")
        else:
            legend_labels_list = legend_labels
        ax.legend(legend_labels_list, ncol=1,
                    markerscale=leg_marker_sc, fontsize=leg_text_s, loc=legend_position)

    ax.set_title(title)
    ax.set_xlabel(xlabel, labelpad=setts.axis_labelpad, weight=setts.axis_weight, size=setts.axis_case_size[fig_size_shape])
    ax.set_ylabel(ylabel, labelpad=setts.axis_labelpad, weight=setts.axis_weight, size=setts.axis_case_size[fig_size_shape])

    plt.xticks(fontsize=setts.axis_thicks_case_size[fig_size_shape])
    plt.yticks(fontsize=setts.axis_thicks_case_size[fig_size_shape])

    if xmin != 0 and xmax != 0:
        plt.xlim(xmin, xmax)

    if label_cms_status:
        add_extra_text(ax, fig_size_shape, energy_year_label=energy_year_label,
                       experiment=setts.experiment, work_status=setts.work_status)

    return fig


def plot_bad_fill_info(data_frame, x_data_label, y_data_label, z_data_label, title='', xlabel='', ylabel='',
                                   ymin=0.0, ymax=0.0, label_cms_status=True,
                                   energy_year_label='', mean=None, stdv=None,
                                   fig_size_shape='nsq', ratio_acceptance=0.01, filePath='', txtfileName="NoName"):
    fig_size = get_fig_size(fig_size_shape)
    fig, ax = plt.subplots(figsize=fig_size)

    x_data_list = list(data_frame[x_data_label])
    y_data_list = list(data_frame[y_data_label])
    z_data_list = list(data_frame[z_data_label])

    data = pd.DataFrame()

    for i in range(0,len(x_data_list), 1):
        if str(x_data_list[i]) in data.keys():
            data[str(x_data_list[i])][0] = data[str(x_data_list[i])][0] + y_data_list[i]*z_data_list[i]
            data[str(x_data_list[i])][1] = data[str(x_data_list[i])][1] + z_data_list[i]
        else:
            data[str(x_data_list[i])] = (y_data_list[i]*z_data_list[i], z_data_list[i])

    total_lumi = 0
    for i in list(data):
        data[i][0] = data[i][0]/data[i][1]
        total_lumi = total_lumi + data[i][1]

    x_data_reduced = []
    y_data_mean = []
    z_data_acum = []
    ey = []
    colors = []

    for i in list(data):
        x_data_reduced.append(int(i))
        y_data_mean.append(data[i][0])
        z_data_acum.append(data[i][1])
        ey.append(1000*data[i][1]/total_lumi)
        colors.append(100*data[i][1]/total_lumi)

    lumi_porcent = max(colors)*setts.lumisensitivity
    limratioUP = mean + stdv*ratio_acceptance
    limratioDOWN = mean - stdv*ratio_acceptance

    bad_fills = []
    bad_fillsPos = []

    for i in range(0,len(x_data_reduced), 1):
        if (y_data_mean[i] > limratioUP or y_data_mean[i] < limratioDOWN) and (100*z_data_acum[i]/total_lumi) > lumi_porcent:
            bad_fills.append(x_data_reduced[i])
            bad_fillsPos.append((x_data_reduced[i],y_data_mean[i]))

    print(x_data_label + "s to analize", bad_fills)
    ##write file:
    ltools.check_and_create_folder(filePath)
    fileout = open(filePath + txtfileName + ".txt", "w+")
    for i in bad_fills:
        fileout.write(str(i) + ' ')
    fileout.close()

    sc = plt.scatter(x_data_reduced, y_data_mean, s=ey, c=colors, vmin=0, vmax=max(colors))

    plt.plot(x_data_reduced, y_data_mean, linewidth=0.5, alpha=0.8, linestyle="--")
    plt.axhline(mean, color='black', lw=0.8, alpha=0.7, linestyle="--")
    plt.axhline(limratioUP, color='red', lw=0.8, alpha=0.7, linestyle="--")
    plt.axhline(limratioDOWN, color='red', lw=0.8, alpha=0.7, linestyle="--")

    for i, txt in enumerate(bad_fills):
        plt.annotate(txt, bad_fillsPos[i], xytext=(-25, 25),
                     textcoords='offset points', ha='right', va='top',
                     bbox=dict(boxstyle='round,pad=0.3', fc='yellow', alpha=0.5),
                     arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))
    plt.colormaps()
    cbar = plt.colorbar(sc)
    cbar.set_label('% of the total integrated luminosity')

    # plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    if ymin != 0 and ymax != 0:
        plt.ylim(ymin, ymax)

    if label_cms_status:
        add_extra_text(ax, fig_size_shape, energy_year_label=energy_year_label,
                       experiment=setts.experiment, work_status=setts.work_status)

    return fig


## All/excluded plots
def snsplot_detector_all_and_excluded(data_frame, x_data_label, y_data_label, conditional_label, title='', xlabel='', ylabel='',
                                      ymin=None, ymax=None, xmin=None, xmax=None, xlabel_rotation = None,
                                      label_cms_status=True, energy_year_label='', fig_size_shape='nsq', use_pts_white_border=False,
                                      marker_size=5, leg_col=None):

    fig_size = get_fig_size(fig_size_shape)
    fig, ax = plt.subplots(figsize=fig_size)

    assert (type(y_data_label) == type(conditional_label))

    if type(y_data_label) == list and type(conditional_label) == list:
        for excl_label in conditional_label:
            if 'by nls exclusion (' not in excl_label:
                raise AssertionError("Only binary exclusion label [by_nls_binary_exclusion_info_label] "
                                     "have been implemented for this plot function. Input label: " + excl_label)

        colors = ('blue', 'darkorange', 'magneta')
        excluded_color = 'red'
        markers = ("o", "v", "s")
        marker_size = 13

        detcs_excld_data = {}
        detcs__data = {}

        for key in y_data_label:
            detcs_excld_data[key] = {"x": [], "y": []}
            detcs__data[key] = {"x": [], "y": []}

        for index_data in range(0, len(data_frame)):
            for detector_pair_index in range(0, len(y_data_label)):
                if data_frame[conditional_label[detector_pair_index]][index_data] == "included":
                    detcs__data[y_data_label[detector_pair_index]]["x"].append(data_frame[x_data_label][index_data])
                    detcs__data[y_data_label[detector_pair_index]]["y"].append(data_frame[y_data_label[detector_pair_index]][index_data])
                elif data_frame[conditional_label[detector_pair_index]][index_data] == "excluded":
                    detcs_excld_data[y_data_label[detector_pair_index]]["x"].append(data_frame[x_data_label][index_data])
                    detcs_excld_data[y_data_label[detector_pair_index]]["y"].append(data_frame[y_data_label[detector_pair_index]][index_data])
                else:
                    raise AssertionError("something wrong!")

        for detector_pair_index in range(0, len(y_data_label)-1):
            data_label = y_data_label[detector_pair_index]
            label_in_plot = data_label.split("_")[-1]
            sns.scatterplot(y=detcs__data[data_label]["y"], x=detcs__data[data_label]["x"], color=colors[detector_pair_index], ax=ax, s=marker_size, linewidth=0,
                            marker=markers[detector_pair_index], label=label_in_plot)
            sns.scatterplot(y=detcs_excld_data[data_label]["y"], x=detcs_excld_data[data_label]["x"], color=excluded_color,
                            ax=ax, s=marker_size, linewidth=0, marker=markers[detector_pair_index], label=label_in_plot + " excl.")
        plt.legend(ncol=2)
    else:
        if not use_pts_white_border:
            sns.scatterplot(x=x_data_label, y=y_data_label, hue = conditional_label, data = data_frame,
                                   ax=ax, s=marker_size, linewidth=0)
        else:
            sns.scatterplot(x=x_data_label, y=y_data_label, hue=conditional_label, style=conditional_style, data=data_frame,
                                   ax=ax, s=marker_size)


    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    if ymin is not None:
        if ymax is not None:
            plt.ylim(ymin, ymax)

    if xmin is not None:
        if xmax is not None:
            plt.xlim(xmin, xmax)

    if label_cms_status:
        add_extra_text(ax, fig_size_shape, energy_year_label=energy_year_label,
                       experiment=setts.experiment, work_status=setts.work_status)

    ax.set_title(title)
    ax.set_ylabel(ylabel, labelpad=setts.axis_labelpad, weight=setts.axis_weight, size=setts.axis_case_size[fig_size_shape])
    ax.set_xlabel(xlabel, labelpad=setts.axis_labelpad, weight=setts.axis_weight, size=setts.axis_case_size[fig_size_shape])

    plt.xticks(fontsize=setts.axis_thicks_case_size[fig_size_shape])
    plt.yticks(fontsize=setts.axis_thicks_case_size[fig_size_shape])

    # margin optimization
    if fig_size == (12, 4):
        plt.subplots_adjust(left=0.1, right=0.97, top=0.9, bottom=0.2)

    if xlabel_rotation:
        plt.xticks(rotation=xlabel_rotation)

    if leg_col is not None:
        plt.legend(ncol=leg_col)


    return fig



def snsplot_hist_all_and_excluded(data_frame, x_data_label, conditional_label, bins, xmin, xmax,
                                  title='', xlabel='', ylabel='',
                                  ymin=None, ymax=None,
                                  label_cms_status=True, energy_year_label='', fig_size_shape='sq'):

    fig_size = get_fig_size(fig_size_shape)
    fig, ax = plt.subplots(figsize=fig_size)

    keys = np.unique(data_frame[conditional_label])
    histos_data = {}

    for key in keys:
        histos_data[key] = []

    for index_data in range(0, len(data_frame)):
        key = data_frame[conditional_label][index_data]
        histos_data[key].append(data_frame[x_data_label][index_data])

    for key in keys:
        data_to_plot = np.array(histos_data[key])
        data_to_plot = data_to_plot[~np.isnan(data_to_plot)]
        data_to_plot = data_to_plot[(data_to_plot >= xmin) & (data_to_plot <= xmax)]
        sns.distplot(data_to_plot, ax=ax, kde=False, bins=bins, label=key)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    if ymin is not None:
        if ymax is not None:
            plt.ylim(ymin, ymax)

    plt.xlim(xmin, xmax)

    if label_cms_status:
        add_extra_text(ax, fig_size_shape, energy_year_label=energy_year_label,
                       experiment=setts.experiment, work_status=setts.work_status)

    ax.set_title(title)
    ax.set_ylabel(ylabel, labelpad=setts.axis_labelpad, weight=setts.axis_weight, size=setts.axis_case_size[fig_size_shape])
    ax.set_xlabel(xlabel, labelpad=setts.axis_labelpad, weight=setts.axis_weight, size=setts.axis_case_size[fig_size_shape])

    plt.xticks(fontsize=setts.axis_thicks_case_size[fig_size_shape])
    plt.yticks(fontsize=setts.axis_thicks_case_size[fig_size_shape])

    plt.legend()

    return fig


def add_extra_text(ax, plot_frame_ratio, energy_year_label='', experiment='', work_status=''):
    if plot_frame_ratio == 'nsq':
        plt.text(__year_energy_label_pos_nsq[0], __year_energy_label_pos_nsq[1], str(energy_year_label), ha='left',
                 fontsize=setts.year_energy_font_size, fontweight='bold', transform=ax.transAxes)
        plt.text(__cms_label_pos_nsq[0], __cms_label_pos_nsq[1], str(experiment), ha='left',
                 fontsize=setts.experiment_font_size, fontweight='bold', transform=ax.transAxes)
        plt.text(__cms_label_pos_nsq[0] + __delta_y_pos_nsq, __cms_label_pos_nsq[1], str(work_status), ha='left',
                 fontsize=setts.experiment_font_size, fontweight='light', transform=ax.transAxes)

    elif plot_frame_ratio == 'sq':
        plt.text(__year_energy_label_pos_sq[0], __year_energy_label_pos_sq[1], str(energy_year_label), ha='left',
                 fontsize=setts.year_energy_font_size, fontweight='bold', transform=ax.transAxes)
        plt.text(__cms_label_pos_sq[0], __cms_label_pos_sq[1], str(experiment), ha='left',
                 fontsize=setts.experiment_font_size, fontweight='bold', transform=ax.transAxes)
        plt.text(__cms_label_pos_sq[0] + __delta_y_pos_sq, __cms_label_pos_sq[1], str(work_status), ha='left',
                 fontsize=setts.experiment_font_size, fontweight='light', transform=ax.transAxes)
    else:
        raise ValueError('frame ratio value not implemented')


def plot_plt_scatter(x, y):
    fig, ax = plt.subplots()
    plt.scatter(x, y)

    return fig


def sns_regresion_plot(data, x_data_label, y_data_label):
    sns.set(color_codes=True)
    fig = sns.lmplot(x=x_data_label, y=y_data_label, data=data)

    return fig


def get_fig_size(shape: str) -> tuple:
    return setts.fig_sizes[shape]