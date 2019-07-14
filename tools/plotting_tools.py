import matplotlib.pyplot as plt
import settings as setts
import tools.lumi_tools as ltools
from statsmodels.graphics.api import abline_plot
import seaborn as sns
import numpy as np

__cms_label_pos_sq = setts.cms_label_pos_sq
__cms_label_pos_nsq = setts.cms_label_pos_nsq
__delta_y_pos_sq = setts.delta_y_pos_sq
__delta_y_pos_nsq = setts.delta_y_pos_nsq

__year_energy_label_pos_sq = setts.year_energy_label_pos_sq
__year_energy_label_pos_nsq = setts.year_energy_label_pos_nsq


def save_py_fig_to_file(fig, output_name):
    fig.savefig(output_name)


def save_plots(names_and_plots, output_name):
    ltools.check_and_create_folder(output_name)
    for plot in names_and_plots:
        plot_name = plot[0]
        plot_object = plot[1]
        for plot_format in setts.plots_formats:
            plot_full_path = output_name + plot_name + '.' + plot_format
            save_py_fig_to_file(plot_object, plot_full_path)
            print('saved plot: ' + plot_full_path)


# TODO: set histo range for same binning
def hist_from_pandas_frame(data_frame, col_label, nbins, title='', xlabel='', ylabel='', xmin=0.0, xmax=0.0,
                           weight_label=None,
                           mean=None, stdv=None, err_mean=None, color='#86bf91',
                           label_cms_status=True, energy_year_label='',
                           fig_size_shape='sq'):
    fig_size = get_fig_size(fig_size_shape)
    if weight_label:
        ratio_hist = data_frame.hist(bins=nbins, column=col_label, grid=False, weights=data_frame[weight_label],
                                     figsize=fig_size, sharex=True, color=color, zorder=2, rwidth=0.9)
    else:
        ratio_hist = data_frame.hist(bins=nbins, column=col_label, grid=False,
                                     figsize=fig_size, sharex=True, color=color, zorder=2, rwidth=0.9)

    ratio_hist_ax = ratio_hist[0][0]

    ratio_hist_ax.set_title(title)
    ratio_hist_ax.set_xlabel(xlabel, labelpad=setts.axis_labelpad, weight=setts.axis_weight, size=setts.axis_case_size)
    ratio_hist_ax.set_ylabel(ylabel, labelpad=setts.axis_labelpad, weight=setts.axis_weight, size=setts.axis_case_size)
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
    ax.set_xlabel(xlabel, labelpad=setts.axis_labelpad, weight=setts.axis_weight, size=setts.axis_case_size)
    ax.set_ylabel(ylabel, labelpad=setts.axis_labelpad, weight=setts.axis_weight, size=setts.axis_case_size)
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
                                   leg_text_s=setts.leg_vs_plots_text_s,
                                   fig_size_shape='nsq'):
    fig_size = get_fig_size(fig_size_shape)
    if xlabel == '':
        xlabel = x_data_label
    if ylabel == '':
        ylabel = y_data_label
    fig, ax = plt.subplots()
    if type(y_data_label) == list:
        plot = data_frame.plot(x=x_data_label, y=y_data_label, style='o', figsize=fig_size,
                               markersize=0.5, ax=ax)
        if legend_labels is None:
            legend_labels = y_data_label
        else:
            assert len(legend_labels) == len(y_data_label)
        ax.legend(legend_labels, ncol=len(legend_labels),
                  markerscale=leg_marker_sc, fontsize=leg_text_s, loc=legend_position)
    else:
        plot = data_frame.plot(x=x_data_label, y=y_data_label, style='o', figsize=fig_size,
                               markersize=0.5, legend=None, ax=ax)
    plot.set_title(title)
    plot.set_ylabel(ylabel, labelpad=setts.axis_labelpad, weight=setts.axis_weight, size=setts.axis_case_size)
    plot.set_xlabel(xlabel, labelpad=setts.axis_labelpad, weight=setts.axis_weight, size=setts.axis_case_size)

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

    ax.set_xlabel(xlabel, labelpad=setts.axis_labelpad, weight=setts.axis_weight, size=setts.axis_case_size)
    ax.set_ylabel(ylabel, labelpad=setts.axis_labelpad, weight=setts.axis_weight, size=setts.axis_case_size)

    add_extra_text(ax, fig_size_shape, energy_year_label=energy_year_label,
                   experiment=setts.experiment, work_status=setts.work_status)

    # Adding fitting info
    plt.text(setts.pos_by_fill_fit_info[0], setts.pos_by_fill_fit_info[1],
             r"slope[hz/" + r"$\mu$" + "b]: " + str(float("{0:.4f}".format(fitted_slope))) +
             r' $\pm$ ' + str(float("{0:.4f}".format(fitted_slope_err))),
             ha='left',
             fontsize=setts.leg_font_size, fontweight='bold', transform=ax.transAxes)

    return fig


def add_extra_text(ax, plot_frame_ratio, energy_year_label='', experiment='', work_status=''):
    if plot_frame_ratio == 'nsq':
        plt.text(__year_energy_label_pos_nsq[0], __year_energy_label_pos_nsq[1], str(energy_year_label), ha='left',
                 fontsize=setts.leg_font_size, fontweight='bold', transform=ax.transAxes)
        plt.text(__cms_label_pos_nsq[0], __cms_label_pos_nsq[1], str(experiment), ha='left',
                 fontsize=setts.leg_font_size, fontweight='bold', transform=ax.transAxes)
        plt.text(__cms_label_pos_nsq[0], __cms_label_pos_nsq[1] - __delta_y_pos_nsq, str(work_status), ha='left',
                 fontsize=setts.leg_font_size, fontweight='light', transform=ax.transAxes)

    elif plot_frame_ratio == 'sq':
        plt.text(__year_energy_label_pos_sq[0], __year_energy_label_pos_sq[1], str(energy_year_label), ha='left',
                 fontsize=setts.leg_font_size, fontweight='bold', transform=ax.transAxes)
        plt.text(__cms_label_pos_sq[0], __cms_label_pos_sq[1], str(experiment), ha='left',
                 fontsize=setts.leg_font_size, fontweight='bold', transform=ax.transAxes)
        plt.text(__cms_label_pos_sq[0], __cms_label_pos_sq[1] - __delta_y_pos_sq, str(work_status), ha='left',
                 fontsize=setts.leg_font_size, fontweight='light', transform=ax.transAxes)
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
