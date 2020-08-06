import numpy as np
import pandas as pd
import plotly.express as px
import settings as setts
import tools.lumi_tools as ltools


def live_scatter_from_pandas(data_frame: pd.DataFrame, x_data_label, y_data_label, show_also_info_in=[],
                             title='', xlabel='', ylabel='',
                             ymin=None, ymax=None, xmin=None, xmax=None):
    fig = px.scatter(data_frame, x=x_data_label, y=y_data_label, hover_data=show_also_info_in,
                     title=title)
    fig.update_yaxes(title=ylabel)
    fig.update_xaxes(title=xlabel)
    fig.show()


def live_line_from_pandas(data_frame: pd.DataFrame, x_data_label, y_data_label, show_also_info_in=[],
                          title='', xlabel='', ylabel='',
                          ymin=None, ymax=None, xmin=None, xmax=None):
    fig = px.line(data_frame, x=x_data_label, y=y_data_label, hover_data=show_also_info_in,
                  title=title)
    fig.update_yaxes(title=ylabel)
    fig.update_xaxes(title=xlabel)

    fig.show()
