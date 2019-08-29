import settings as setts
import os
import numpy as np
from datetime import datetime
from statsmodels.stats.weightstats import DescrStatsW
import statsmodels.api as sm
import pandas as pd


def get_year_and_energy(fill):
    range_fill_yearly = setts.rangeFillYearly
    years_energy = list(range_fill_yearly)
    years_energy.sort()
    year = 0
    energy = 0

    for tag in years_energy:
        if in_range(int(fill), range_fill_yearly[tag][0], range_fill_yearly[tag][1]):
            year = tag[0]
            energy = tag[1]
            break
    if year == 0:
        print("Initial fill not found in official fill ranges. Automatic yearly functions wont work!!")
    if energy == 0:
        print("Initial fill not found in official fill ranges. Automatic energy functions wont work!!")

    return year, energy


def in_range(val, minval, maxval):
    in_range_out = False
    if minval <= val <= maxval:
        in_range_out = True
    return in_range_out


def copy_data_columns(pd_copy_from: object, pd_copy_to: object, labels_list: object) -> None:
    for label in labels_list:
        pd_copy_to[label] = pd_copy_from[label]


def check_and_create_folder(folder_path, creation_info=True):
    try:
        os.makedirs(folder_path)
        print('output folder has been created: ' + folder_path)
    except:
        if creation_info:
            print(folder_path + ' Already exists -->> Content will be overwritten.')


def get_lumi_unit_from_csv(file_path):
    csv_file = open(file_path, 'r')
    header_line = csv_file.readlines()[1]
    header_line_array = header_line.split(',')
    if header_line_array[0] == '#run:fill':
        lumi_header = header_line_array[5]
        return lumi_header.split('(')[1].split(')')[0]
    else:
        raise AssertionError('header line not found!')


def get_lumi_unit_convertion_factor_to_inv_fb(unit_from):
    conv_factor = 1.00
    # convert all to 1/fb
    if unit_from == 'hz/ub':
        conv_factor = conv_factor / 1000000000.0

    return np.float64(conv_factor * setts.xLS)


def add_date_column(pandas_frm):
    date_array = []
    for i_time in pandas_frm['time']:
        date_array.append(np.datetime64(datetime.utcfromtimestamp(i_time)))

    assert len(date_array) == len(pandas_frm['time'])

    pandas_frm['date'] = np.array(date_array)


def check_and_clean_after_merging(pd_df):
    cols_x = []
    cols_y = []
    all_cols = list(pd_df)
    for col in all_cols:
        if '_x' in col:
            cols_x.append(col)
            col_y_name = col.replace('_x', '_y')
            col_name = col.replace('_x', '')
            if col_y_name in all_cols:
                nan_vals_y = pd_df[col_y_name].isnull().values
                full_vals = []
                for iii in range(0, len(pd_df)):
                    if nan_vals_y[iii]:
                        full_vals.append(pd_df[col][iii])
                    else:
                        full_vals.append(pd_df[col_y_name][iii])
                pd_df[col_name] = np.array(full_vals)
                not_needed_columns = [col_y_name, col]
                pd_df.drop(columns=not_needed_columns, inplace=True)
            else:
                raise AssertionError('Cleaning after merging not possible!')


def get_w_stats(vals_array, w_array, min_val=None, max_val=None):
    # dropping NaN values:
    temp_list = []
    temp_list_w = []
    if min_val is not None and max_val is not None:
        for ratio_index in range(0, len(vals_array)):
            ratio_val = vals_array[ratio_index]
            if ratio_val >= min_val and ratio_val <= max_val:
                temp_list.append(ratio_val)
                temp_list_w.append(w_array[ratio_index])
    vals_array=np.array(temp_list)
    w_array=np.array(temp_list_w)
    w_array = w_array[~np.isnan(w_array)]
    vals_array = vals_array[~np.isnan(vals_array)]
    assert len(vals_array) == len(w_array)

    w_stats = DescrStatsW(vals_array, weights=w_array)

    return w_stats


def get_linear_model_from_pd_cols(data: pd.DataFrame, x_col_name: str, y_col_name: str) -> sm.OLS:
    x = data[x_col_name]
    y = data[y_col_name]

    return sm.OLS(y, x)


def convert_detector_name(name : str) -> str:
    exit_label = name
    if name == "PXL":
        exit_label = "pcc"
    elif name == "PLTZERO":
        exit_label = "plt"
    elif name == "RAMSES":
        exit_label = "ram"

    return exit_label

def keys_exists(element, *keys):
    '''
    Check if *keys (nested) exists in `element` (dict).
    '''
    if not isinstance(element, dict):
        raise AttributeError('keys_exists() expects dict as first argument.')
    if len(keys) == 0:
        raise AttributeError('keys_exists() expects at least two arguments, one given.')

    _element = element
    for key in keys:
        try:
            _element = _element[key]
        except KeyError:
            return False
    return True
