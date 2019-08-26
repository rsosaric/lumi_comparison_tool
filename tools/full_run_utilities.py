import numpy as np
from tools.luminometer import Luminometer as L
import pandas as pd


def merge_lumies(name: str, data_file_name: list, all_data: bool = False,
                 mixed_data=False) -> L:

    lumi_merge:L = L(name, data_file_name[0] + name + ".csv", all_data=all_data, mixed_data=mixed_data)
    data_frame_list = [lumi_merge.data]

    for i in range(1, len(data_file_name), 1):
        data_frame_list.append(L(name, data_file_name[i] + name + ".csv", all_data=all_data, mixed_data=mixed_data).data)

    lumi_merge.data = pd.concat(data_frame_list)

    return lumi_merge


