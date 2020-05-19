data_base_folder = "beamData/"
conf_label_data_file_path = "data_file_path"
conf_label_offset_values = "offset_values"
conf_label_offset_time = "offset_time"
conf_label_min_time = "min_time"
conf_label_max_time = "max_time"
conf_label_zero_time = "zero_time"
conf_label_special_time_intervals = "special_time_intervals"
conf_label_y_range = "y_range"
conf_label_y_diff_range = "y_diff_range"

l_min_plot = -100.
l_max_plot = 100.
l_diff_min_plot = -20.
l_diff_max_plot = 20.
l_diff_zoom_min_plot = -5.
l_diff_zoom_max_plot = 5.

save_figures_as_pickle = True

col_pos_offset_array = {
    "B1_H": 0,
    "B1_V": 1,
    "B2_H": 2,
    "B2_V": 3
}

config_dict = {
    "DOROS": {
        4954: {
            conf_label_data_file_path: "DOROS_Fill4954_Timber.csv",
            conf_label_offset_values: [[204.35, -25.65, 68.65, -231.2], [206.9, -20.55, 67.45, -229],
                                       [205.45, -19.75, 69.45, -227.9]],
            conf_label_offset_time: [1464336119, 1464348464, 1464350758],
            conf_label_min_time: 1464336193,
            conf_label_max_time: 1464353223,
            # conf_label_zero_time: 1464335770,

            # conf_label_special_time_intervals: [[75, 100], [110, 180], [180, 215]]

            conf_label_y_range: [-600., 600.],
            # conf_label_y_diff_range: [-30., 30.]
            conf_label_y_diff_range: [-7., 7.]

        },
        7442: {
            conf_label_data_file_path: "DOROS_Fill7442_Timber_v02.csv",
            conf_label_offset_values: [[318.65, -1710.15, 172.35, -1722.15], [324.65, -1701.95, 172.95, -1722.15],
                                       [324.60, -1702.00, 172.55, -1722.85], [314.55, -1703.10, 172.00, -1721.25],
                                       [316.95, -1704.15, 172.10, -1721.25], [317.45, -1704.15, 152.05, -1713.45],
                                       [316.6, -1704.35, 152.70, -1714.00], [315.90, -1704.35, 152.70, -1714.00],
                                       [315.90, -1704.35, 152.90, -1714.85], [315.80, -1702.95, 152.75, -1712.70],
                                       [314.25, -1701.75, 152.65, -1712.00]],
            conf_label_offset_time: [1542135000, 1542135708, 1542136173, 1542139644, 1542140214, 1542140579,
                                     1542141318, 1542143712, 1542144289, 1542145836, 1542148211],
            conf_label_min_time: 1542135000,
            conf_label_max_time: 1542150000,
            # conf_label_zero_time: 1542135660,

            conf_label_special_time_intervals: [[75, 100], [110, 180], [180, 215]]

        },
        7443: {
            conf_label_data_file_path: "DOROS_Fill7443_Timber_v02.csv",
            conf_label_offset_values: [[320.40, -1662.80, 161.20, -1682.70], [320.40, -1662.90, 163.10, -1683.55],
                                       [324.55, -1664.90, 162.80, -1683.50], [316.80, -1666.55, 160.10, -1680.95],
                                       [319.05, -1660.20, 160.50, -1680.95], [318.95, -1660.25, 159.00, -1682.15],
                                       [318.10, -1659.50, 158.60, -1681.65], [321.90, -1654.20, 158.05, -1679.75],
                                       [326.45, -1656.15, 157.85, -1679.40], [326.35, -1655.05, 157.75, -1679.50],
                                       [326.80, -1655.05, 139.55, -1673.00], [323.10, -1653.00, 141.10, -1670.75],
                                       [321.85, -1652.40, 142.10, -1671.00], [320.80, -1651.70, 143.30, -1670.40],
                                       [309.75, -1640.15, 143.40, -1670.35]],
            conf_label_offset_time: [1542211800, 1542213059, 1542213602, 1542215346, 1542215382, 1542216031,
                                     1542216740, 1542220016, 1542221586, 1542221877, 1542221936, 1542228161,
                                     1542229709, 1542231642, 1542232040],
            conf_label_min_time: 1542211800,
            conf_label_max_time: 1542233400,
            # conf_label_zero_time: 1542211800
        }

    },

    "Nominal": {
        4954: {
            conf_label_data_file_path: "vdm_Nominal_Fill4954_n02.csv",
            # conf_label_special_time_intervals: [[75, 100], [110, 180], [180, 215]],

            conf_label_min_time: 1464336193,
            conf_label_max_time: 1464353223,
            # conf_label_zero_time: 1464335770,

            conf_label_y_range: [-600., 600.],
        },
        7442: {
            conf_label_data_file_path: "vdm_Nominal_Fill7442_n03.csv",
            conf_label_special_time_intervals: [[75, 100], [110, 180], [180, 215]],

            conf_label_min_time: 1542135000,
            conf_label_max_time: 1542150000,
            # conf_label_zero_time: 1542135660,
        },
        7443: {
            conf_label_data_file_path: "vdm_Nominal_Fill7443_n01.csv",

            conf_label_min_time: 1542211800,
            conf_label_max_time: 1542233400,
            # conf_label_zero_time: 1542211800
        }

    },

    "arcBPM": {
        7442: {
            conf_label_data_file_path: "CMS_arc_Fill7442_pos_merged.csv",
        }
    }
}

# -*-*-*-*-*-*-*-*-*-*-*- > > > Other Settings and helpers < < < -*-*-*-*-*-*-*-*-*-*-*-
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'