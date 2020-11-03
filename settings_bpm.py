data_base_folder = "beamData/"
conf_label_data_file_path = "data_file_path"
conf_label_offset_values = "offset_values"
conf_label_offset_time = "offset_time"
conf_label_optimization_scans_times = "optimization scans times"
conf_label_min_time = "min_time"
conf_label_max_time = "max_time"
conf_label_zero_time = "zero_time"
conf_label_special_time_intervals = "special_time_intervals"
conf_label_y_range = "y_range"
conf_label_y_diff_range = "y_diff_range"
conf_label_y_diff_smaller_range = "y_diff_range_smaller"
conf_label_length_scale = "length_scale"  # [x,y] corrections
conf_label_length_scale_per_beam_epsilon = "length_scale per_beam_epsilon"  # [[B1_x, B2_x],[B1_y, B2_y]]
conf_label_apply_length_scale_per_beam_epsilon_from_file = "length_scale per_beam_epsilon loaded from file"
conf_limits_for_length_scale_per_beam_epsilon = "limits_for_length_scale_per_beam_epsilon"  # [min, max]
conf_exclusion_from_length_scale_per_beam_epsilon = "conf_exclusion_from_length_scale_per_beam_epsilon"  # {scan:beam,.}
conf_label_beam_deflection = "beam_deflection"
conf_label_beam_deflection_gfactor = "beam_deflection_gfactor"  # geometrical factor for beam_deflection correction
conf_label_special_col_sign = "special col sign"  # it can be used to invert input data column sign or scale it
conf_label_orbit_drift = "orbit_drift"  # Orbit drift input file
conf_label_deep_studies = "deep_studies"  # Do extra studies and plots: LR studies
conf_label_compute_offsets = "compute offsets"  # Compute offsets internally
conf_label_compute_offsets_method = "compute offsets method"  # Compute offsets method: only after optimization (more
# precised)(default), all head-on (better for initial plots)
conf_label_special_input_col_names = "special_input_col_names"  # Use if the format of the input file is different
conf_label_unit_position_factor_to_um = "factor to convert input positions to [um]"
conf_label_timestamp_in_ms = "timestamp_in_ms"  # Use if the input times are in ms
conf_label_scans_names = "scans names"  # Scans name list. Will be used in plots
conf_label_scans_time_windows = "scans times"  # Scans times limits.
conf_label_exclusion_times = "times to exclude from results"  # Optimization scans, not useful times during VDM, etc
conf_label_time_window_for_offsets = "time_window_for_offsets"
conf_label_orbit_drift_plot_range = "orbit_drift_plot_range"
conf_compute_length_scale_epsilon_fits = "compute_length_scale_epsilon_fits"  # per scan
conf_label_scans_with_inverted_beam_sign = "list of scans where the corrections like deflection or orbit drift " \
                                           "needs an opposite sign"

debugging = False

l_min_plot = -100.
l_max_plot = 100.
l_diff_min_plot = -20.
l_diff_max_plot = 20.
l_diff_orbit_drift_min_plot = -40.
l_diff_orbit_drift_max_plot = 40.
l_diff_zoom_min_plot = -5.
l_diff_zoom_max_plot = 5.

zero_position_window = 0.1     # in um
zero_position_window_for_offsets = 0.1     # in um
delta_scan_time_window = 100   # seconds
time_window_for_offsets = 150  # seconds

delta_pos_for_scan_labels = 0.15  # pos_y = min_y + min_y * delta_pos_for_scan_labels

save_figures_as_pickle = False

col_pos_offset_array = {
    "B1_H": 0,
    "B1_V": 1,
    "B2_H": 2,
    "B2_V": 3
}

config_dict = {
    "Nominal": {
        4266: {
            conf_label_data_file_path: "vdm_Nominal_Fill4266_n02.csv",
            conf_label_length_scale: [0.9952, 0.9969],
            conf_label_beam_deflection: "beamData/extra_corrections/deflectioncorrection4266.csv",
            conf_label_beam_deflection_gfactor: [2.65, 2.76],  # [x, y]

            conf_label_scans_names: ["X1", "Y1", "Y2", "X2", "L1X", "L1Y", "X3", "Y3", "X4", "Y4", "L2X", "X5", "Y5"],
            conf_label_scans_time_windows: [[1440463385, 1440464388], [1440464725, 1440465732],
                                            [1440466016, 1440467017], [1440467318, 1440468322],
                                            [1440468604, 1440470118], [1440470347, 1440471670],
                                            [1440472571, 1440473372], [1440473655, 1440474470],
                                            [1440474775, 1440475585], [1440475865, 1440476670],
                                            [1440476953, 1440478244], [1440478566, 1440479567],
                                            [1440479880, 1440480885]],

            conf_label_optimization_scans_times: [1440462850, 1440472132],

            # conf_label_exclusion_times: [[1440471893, 1440472149]],

            conf_label_min_time: 1440462850,
            conf_label_max_time: 1440480981,
            conf_label_zero_time: 1440462850,
            conf_label_special_time_intervals: [[5, 95]],

            conf_label_y_range: [-650., 650.],
            conf_label_y_diff_range: [-7., 7.]
        },
        4689: {
            conf_label_data_file_path: "vdm_Nominal_Fill4689_n01.csv",

            conf_label_scans_names: ["X1", "Y1", "Y2", "X2", "X3", "Y3", "Y4", "X4", "X5", "Y5"],
            conf_label_scans_time_windows: [[1449134119, 1449135013], [1449135216, 1449136109],
                                            [1449136443, 1449137337], [1449137529, 1449138424],
                                            [1449154166, 1449154670], [1449154731, 1449155279],
                                            [1449155460, 1449155961], [1449156079, 1449156581],
                                            [1449158467, 1449158972], [1449159125, 1449159628]],

            conf_label_optimization_scans_times: [1449133477, 1449153809],

            conf_label_exclusion_times: [[1449133337, 1449133477], [1449153552, 1449153809], [1449156746, 1449158261],
                                         [1449138606, 1449153552]],

            conf_label_min_time: 1449134042,
            conf_label_max_time: 1449163025,
            conf_label_zero_time: 1449134042,
            # conf_label_special_time_intervals: [[5, 95]],

            conf_label_y_range: [-100., 100.],
            conf_label_y_diff_range: [-7., 7.]
        },
        4954: {
            conf_label_data_file_path: "vdm_Nominal_Fill4954_n02.csv",
            # conf_label_special_time_intervals: [[5, 80]],
            conf_label_length_scale: [0.9899, 0.9950],
            # conf_label_length_scale: [1.0, 1.0],

            conf_label_beam_deflection: "beamData/extra_corrections/deflectioncorrection4954.csv",
            conf_label_beam_deflection_gfactor: [2.65, 2.76],  # [x, y]

            conf_label_scans_names: ["X1", "Y1", "Y2", "X2", "I1X", "I1Y", "I2X", "I2Y", "L1X", "L1Y", "X3", "Y3"],
            # conf_label_scans_with_inverted_beam_sign: ["L1X", "L1Y"],
            conf_label_scans_time_windows: [[1464336285, 1464337268], [1464337522, 1464338508],
                                            [1464338769, 1464339757], [1464339948, 1464340930],
                                            [1464341247, 1464342231], [1464342677, 1464343658],
                                            [1464344429, 1464345408], [1464345652, 1464346642],
                                            [1464347152, 1464348181], [1464348634, 1464349677],
                                            [1464350869, 1464351851], [1464352142, 1464353130]],

            conf_label_optimization_scans_times: [1464336119, 1464348464, 1464350758],

            conf_label_min_time: 1464336193,
            conf_label_max_time: 1464353223,
            conf_label_zero_time: 1464336193,
            # conf_label_zero_time: 1464335770,

            conf_label_y_range: [-600., 600.],
        },
        7442: {
            conf_label_data_file_path: "vdm_Nominal_Fill7442_n03.csv",
            conf_label_special_time_intervals: [[75, 100], [110, 180], [180, 215]],
            conf_label_scans_names: ["emX1", "emY1", "L1Y", "L1X", "X1", "Y1", "offY1", "offX1", "offX2", "XY1", "XY2"],
            conf_label_scans_time_windows: [[1542136154, 1542136383], [1542136391, 1542136620],
                                            [1542137597, 1542138346], [1542138422, 1542139170],
                                            [1542139705, 1542140192], [1542140262, 1542140752],
                                            [1542141840, 1542143099], [1542143108, 1542144353],
                                            [1542144358, 1542145619], [1542146035, 1542146905],
                                            [1542146980, 1542147848]],


            conf_label_min_time: 1542135000,
            conf_label_max_time: 1542150000,
            # conf_label_zero_time: 1542135660,
        },
        7443: {
            conf_label_data_file_path: "vdm_Nominal_Fill7443_n01.csv",
            conf_label_scans_names: ["emX2", "emY2", "X2", "Y2", "X3", "Y3", "offY3", "offX3",
                                     "XY3", "XY4", "X4", "Y4", "X5", "Y5", "XY5", "XY6", "emX3", "emY3", "X6", "6Y"],
            conf_label_scans_time_windows: [[1542213090, 1542213255], [1542213263, 1542213428],
                                            [1542213934, 1542214423], [1542214521, 1542215000],
                                            [1542215492, 1542215981], [1542216058, 1542216536],
                                            [1542217148, 1542218409], [1542218434, 1542219682],
                                            [1542220101, 1542220579], [1542221107, 1542221585],
                                            [1542222008, 1542222498], [1542222576, 1542223054],
                                            [1542228435, 1542228923], [1542229018, 1542229498],
                                            [1542229911, 1542230390], [1542230465, 1542230944],
                                            [1542231672, 1542231834], [1542231840, 1542231995],
                                            [1542232313, 1542232543], [1542232550, 1542232780]],

            conf_label_min_time: 1542211800,
            conf_label_max_time: 1542233400,
            # conf_label_zero_time: 1542211800
        },
        7483: {
            conf_label_data_file_path: "vdm_Nominal_Fill7483_n01.csv",
            conf_label_scans_names: ["L1Y", "L1X", "X1", "Y1"],
            conf_label_scans_time_windows: [[1543464250, 1543465002], [1543465076, 1543465827],
                                            [1543466381, 1543466611], [1543466693, 1543466924]],

            conf_label_min_time: 1543463700,
            conf_label_max_time: 1543467300
        },

    },

    "DOROS": {
        4266: {
            conf_label_data_file_path: "DOROS_Fill4266_Timber.csv",
            conf_label_offset_time: [1440462850, 1440472132],
            conf_label_offset_values: [[-87.8, -157.85, -291.55, -457.7],
                                       [-84.05, -162.2, -281.85, -461.7]],

            conf_label_min_time: 1440462850,
            conf_label_max_time: 1440480981,

            # conf_label_special_time_intervals: [[5, 95], [5, 27], [28, 50]],

            conf_label_length_scale: [0.9941, 0.9943],
            conf_label_apply_length_scale_per_beam_epsilon_from_file: True,
            conf_exclusion_from_length_scale_per_beam_epsilon: {"L1X": ["B1", "B2"], "L1Y": ["B1", "B2"],
                                                                "L2X": ["B1", "B2"],
                                                                "X3": ["B2"], "Y3": ["B2"],
                                                                "X4": ["B1"], "Y4": ["B1"]},
            conf_limits_for_length_scale_per_beam_epsilon: [-0.03, 0.03],
            # conf_label_beam_deflection: "beamData/extra_corrections/deflectioncorrection4266.csv",
            # conf_label_beam_deflection_gfactor: [2.65, 2.76],  # [x, y]

            # conf_label_orbit_drift: "beamData/extra_corrections/orbitdriftcorrection4266average.csv",
            conf_label_orbit_drift: "beamData/extra_corrections/orbitdriftcorrection4266doros.csv",

            conf_label_y_range: [-650., 650.],
            # conf_label_y_diff_range: [-8., 8.],
            conf_label_y_diff_range: [-10., 10.],
            conf_label_y_diff_smaller_range: [-8, 8],

            conf_label_special_col_sign: {
                "B1_H_L": -1, "B1_H_R": -1, "B2_H_R": -1
            },

            conf_compute_length_scale_epsilon_fits: True,

            # conf_label_data_file_path: "4266_P5_data_exper.csv",
            # conf_label_special_input_col_names: {"time": "timestamp",
            #                                      "pc H.B1.L [um]": "B1_H_L",
            #                                      "pc V.B1.L [um]": "B1_V_L",
            #                                      "pc H.B2.L [um]": "B2_H_L",
            #                                      "pc V.B2.L [um]": "B2_V_L",
            #                                      "pc H.B1.R [um]": "B1_H_R",
            #                                      "pc V.B1.R [um]": "B1_V_R",
            #                                      "pc H.B2.R [um]": "B2_H_R",
            #                                      "pc V.B2.R [um]": "B2_V_R"},
            #
            # conf_label_unit_position_factor_to_um: 1.,
            # conf_label_timestamp_in_ms: True,
            # conf_label_compute_offsets: True,
            # conf_label_time_window_for_offsets: 600,  #s

            # conf_label_deep_studies: True,

        },
        4689: {
            conf_label_data_file_path: "4689_DOROS_raw_data_from_Timber.csv",

            conf_label_offset_time: [1449133477, 1449153809],
            conf_label_offset_values: [[-83.85, -191.55, 2072.3, -339.3], [-74.8, -193.85, 2120.8, -343.55]],

            conf_label_min_time: 1449134042,
            conf_label_max_time: 1449163025,
            conf_label_zero_time: 1449134042,
            # conf_label_special_time_intervals: [[5, 95]],

            conf_label_y_range: [-100., 100.],
            conf_label_y_diff_range: [-20., 20.],

            conf_label_special_col_sign: {
                "B1_H_L": -1, "B1_V_L": -1, "B2_H_L": -1, "B2_V_L": -1
            },

            conf_label_orbit_drift_plot_range: [-15., 40.],

            # conf_label_deep_studies: True,

            # conf_label_compute_offsets: True,
            # conf_label_time_window_for_offsets: 600,  #s
        },
        4954: {
            conf_label_data_file_path: "DOROS_Fill4954_Timber.csv",
            conf_label_offset_time: [1464336119, 1464348464, 1464350758],
            # conf_label_compute_offsets: True,
            conf_label_offset_values: [[204.35, -25.65, 68.65, -231.2], [206.9, -20.55, 67.45, -229],
                                       [205.45, -19.75, 69.45, -227.9]],

            conf_label_min_time: 1464336193,
            conf_label_max_time: 1464353223,
            conf_label_zero_time: 1464336193,

            conf_label_length_scale: [1.0315, 1.0314],
            conf_label_apply_length_scale_per_beam_epsilon_from_file: True,
            conf_exclusion_from_length_scale_per_beam_epsilon: {"I1X": ["B1"], "I2X": ["B2"],
                                                                "I1Y": ["B1"], "I2Y": ["B2"]},
            conf_limits_for_length_scale_per_beam_epsilon: [-0.03, 0.03],
            # conf_label_length_scale_per_beam_epsilon: [[-0.0044, 0.0118], [-0.0017, 0.0095]],

            # conf_label_orbit_drift: "beamData/extra_corrections/orbitdriftcorrection4954average.csv",
            conf_label_orbit_drift: "beamData/extra_corrections/orbitdriftcorrection4954doros.csv",

            # conf_label_zero_time: 1464335770,

            # conf_label_special_time_intervals: [[0, 80], [0, 20]],

            conf_label_y_range: [-600., 600.],
            # conf_label_y_diff_range: [-30., 30.]
            conf_label_y_diff_range: [-5., 8.],
            conf_label_y_diff_smaller_range: [-5.5, 5.5],
            conf_compute_length_scale_epsilon_fits: True

        },
        7442: {
            conf_label_data_file_path: "DOROS_Fill7442_Timber_v02.csv",
            conf_label_offset_time: [1542135000, 1542135708, 1542136173, 1542139644, 1542140214, 1542140579,
                                     1542141318, 1542143712, 1542144289, 1542145836, 1542148211],
            conf_label_offset_values: [[318.65, -1710.15, 172.35, -1722.15], [324.65, -1701.95, 172.95, -1722.15],
                                       [324.60, -1702.00, 172.55, -1722.85], [314.55, -1703.10, 172.00, -1721.25],
                                       [316.95, -1704.15, 172.10, -1721.25], [317.45, -1704.15, 152.05, -1713.45],
                                       [316.6, -1704.35, 152.70, -1714.00], [315.90, -1704.35, 152.70, -1714.00],
                                       [315.90, -1704.35, 152.90, -1714.85], [315.80, -1702.95, 152.75, -1712.70],
                                       [314.25, -1701.75, 152.65, -1712.00]],
            conf_label_min_time: 1542135000,
            conf_label_max_time: 1542150000,
            # conf_label_zero_time: 1542135660,

            conf_label_special_time_intervals: [[75, 100], [110, 180], [180, 215]]

        },
        7443: {
            conf_label_data_file_path: "DOROS_Fill7443_Timber_v02.csv",
            conf_label_offset_time: [1542211800, 1542213059, 1542213602, 1542215346, 1542215382, 1542216031,
                                     1542216740, 1542220016, 1542221586, 1542221877, 1542221936, 1542228161,
                                     1542229709, 1542231642, 1542232040],
            conf_label_offset_values: [[320.40, -1662.80, 161.20, -1682.70], [320.40, -1662.90, 163.10, -1683.55],
                                       [324.55, -1664.90, 162.80, -1683.50], [316.80, -1666.55, 160.10, -1680.95],
                                       [319.05, -1660.20, 160.50, -1680.95], [318.95, -1660.25, 159.00, -1682.15],
                                       [318.10, -1659.50, 158.60, -1681.65], [321.90, -1654.20, 158.05, -1679.75],
                                       [326.45, -1656.15, 157.85, -1679.40], [326.35, -1655.05, 157.75, -1679.50],
                                       [326.80, -1655.05, 139.55, -1673.00], [323.10, -1653.00, 141.10, -1670.75],
                                       [321.85, -1652.40, 142.10, -1671.00], [320.80, -1651.70, 143.30, -1670.40],
                                       [309.75, -1640.15, 143.40, -1670.35]],
            conf_label_min_time: 1542211800,
            conf_label_max_time: 1542233400,
            conf_label_zero_time: 1542211800
        },
        7483: {
            conf_label_data_file_path: "DOROS_Fill7483_Timber.csv",
            conf_label_offset_time: [1543463700, 1543463894, 1543464276, 1543466092, 1543466233],
            conf_label_offset_values: [[340.70, -1641.80, 164.00, -1696.40], [342.00, -1642.80, 164.00, -1696.40],
                                       [341.95, -1642.75, 165.20, -1697.00], [343.30, -1630.35, 164.40, -1692.55],
                                       [343.40, -1630.40, 163.45, -1693.05]],
            conf_label_min_time: 1543463700,
            conf_label_max_time: 1543467300
        }

    },

    "arcBPM": {
        4689: {
            conf_label_data_file_path: "CMS_arc_Fill4689_pos_merged.csv",
            conf_label_offset_time: [1449133477, 1449153809],
            conf_label_offset_values: [[2.2804070000000003, 0.09700449999999997, -4.053462, -1.1241844999999997],
                                       [-2.469151272727275, -1.8924681363636366, -4.418165909090907, -2.768781363636364]],

            conf_label_min_time: 1449134042,
            conf_label_max_time: 1449163025,
            conf_label_zero_time: 1449134042,
            # conf_label_special_time_intervals: [[5, 95]],

            conf_label_y_range: [-10., 10.],
            conf_label_y_diff_range: [-20., 20.],

            conf_label_orbit_drift_plot_range: [-15., 40.],

            # conf_label_deep_studies: True,

            # conf_label_compute_offsets: True,
            # conf_label_time_window_for_offsets: 600,  # s

            # conf_label_compute_offsets_method: "all head-on"
        },
        7442: {
            conf_label_data_file_path: "CMS_arc_Fill7442_pos_merged.csv",
            conf_label_offset_time: [1542135000, 1542139644, 1542141318, 1542145836, 1542148211],
            conf_label_offset_values: [[1.49775, -2.32646, -1.12546, -1.00107],
                                       [-3.6725, -0.861591, 0.768867, -0.114422],
                                       [-3.9657, -0.732589, 1.30359, 0.021505],
                                       [-5.33691, 0.910757, 1.11751, 1.8241],
                                       [-7.0685, 0.110472, 0.856009, 2.08142]],

            conf_label_y_range: [-10., 10.],
            conf_label_y_diff_range: [-20., 20.],
        },
        7443: {
            conf_label_data_file_path: "CMS_arc_Fill7443_pos_merged.csv",
            conf_label_offset_time: [1542212125, 1542215346, 1542216740, 1542220016, 1542221877,
                                     1542228161, 1542229709, 1542231642],
            conf_label_offset_values: [[-1.94542, 3.61686, 0.727639, 1.81203], [-4.88406, 3.41358, 1.69413, 2.33605],
                                       [-5.93267, 3.48517, 2.21978, 2.74673], [-2.95849, 3.88697, 0.733709, 2.21389],
                                       [-3.21834, 4.05577, 0.984082, 2.25003], [-2.35671, 3.04677, 0.951105, 2.69318],
                                       [-2.12264, 2.50401, 1.31455, 2.74764], [-1.7711, 1.40572, 1.7192, 2.68774]],

            conf_label_y_range: [-10., 10.],
            conf_label_y_diff_range: [-20., 20.],

            conf_label_min_time: 1542211800,
            conf_label_max_time: 1542233400,
            conf_label_zero_time: 1542211800
        },
        7483: {
            conf_label_data_file_path: "CMS_arc_Fill7483_pos_merged.csv",
            conf_label_offset_time: [1543463700, 1543466092],
            conf_label_offset_values: [[0.0257093, -0.0330421, 0.0131595, -0.0249084],
                                       [-1.34772, 0.132119, -0.814193, 0.762369]],

            conf_label_y_range: [-10., 10.],
            conf_label_y_diff_range: [-20., 20.],

            conf_label_min_time: 1543463700,
            conf_label_max_time: 1543467300,
            conf_label_zero_time: 1543463700
        },

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