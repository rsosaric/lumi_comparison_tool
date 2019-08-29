# -*-*-*-*-*-*-*-*-*-*-*- > > > General settings < < < -*-*-*-*-*-*-*-*-*-*-*-
work_status = 'Preliminary'
# work_status = 'Work in Progress'
experiment = 'CMS'
csv_input_base_dir = 'csv_input_files/'
clean_run = False  # if [True]: all needed files will be produced even if they exist
# -*-*-*-*-*-*-*-*-*-*-*- > > > Plotting settings < < < -*-*-*-*-*-*-*-*-*-*-*-

cms_label_pos_sq = (0.008, 1.02)
cms_label_pos_nsq = (0.008, 1.02)
experiment_font_size = 12
delta_y_pos_sq = 0.08
delta_y_pos_nsq = 0.045

year_energy_label_pos_sq = (0.82, 1.01)
year_energy_label_pos_nsq = (0.87, 1.02)
year_energy_font_size = 12

# canvas sizes
fig_sizes = {
    'sq': (8, 8),
    'nsq': (12, 4)
}

# legends position for histograms mean values, etc (axis converted format x,y: [0->1][0->1])
leg_font_size = 15
pos_leg_1 = (0.75, 0.95)
pos_leg_2 = (0.75, 0.90)
pos_leg_3 = (0.75, 0.85)

# fitting info position in by fill plots
pos_by_fill_fit_info = (0.52, 0.97)

# legends position for ratio vs X(time, date, lumi) plots
leg_vs_plots_pos = 9  # -> 9 for 'center top'
leg_vs_plots_text_s = 12
leg_vs_plots_marker_scale = 7  # -> marker_in_legend/...in_plot scale

# axis configuration:
axis_labelpad = 20
axis_weight = 'bold'

axis_case_size = {
    'sq': 14,
    'nsq': 14
}
axis_thicks_case_size = {
    'sq': 14,
    'nsq': 13
}

ratio_min = 0.95
ratio_max = 1.05
nbins = 300

date_format = "%Y-%m-%d %H:%M"

# -*-*-*-*-*-*-*-*-*-*-*- > > > Luminosity related configuration < < < -*-*-*-*-*-*-*-*-*-*-*-
xLS = 23.31

# -*-*-*-*-*-*-*-*-*-*-*- > > > by year configuration < < < -*-*-*-*-*-*-*-*-*-*-*-
# Always put 5TeV first!!! (less number of fills)
rangeFillYearly = {
    (2015, 5): [4634, 4647],
    (2015, 13): [3829, 4720],
    (2016, 13): [4856, 5575],
    (2017, 13): [5718, 6417],
    (2018, 13): [6570, 7407]
}

energyLabel_list = {
    0: "nanTeV",
    5: "5TeV",
    13: "13TeV"
}

# To be used if no configuration for year is found
nls_default = 20

nls_year = {
    (2015, 5): 15,
    (2015, 13): 10,
    (2016, 13): 50,
    (2017, 13): 50,
    (2018, 13): 50
}
nls_for_lin_year = {
    (2015, 5): 5,
    (2015, 13): 5,
    (2016, 13): 15,
    (2017, 13): 15,
    (2018, 13): 15
}

# -*-*-*-*-*-*-*-*-*-*-*- > > > Input, Output < < < -*-*-*-*-*-*-*-*-*-*-*-
plots_formats = ('pdf', 'png')
default_output_dir = 'plots/'

# -*-*-*-*-*-*-*-*-*-*-*- > > > Data structure and format < < < -*-*-*-*-*-*-*-*-*-*-*-

# TODO: create plot to check this config
allowed_ratio_stdv_factor = 3
remove_outliers_in_lin_by_fill = True
min_number_of_points_req = 5

bad_ratio_to_plot_stdv_factor = 2
lumisensitivity = 0.05                  # Multiplicity of the total luminosity percent of each pint for acceptance.
                                        # Higher means less strict.

# -*-*-*-*-*-*-*-*-*-*-*- > > > all/excluded < < < -*-*-*-*-*-*-*-*-*-*-*-
by_nls_exclusion_threshold = 50.  # percent of excluded points in the nls interval

# -*-*-*-*-*-*-*-*-*-*-*- > > > best/second best < < < -*-*-*-*-*-*-*-*-*-*-*-
# normalization_factor = {
#     2017: {
#         "HFET/DT": 0.995,
#         "HFET/HFOC": 0.981,
#         "HFET/PLTZERO": 1.000,
#         "HFET/PXL": 0.997,
#         "HFOC/DT": 1.014,
#         "HFOC/PXL": 1.017,
#         "PXL/DT": 0.998
#     }
# }

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'