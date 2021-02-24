# -*-*-*-*-*-*-*-*-*-*-*- > > > General settings < < < -*-*-*-*-*-*-*-*-*-*-*-
show_live_plots = False
# work_status = 'Preliminary'
# work_status = 'Work in Progress'
work_status = ""
experiment = 'CMS'
csv_input_base_dir = 'csv_input_files/'
clean_run = True  # Currently only used in best data comparison if [True]: all needed files will be produced even if they exist.
get_extra_plots = False
get_by_fill_by_run_plots = False
get_normalized_plots = False
# skipping lines when reading input file
rows_to_skip_in_csv = 1
rows_to_skip_at_the_end = 5
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
pos_by_fill_fit_info = (0.40, 0.96)
pos_by_fill_fit_chi2_info = (0.60, 0.90)

# legends position for ratio vs X(time, date, lumi) plots
leg_vs_plots_pos = 9  # -> 9 for 'center top', 1: 'upper right'
leg_vs_plots_text_s = 12
leg_vs_plots_marker_scale = 7  # -> marker_in_legend/...in_plot scale

# axis configuration:
axis_labelpad = 10
axis_labelpad_x = 10
axis_labelpad_y = 8
axis_weight = 'bold'

axis_case_size = {
    'sq': 14,
    'nsq': 14
}
axis_thicks_case_size = {
    'sq': 14,
    'nsq': 13
}

# ratio_min = 0.995
# ratio_max = 1.005

# ratio_min = 0.98
# ratio_max = 1.02

ratio_min = 0.95
ratio_max = 1.05

# ratio_min = 0.90
# ratio_max = 1.10

# ratio_min = 0.80
# ratio_max = 1.20

# ratio_min = 0.60
# ratio_max = 1.40


ratio_interval_reduction_factor = 0.0
ratio_min_reduced = ratio_min + ratio_interval_reduction_factor
ratio_max_reduced = ratio_max - ratio_interval_reduction_factor
nbins = 200
nbins_by_run = 150
nbins_by_fill = 50
nbins_linearity = 150
nbins_sbil_histo = 50

date_format = "%Y-%m-%d %H:%M"

# -*-*-*-*-*-*-*-*-*-*-*- > > > Luminosity related configuration < < < -*-*-*-*-*-*-*-*-*-*-*-
xLS = 23.31

# -*-*-*-*-*-*-*-*-*-*-*- > > > by year configuration < < < -*-*-*-*-*-*-*-*-*-*-*-
# Always put 5TeV first!!! (less number of fills)
rangeFillYearly = {
    (2015, 5): [4634, 4647],
    (2015, 13): [3829, 4720],
    (2016, 13): [4856, 5575],
    (2017, 5): [6380, 6398],
    (2017, 13): [5698, 6417],
    (2018, 5): [7427, 7492],
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
    (2017, 5): 15,
    (2017, 13): 50,
    (2018, 5): 20,
    (2018, 13): 50
}
nls_for_lin_year = {
    (2015, 5): 5,
    (2015, 13): 5,
    (2016, 13): 5,
    (2017, 5): 5,
    (2017, 13): 15,
    (2018, 13): 15
}

points_in_summary_lin_year = {
    (2015, 5): 4,
    (2015, 13): 4,
    (2016, 13): 5,
    (2017, 5): 4,
    (2017, 13): 5,
    (2018, 13): 5
}

# limit for the slope values in the linearity analysis: (min slope, max slope) [hz/muB]
limits_linearity_slopes_year = {
    (2015, 5): (-0.05, 0.05),
    (2015, 13): (-0.05, 0.05),
    (2016, 13): (-0.01, 0.01),
    (2017, 13): (-0.01, 0.01),
    (2017, 5): (-0.01, 0.01),
    (2018, 13): (-0.05, 0.05)
}

# sbil range for plotting [ub]
sbil_min = 0.0
sbil_max = 10.0

# require min number of points to fit
min_num_data_to_fit = 5

# -*-*-*-*-*-*-*-*-*-*-*- > > > Input, Output < < < -*-*-*-*-*-*-*-*-*-*-*-
# plots_formats = ('pdf', 'png')
plots_formats = 'png'
# plots_formats = 'pdf'
dpi_png_plots = 300
default_output_dir = 'plots/'
linearity_analysis_output_dir = 'linearity_analysis/'

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


# -*-*-*-*-*-*-*-*-*-*-*- > > > tests options < < < -*-*-*-*-*-*-*-*-*-*-*-
# nls_list_to_test = (5, 10, 15, 50, 100, 150, 200, 400)
# nls_list_to_test = (5, 10, 15, 50, 70, 100, 150)
nls_list_to_test = (1, 5, 10, 15, 30, 50)

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