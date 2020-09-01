from tools.beammonitor import BPM
import tools.lumi_tools as ltools
import matplotlib.pyplot as plt
import settings as def_setts
import settings_bpm as setts
import tools.plotting_tools as plotting


def get_orbit_drifts(dets_labels: list, fill: int) -> None:
    detectors = []
    output_folder = "plots_bpm/" + str(fill) + "/OD_results/"

    if "Nominal" in dets_labels:
        dets_labels.remove("Nominal")

    nominal_detector = BPM("Nominal", fill=fill, get_only_data=True)
    for det_name in dets_labels:
        detectors.append(BPM(det_name, fill=fill, nominal_data=nominal_detector, get_only_data=True,
                             compute_orbit_drift=True, get_orbit_drift_plots=True))

    plot_orbit_drift_result(detectors, output_folder)

    number_of_dets = len(detectors)
    get_avg_od = False
    if number_of_dets >= 2:
        get_avg_od = True

    ltools.check_and_create_folder(output_folder)
    scan_info = nominal_detector.scans_info_dict
    scan_names = list(scan_info)

    draw_labels_pos_dict = {}
    draw_lines_scans_limits = []
    scan_times = []

    if get_avg_od:
        avg_od = {}
        for scan in scan_names:

            avg_ini_x = 0.0
            avg_mid_x = 0.0
            avg_end_x = 0.0

            avg_ini_y = 0.0
            avg_mid_y = 0.0
            avg_end_y = 0.0

            for det in detectors:
                od_info = det.orbit_drifts_per_scan[scan]
                od_info_x = od_info['OrbitDrifts_X']
                od_info_y = od_info['OrbitDrifts_Y']

                avg_ini_x += od_info_x[0]
                avg_mid_x += od_info_x[1]
                avg_end_x += od_info_x[2]

                avg_ini_y += od_info_y[0]
                avg_mid_y += od_info_y[1]
                avg_end_y += od_info_y[2]

            avg_ini_x = avg_ini_x/number_of_dets
            avg_mid_x = avg_mid_x/number_of_dets
            avg_end_x = avg_end_x/number_of_dets

            avg_ini_y = avg_ini_y/number_of_dets
            avg_mid_y = avg_mid_y/number_of_dets
            avg_end_y = avg_end_y/number_of_dets

            avg_od[scan] = {
                "TimeWindows": detectors[0].orbit_drifts_per_scan[scan]["TimeWindows"],
                "TimeWindows_min": detectors[0].orbit_drifts_per_scan[scan]["TimeWindows_min"],
                'OrbitDrifts_X': [avg_ini_x, avg_mid_x, avg_end_x],
                'OrbitDrifts_Y': [avg_ini_y, avg_mid_y, avg_end_y]
            }
        store_orbit_drift_in_file(orbit_drifts_per_scan=avg_od, output_dir=output_folder, file_name_label="average")

        plot_orbit_drift_result(detectors, output_folder, special_orbit_drifts=avg_od, special_output_name="_avg")


def store_orbit_drift_in_file(orbit_drifts_per_scan, output_dir, file_name_label):
    scans = list(orbit_drifts_per_scan)
    times = []
    od_x = []
    od_y = []

    for scan in scans:
        times.append(orbit_drifts_per_scan[scan]["TimeWindows"])
        od_x.append(orbit_drifts_per_scan[scan]["OrbitDrifts_X"])
        od_y.append(orbit_drifts_per_scan[scan]["OrbitDrifts_Y"])

    with open(output_dir + "OrbitDrifts_" + file_name_label + ".json", 'w') as file:
        file.write("{\n")
        file.write('"Names": ' + str(scans) + ",\n")
        file.write('"TimeWindows": ' + str(times) + ",\n")
        file.write('"OrbitDrifts_X": ' + str(od_x) + ",\n")
        file.write('"OrbitDrifts_Y": ' + str(od_y) + "\n")
        file.write("}")

    ltools.color_print("\n (OUTPUT FILE) Orbit drifts saved in " +
                       output_dir + "OrbitDrifts_" + file_name_label + ".json \n", "yellow")


def plot_orbit_drift_result(detectors, out_dir, special_orbit_drifts=None, special_output_name=""):
    ref_detc = detectors[0]
    x_data_label = ref_detc.col_time_min
    ymin = ref_detc.y_range_orbit_drift[0]
    ymax = ref_detc.y_range_orbit_drift[1]

    title = "LHC Beam Position Monitors"
    plot_file_name = "OrbitDrift_summary"
    ylabel = "Orbit Drift [" + ref_detc.distance_unit + "]"
    # ymin = self.__y_range_orbit_drift[0]
    # ymax = self.__y_range_orbit_drift[1]
    xlabel = "time [min]"
    marker_size = 2.0

    colors_for_detcs = {
        "DOROS": "blue",
        "arcBPM": "green"
    }

    marker_for_detcs = {
        "DOROS": 'o',
        "arcBPM": '^'
    }

    fig, (ax_H, ax_V) = plt.subplots(2, sharex=True, figsize=(12, 7))
    fig.suptitle(title, x=0.24, y=0.98, fontsize=14, fontweight='bold')
    all_plots = (ax_H, ax_V)
    all_plots_y_label = ("H", "V")

    sample_plots_for_legend = []
    sample_labels_for_legend = []

    scans_names = list(ref_detc.orbit_drifts_per_scan)
    y_pos_scans = ymin + abs(ymin * setts.delta_pos_for_scan_labels)
    draw_vertical_line_pos = []
    draw_labels_pos_dict = {}
    for scan in scans_names:
        i_x = ref_detc.orbit_drifts_per_scan[scan]["TimeWindows_min"]

        draw_vertical_line_pos.extend([i_x[0], i_x[2]])
        draw_labels_pos_dict[scan] = [i_x[1], y_pos_scans]

    for det in detectors:
        data_frame = det.data_in_zero_beam_position
        color = colors_for_detcs[det.name]
        marker = marker_for_detcs[det.name]
        scatter_ax_H = ax_H.scatter(data_frame[x_data_label], data_frame[det.col_H_diff],
                                    s=marker_size, c=color, marker=marker)
        ax_V.scatter(data_frame[x_data_label], data_frame[det.col_V_diff],
                     s=marker_size, c=color, marker=marker)

        sample_plots_for_legend.append(scatter_ax_H)
        sample_labels_for_legend.append(det.name)

    if special_orbit_drifts is None:
        for det in detectors:
            for scan in scans_names:
                color = colors_for_detcs[det.name]
                i_x = det.orbit_drifts_per_scan[scan]["TimeWindows_min"]
                i_y_h = det.orbit_drifts_per_scan[scan]["OrbitDrifts_X"]
                i_y_v = det.orbit_drifts_per_scan[scan]["OrbitDrifts_Y"]

                ax_H.plot(i_x, i_y_h, c=color)
                ax_V.plot(i_x, i_y_v, c=color)
    else:
        line_ax_H = None
        for scan in scans_names:
            color = "red"
            i_x = special_orbit_drifts[scan]["TimeWindows_min"]
            i_y_h = special_orbit_drifts[scan]["OrbitDrifts_X"]
            i_y_v = special_orbit_drifts[scan]["OrbitDrifts_Y"]

            line_ax_H = ax_H.plot(i_x, i_y_h, c=color)[0]
            ax_V.plot(i_x, i_y_v, c=color)
        sample_plots_for_legend.append(line_ax_H)
        sample_labels_for_legend.append("average")

    for i in range(0, len(all_plots)):
        ax = all_plots[i]
        ylabel = all_plots_y_label[i] + " [" + ref_detc.distance_unit + "]"

        ax.set_ylabel(ylabel, labelpad=def_setts.axis_labelpad_y, weight=def_setts.axis_weight,
                      size=14)

        ax.set_ylim([ymin, ymax])

        ax.tick_params(axis="x", labelsize=14)
        ax.tick_params(axis="y", labelsize=14)

        ax.xaxis.grid(False)
        ax.yaxis.grid(True)

        # Draw vertical lines
        for x_pos in draw_vertical_line_pos:
            ax.axvline(x=x_pos, linestyle='dashed', alpha=0.3)

        # Draw text from draw_labels_pos_dict -> {"text": [x, y], ... }

        labels_text = list(draw_labels_pos_dict)
        for i_label in labels_text:
            ax.text(draw_labels_pos_dict[i_label][0], draw_labels_pos_dict[i_label][1], i_label,
                    fontweight='bold', alpha=0.5, horizontalalignment='center', verticalalignment='center')

    ax_V.set_xlabel(xlabel, labelpad=def_setts.axis_labelpad_x, weight=def_setts.axis_weight,
                    size=14)

    fig.legend(handles=sample_plots_for_legend,  # The line objects
               labels=sample_labels_for_legend,  # The labels for each line
               loc="center right",  # Position of legend
               borderaxespad=1.0,  # Small spacing around legend box
               markerscale=def_setts.leg_vs_plots_marker_scale,
               fontsize=14
               )

    plt.subplots_adjust(left=0.1, right=0.85, top=0.94, bottom=0.1)
    plotting.save_py_fig_to_file(fig, out_dir + "OD_summary" + special_output_name + ".png",
                                 plot_dpi=def_setts.dpi_png_plots)



