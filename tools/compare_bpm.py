from tools.beammonitor import BPM
# import pandas as pd


def compare_bpms(dets_labels: list, fill: int) -> None:
    detectors = []
    get_only_nominal_data = True

    if "Nominal" not in dets_labels:
        dets_labels.insert(0, "Nominal")
    elif "Nominal" in dets_labels and dets_labels[0] != "Nominal":
        dets_labels.remove("Nominal")
        dets_labels.insert(0, "Nominal")
        get_only_nominal_data = False
    elif "Nominal" in dets_labels:
        get_only_nominal_data = False

    n_det = 0
    for det_name in dets_labels:
        # producing Nominal
        if n_det == 0:
            detectors.append(BPM(det_name, fill=fill, get_only_data=get_only_nominal_data))
        else:
            detectors.append(BPM(det_name, fill=fill, nominal_data=detectors[0]))
        n_det += 1
