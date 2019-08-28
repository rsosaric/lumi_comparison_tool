from tools.lumianalysis import LAnalysis as LA
import numpy as np


def exclusion_fill(dets_file_labels: list, input_dir: str, mixed_data=False, run_stddev_test=False,
                   c_years=False) -> None:

    print('Executing exclusion mode ...\n')
    n_files = len(dets_file_labels)
    if n_files < 3:
        raise AssertionError("Number of detectors must be equal or higher than 3.")
    det_read = []
    years_and_dir = input_dir.split('/')
    year = years_and_dir[1]
    intersection = {}
    bad_runs_in_detectors = {}

    for k in range(0, n_files, 1):
        intersection[k] = []

    counts = 0
    datatype = type(1)

    for i in range(0, n_files - 1, 1):
        for j in range(i + 1, n_files, 1):
            input_path = 'plots/' + year + '/' + dets_file_labels[i] + '-' + dets_file_labels[j] + \
                         '/txt/Bad_runs.txt'
            try:
                det_read.append(np.loadtxt(input_path, dtype=datatype))
            except IOError as errIO:
                print('\n Missing ' + dets_file_labels[i] + '-' + dets_file_labels[j] + ' bad runs \n')
                det_missing = [dets_file_labels[i], dets_file_labels[j]]
                LA(det_missing, input_dir=input_dir, mixed_data=mixed_data,
                          run_stddev_test=run_stddev_test, c_years=c_years, exclusion=False)
                det_read.append(np.loadtxt(input_path, dtype=datatype))

            intersection[i].append(counts)
            intersection[j].append(counts)
            counts = counts + 1

    for i in range(0, n_files, 1):
        bad_runs_in_detectors[i] = np.intersect1d(det_read[intersection[i][0]], det_read[intersection[i][1]])
        for j in range(2, len(intersection[i]), 1):
            bad_runs_in_detectors[i] = np.intersect1d(bad_runs_in_detectors[i], det_read[intersection[i][j]])

    filePath = 'plots/'
    txtfileName = 'Bad_runs_in_detectors'
    fileout = open(filePath + txtfileName + ".txt", "w+")
    for i in range(0, n_files, 1):
        fileout.write(dets_file_labels[i] + ': ' + str(bad_runs_in_detectors[i]) + '\n')
    fileout.close()