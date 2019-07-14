from tools.luminometer import Luminometer as L
from tools.detectorsratio import DetectorsRatio as Ratios
from tools.detectorsratio import MultipleDetectorsRatio as MRatios
import optparse
import settings as setts
import tools.lumi_tools as ltools
from tools.linearityanalysis import LinearityAnalysis

# TODO: develop a proper input option interface here
p = optparse.OptionParser()
usage = "usage: %prog [options] arg1 arg2 ..."
# p.add_option("-o", "--outdir", type="string",help="Path to output Dir", dest="outDir", default=".")
p.add_option("-i", "--indir", type="string",help="Path to input Dir", dest="inDir", default=None)
p.add_option("-f", "--plot_format", type="string", help="output format for plots",
             dest="plot_format", default="")

p.add_option("-y", "--year", type="string", help="Year", dest="year", default=None)



(options, args) = p.parse_args()

n_files=len(args)
plot_format = options.plot_format
if plot_format == "":
    plot_format = setts.plots_formats

if n_files == 2:
    print("******** 2 detector comparison choose!! ******")
    det1_label = args[0]
    det2_label = args[1]
    if options.inDir:
        base_input_path = options.inDir
    elif options.year:
        base_input_path = setts.csv_input_base_dir + options.year + '/'
        try:
            detc1 = L(det1_label, base_input_path + det1_label + ".csv")
            detc2 = L(det2_label, base_input_path + det2_label + ".csv")
        except IOError as errIO:
            print(errIO)
            print('Please check if default input folder is correctly created: ' + setts.csv_input_base_dir)
            print('Also check that your .csv file is in the correct year folder: ' + base_input_path)
            raise
    else:
        raise IOError('Please specify input folder or year for the input .csv files')

    ratios12 = Ratios(detc1, detc2)

    # Fill plots
    ratios12.plot_nls_ratio_hist()

    # Save plots
    ratios12.save_plots()


else:
    raise IOError('No detectors names entered')
