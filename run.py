import optparse
import settings as setts
from tools.lumianalysis import LAnalysis as Lumi

# TODO: develop a proper input option interface here
p = optparse.OptionParser()
usage = "usage: %prog [options] arg1 arg2 ..."
# p.add_option("-o", "--outdir", type="string",help="Path to output Dir", dest="outDir", default=".")
p.add_option("-i", "--indir", type="string", help="Path to input Dir", dest="inDir", default=None)
p.add_option("-a", "--all", action="store_true", help="all data vs. selected analysis",
             dest="all", default=False)
p.add_option("-m", "--mixed", action="store_true", help="for handling .csv with more than one detector",
             dest="mixed", default=False)
p.add_option("-y", "--year", type="string", help="Year", dest="year", default=None)
p.add_option("-l", "--lin", action="store_true", help="linearity analysis", dest="lin_an", default=False)

(options, args) = p.parse_args()

if options.inDir:
    base_input_path = options.inDir
elif options.year:
    base_input_path = setts.csv_input_base_dir + options.year + '/'
else:
    raise IOError('Please specify input folder or year for the input .csv files')

lumi_analysis = Lumi(dets_file_labels=args, input_dir=base_input_path, run_linearity_analysis=options.lin_an,
                     mixed_data=options.mixed)
