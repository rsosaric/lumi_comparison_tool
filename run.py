import optparse
import settings as setts
from tools.lumianalysis import LAnalysis as Lumi
from tools.bestdata import BestDataAnalysis

# TODO: develop a proper input option interface here
p = optparse.OptionParser()
usage = "usage: %prog [options] [detector labels]"
# p.add_option("-o", "--outdir", type="string",help="Path to output Dir", dest="outDir", default=".")
p.add_option("-i", "--indir", type="string", help="Path to input Dir", dest="inDir", default=None)
p.add_option("-a", "--all", action="store_true", help="all data vs. selected analysis",
             dest="all", default=False)
p.add_option("-m", "--mixed", action="store_true", help="for handling .csv with more than one detector",
             dest="mixed", default=False)
p.add_option("-y", "--year", type="string", help="Year", dest="year", default=None)
p.add_option("-c", "--combined_years", type="string", help="Years - Introduce separated by comma", dest="years", default=None)
p.add_option("-l", "--lin", action="store_true", help="linearity analysis", dest="lin_an", default=False)
p.add_option("-t", "--test", action="store_true", help="test class", dest="test", default=False)
p.add_option("-p", "--physics", action="store_true", help="physics selected data analysis", dest="physics", default=False)
p.add_option("-e", "--exclusion", action="store_true", help="Find bad detector behavior", dest="exclusion", default=False)

(options, args) = p.parse_args()

several_years = False

if options.inDir:
    base_input_path = options.inDir
elif options.year:
    base_input_path = setts.csv_input_base_dir + options.year + '/'
elif options.years:
    base_input_path = setts.csv_input_base_dir + ',' + options.years
    several_years = True
else:
    raise IOError('Please specify input folder or year for the input .csv files')

if options.test:
    print("test is being done")

if options.physics:
    physics_analysis = BestDataAnalysis(dets_file_labels=args, input_dir=base_input_path, c_years=several_years)
else:
    lumi_analysis = Lumi(dets_file_labels=args, input_dir=base_input_path, run_linearity_analysis=options.lin_an,
                         mixed_data=options.mixed, run_stddev_test=options.test, c_years=several_years,
                         all_and_excluded_analysis=options.all, exclusion=options.exclusion)

