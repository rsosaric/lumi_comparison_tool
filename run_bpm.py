import optparse
import tools.compare_bpm as compare_bpm

p = optparse.OptionParser()
usage = "usage: %prog [options] [detector labels]"
p.add_option("-f", "--fill", type="int", help="Year", dest="fill", default=None)

(options, args) = p.parse_args()

compare_bpm.compare_bpms(dets_labels=args, fill=options.fill)
