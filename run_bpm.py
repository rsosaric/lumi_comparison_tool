import optparse
import tools.compare_bpm as compare_bpm
import tools.overall_orbit_drift as orbit_drifts

p = optparse.OptionParser()
usage = "usage: %prog [options] [detector labels]"
p.add_option("-f", "--fill", type="int", help="Year", dest="fill", default=None)
p.add_option("--do_od", action="store_true", help="test class", dest="do_od", default=False)

(options, args) = p.parse_args()

if options.do_od:
    orbit_drifts.get_orbit_drifts(dets_labels=args, fill=options.fill)
else:
    compare_bpm.compare_bpms(dets_labels=args, fill=options.fill)
