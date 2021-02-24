import numpy as np
import pandas as pd
import optparse

p = optparse.OptionParser()
usage = "usage: %prog [options] [detector labels]"

# p.add_option("-o", "--outdir", type="string", help="Path to output Dir", dest="outDir", default=".")
p.add_option("-c", "--correction", type="float", help="Correction factor", dest="corr_factor", default=None)

(options, args) = p.parse_args()

input_file = args[0]
if options.corr_factor is None:
    raise AssertionError("Correction factor is needed (use -c option)")

corr_factor = options.corr_factor

delivered_lumi_col = "delivered(hz/ub)"
recorded_lumi_col = "recorded(hz/ub)"

rows_to_skip_at_the_end = 5
rows_to_skip_in_csv = 1
data_file_pd_all = pd.read_csv(input_file, index_col=False, skiprows=rows_to_skip_in_csv,
                               skipfooter=rows_to_skip_at_the_end, engine='python')

print("Correcting " + delivered_lumi_col + " and " + recorded_lumi_col + " columns by 1/" + str(corr_factor))
data_file_pd_all[delivered_lumi_col] = data_file_pd_all[delivered_lumi_col]/corr_factor
data_file_pd_all[recorded_lumi_col] = data_file_pd_all[recorded_lumi_col]/corr_factor

data_file_pd_all.to_csv(input_file.replace(".csv", "_norm.csv"), index=False)


