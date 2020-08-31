'''
To be run in lxplus!
Please first do:
export _JAVA_OPTIONS="-Xss2m"
source /cvmfs/sft.cern.ch/lcg/views/LCG_85swan3/x86_64-slc6-gcc49-opt/setup.sh
'''

import pytimber
import datetime
import csv
import pandas as pd
import numpy as np
import optparse

fills_timestamp_ranges = {
    "4689": [1449134042, 1449163025],
    "4266": [1440462850, 1440480981],
}

allowed_fills = list(fills_timestamp_ranges)

p = optparse.OptionParser()
usage = "usage: %prog [options] [detector labels]"
p.add_option("-f", "--fill", type="str", help="Year", dest="fill", default=None)

(options, args) = p.parse_args()

fill_number = options.fill
if fill_number not in allowed_fills:
    raise AssertionError("Fill " + fill_number +  " not in allowed fills: " + str(allowed_fills))

output_filename = fill_number + "_DOROS_raw_data_from_Timber.csv"

start_ts = fills_timestamp_ranges[fill_number][0]
end_ts = fills_timestamp_ranges[fill_number][1]
start = datetime.datetime.fromtimestamp(int(start_ts)).strftime('%Y-%m-%d %H:%M:%S')
end = datetime.datetime.fromtimestamp(int(end_ts)).strftime('%Y-%m-%d %H:%M:%S')

requiered_data = ['LHC.BPM.1L5.B1_DOROS:POS_H',
                  'LHC.BPM.1L5.B1_DOROS:POS_V',
                  'LHC.BPM.1L5.B2_DOROS:POS_H',
                  'LHC.BPM.1L5.B2_DOROS:POS_V',
                  'LHC.BPM.1R5.B1_DOROS:POS_H',
                  'LHC.BPM.1R5.B1_DOROS:POS_V',
                  'LHC.BPM.1R5.B2_DOROS:POS_H',
                  'LHC.BPM.1R5.B2_DOROS:POS_V'
                  ]
column_order = ["sec"] + requiered_data

print("Getting DOROS data for fill " + fill_number + ": ")
print("from " + start + " (timestamp: " + str(start_ts) + ")")
print("until " + end + " (timestamp: " + str(end_ts) + ")")

db = pytimber.LoggingDB()

print("Getting data from Timber ... ")
data = db.get(requiered_data, start, end)

print ("Data collected from Timber:")
for i_data_label in requiered_data:
    i_time, i_values = data[i_data_label]
    print(i_data_label + " : size = " + str(len(i_values)))
    
# Export dict to file
with open("raw_dict_" + output_filename,'w') as f:
    w = csv.writer(f)
    w.writerows(data.items())
print("Raw dict saved in " + "raw_dict_" + output_filename)
    
print("Merging data ... ")
# convert data into pandas df
data_df_dict = {}
for i_name in requiered_data:
    i_time, i_values = data[i_name]
    data_df_dict[i_name] = pd.DataFrame({i_name: i_values, "sec": np.array(i_time).astype(int)})
    
merged_df = pd.merge(data_df_dict[requiered_data[0]], data_df_dict[requiered_data[1]], on='sec')
print(0, requiered_data[0])
print(1, requiered_data[1])

for i in range(2, len(requiered_data)):
    print(i, requiered_data[i])
    merged_df = pd.merge(merged_df, data_df_dict[requiered_data[i]], on='sec')
                    
# Reordering column positions
merged_df = merged_df[column_order]

merged_df.dropna(inplace=True)
print("Size of merged data: " + str(len(merged_df)))

merged_df.to_csv(output_filename, index=False)
print("Data saved into " + output_filename)

print("Done :)")
