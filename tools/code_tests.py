#import pandas as pd
#import numpy as np

#import tools.lumi_tools as ltools

import beammonitor as bmonitor

# df1 = pd.DataFrame([['bird', 'polly'], ['monkey', 'george'], ['conejo', 'lucas'], ['conejo1', 'lucas1']],
#                    columns=['animal', 'name'])
#
# df2 = pd.DataFrame([['bird', 'polly'], ['monkey', 'george1'], ['conejo1', 'lucas1']],
#                    columns=['animal', 'name'])

#print(df1[~df1['animal'].isin(df2['animal'])])

# print(df2)
# df1["animal not in 2"] = ~df1['animal'].isin(df2['animal'])
# print (df1)


# df = pd.DataFrame({"name": ['Alfred', 'Batman', 'Catwoman'],
#                    "toy": [np.nan, 'Batmobile', 'Bullwhip'],
#                    "born": [pd.NaT, pd.Timestamp("1940-04-25"),
#                             pd.NaT]})
# print (df)
#
# print (df.dropna(subset=['toy']))


# df1 = pd.DataFrame([['bird', False, True], ['monkey', False, False], ['conejo', False, True], ['perro', True, True]],
#                    columns=['animal', 'grande', 'bueno'])
#
# df1["Lo quiero"] = np.where((df1["grande"] == True) & (df1["bueno"] == True), True, False)
# print (df1)


# l1 = [0.1, 0.2, 0.3]
# l2 = [0.00002, 0.00004, 0.0003]
# l1 = np.array(l1)
# l2 = np.array(l2)
#
# print(ltools.get_total_weights(l1, l2))


test = bmonitor.bpm_col_b1hL


