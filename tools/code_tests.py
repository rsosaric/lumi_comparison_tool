import pandas as pd
import numpy as np

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


df1 = pd.DataFrame([['bird', False, True], ['monkey', False, False], ['conejo', False, True], ['perro', True, True]],
                   columns=['animal', 'grande', 'bueno'])

df1["Lo quiero"] = np.where((df1["grande"] == True) & (df1["bueno"] == True), True, False)
print (df1)