import pandas as pd
import numpy as np
import tools.lumi_tools as ltools

# df1 = pd.DataFrame({'time': [1, 2, 3, 4, 6, 10],
#                     'ratio_1': [0.1, 0.3, 4.5, 5.5, 2.2, 1.0],
#                     'l_pcc': [0.1123, 0.3123123, 4.512313, 5.5123123, 2.212313, 1.012312],
#                     'r_pcc': [0.1123, 0.3123123, 4.512313, 5.5123123, 2.212313, 1.012312]})
#
# df2 = pd.DataFrame({'time': [1, 2, 3, 5, 6, 9],
#                     'ratio_2': [0.2, 0.5, 4.7, 5.9, 2.3, 1.8],
#                     'l_pcc': [0.1123, 0.3123123, 4.512313, 5.5123123, 2.212313, 3.21212],
#                     'r_pcc': [0.1123, 0.3123123, 4.512313, 5.5123123, 2.212313, 1.012312]})
# df3 = pd.DataFrame({'time': [1.1, 2.3, 3.4, 5.5, 6.6, 9.8],
#                     'ratio_3': [np.nan, 0.5, 4.7, 5.9, 2.3, 1.8],
#                     'l_pcc': [np.nan, 0.3123123, 4.512313, np.nan, 2.212313, 3.21212]})
#
# result_outer = pd.merge(df1, df2, on='time', how='outer', sort=True).merge(df3, on='time', how='outer', sort=True)
#
# print('outer')
# print(result_outer)
#
# nan_vals = result_outer['l_pcc_x'].isnull().values
# full_vals = []
# for iii in range(0,len(result_outer)):
#     if nan_vals[iii]:
#         full_vals.append(result_outer['l_pcc_y'][iii])
#     else:
#         full_vals.append(result_outer['l_pcc_x'][iii])
#
# # for col in result_outer['lumi_pcc_x']:
# #     print(col, type(col))
# #     if col is np.nan:
# #         print('NaN!!!!!!!!!!!!!!!!!!!!!!!!!')
#
# ltools.check_and_clean_after_merging(result_outer)
# # result_outer['lumi_pcc'] = np.array(full_vals)
# print(result_outer)


df1 = pd.DataFrame({'time': [1, 2, 3, 4, 6, 10, 11],
                    'fill': [1122, 1122, 3123, 4123, 4123, 9034, 9035],
                    'ratio_1': [0.1, 0.3, 4.5, 5.5, 2.2, 1.0, 1.1],
                    'l_pcc': [0.1123, 0.3123123, 4.512313, 5.5123123, 2.212313, 1.012312, 1.00002],
                    'r_pcc': [0.1123, 0.3123123, 4.512313, 5.5123123, 2.212313, 1.012312, 1.000002]})

nls = pd.DataFrame({'fill': [2123, 1122, 3123, 4123, 6343, 9034],
                    'nls': [5, 3, 4, 6, 10, 20]})


result_outer = pd.merge(df1, nls, on='fill', how='left', sort=True)

print(result_outer)