from tools.luminometer import Luminometer as L
from tools.detectorsratio import DetectorsRatio as Ratios
from tools.detectorsratio import MultipleDetectorsRatio as MRatios
from tools.linearityanalysis import LinearityAnalysis

# TODO: develop a proper input option interface here
year = 2015
# detc1 = L("hfoc", "csv_input_files/" + str(year) + "/hfoc.csv")
detc1 = L("dt", "csv_input_files/" + str(year) + "/dt.csv")
detc2 = L("ram", "csv_input_files/" + str(year) + "/ram.csv")
detc3 = L("pcc", "csv_input_files/" + str(year) + "/pcc.csv")

detc2all = L("ram", "csv_input_files/" + str(year) + "/ram"+'_all'+".csv")
detc1all = L("hfoc", "csv_input_files/" + str(year) + "/hfoc"+'_all'+".csv")
detc3all = L("pcc", "csv_input_files/" + str(year) + "/pcc"+'_all'+".csv")

ratios13 = Ratios(detc1, detc3)


#print(detc1.data)
#print(detc2.data)
#print(ratios12.common_data)
#print(ratios12.data_exclusion_percent)

# data analysis
# ratios12.plot_nls_ratio_hist()
# # ratios12.plot_ratio_hist()
# # # ratios12.plot_ratio_hist_weighted()
# ratios12.plot_nls_ratio_hist_weighted()
# # ratios12.plot_nls_ratio_vs_date()
# ratios12.plot_nls_ratio_vs_lumi2()
# #
# ratios12.save_plots()



# Linearity
# lin_analysis = LinearityAnalysis(ratios13)
# lin_analysis.plot_hist_sbil()
# lin_analysis.plot_ratio_vs_all_data_sbil()
# lin_analysis.save_plots()

ratios123 = MRatios(detc1, detc2, detc3, full_data=[detc1all, detc2all, detc3all])
ratios123.plot_nls_ratios_vs_date()
ratios123.plot_ratios_vs_date()
ratios123.plot_ratios_vs_lumi3()
ratios123.plot_nls_ratios_vs_lumi3()
#
ratios123.save_plots()
