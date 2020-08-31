import csv
import pandas as pd
import numpy as np

input_folder = "Raw/"

B1_L_name = "B1_L"
B1_R_name = "B1_R"
B2_L_name = "B2_L"
B2_R_name = "B2_R"

cols_file_name_label={
    B1_L_name: "b1L",
    B1_R_name: "b1R",
    B2_L_name: "b2L",
    B2_R_name: "b2R"
}

cols_to_merge = list(cols_file_name_label)
number_of_files = 8

output_filename = "CMS_arc_Fill4689_pos_merged.csv"

# cols to read 
to_read_time = "Timestamp"
to_read_H = ' IP5-H:Pos'
to_read_V = ' IP5-V:Pos'

cols_to_read = [to_read_time, to_read_H, to_read_V]
exclude_cols = [' IP5-H:Residuals', ' IP5-V:Residuals']

all_data = pd.DataFrame()

for i_col in cols_to_merge:
    
    all_data_col = pd.DataFrame()
    modified_H = i_col.replace("_", "_H_")
    modified_V = i_col.replace("_", "_V_")
    rename_dict = {
        to_read_H: modified_H,
        to_read_V: modified_V
    }
    counter = 0
    for i_split in range(1, number_of_files+1):
        input_file_name = input_folder + str(i_split) + "_" +  cols_file_name_label[i_col] + ".csv"
        in_data = pd.read_csv(input_file_name, usecols = cols_to_read,  skiprows=[1])
        
        print(input_file_name)

        counter+=len(in_data)
        if len(all_data_col) == 0:
            all_data_col = in_data
        else:
            all_data_col = all_data_col.append(in_data, ignore_index=True)

    all_data_col[to_read_time] = all_data_col[to_read_time] / 1000.
    all_data_col = all_data_col.astype({to_read_time: int})
    all_data_col.rename(columns=rename_dict, inplace=True)

    # all_data_col.drop_duplicates(inplace=True, subset=to_read_time, keep='first')

    print(" ------------------- Total data for " + i_col + " ------------------- ")
    print(all_data_col)
    print(counter)

    if len(all_data) == 0:
        all_data = all_data_col
    else:
        all_data = pd.merge(all_data, all_data_col, on=to_read_time)
    
print("\n ------------------- Total merged data ------------------- ")
print(all_data)
all_data.to_csv(output_filename, index=False)


    
# merged_df = pd.merge(data_df_dict[requiered_data[0]], data_df_dict[requiered_data[1]], on='sec')
# print(0, requiered_data[0])
# print(1, requiered_data[1])

# for i in range(2, len(requiered_data)):
#     print(i, requiered_data[i])
#     merged_df = pd.merge(merged_df, data_df_dict[requiered_data[i]], on='sec')
                    
# # Reordering column positions
# merged_df = merged_df[column_order]

# merged_df.dropna(inplace=True)
# print("Size of merged data: " + str(len(merged_df)))

# merged_df.to_csv(output_filename, index=False)
# print("Data saved into " + output_filename)

print("Done :)")
