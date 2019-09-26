*********************** INSTALLATION ***********************

Suggested: install miniconda or anaconda. Then use the provided file for create a working conda environment: 

conda env create -f conda_env.yml

Then every time you start a new console you have to load the environment:

source activate lumi_tools

** PYTHON libraries needed **

- numpy=1.16.3
- numpy-base=1.16.3
- numpydoc=0.9.1
- pandas=0.24.2
- patsy=0.5.1
- pip=19.1.1
- python=3.6.8
- scipy=1.2.1
- seaborn=0.9.0
- statsmodels=0.10.0
- tornado=6.0.2

*********************** INSTRUCTIONS ***********************

For running the comparison always use the file: run.py
Example: python run.py -y 2015 pcc dt (python run.py [opts] det1 det2 ...)

The detector label must be the same that the one in the .csv file name -->> pcc.csv, dt.csv, etc.

!!Important: the input folder for .csv files is automatically taken as csv_input_files/year if [-y 'year'] is used.
In any other case input folder must be set using [-i 'folder'].
In the example the input folder would be: csv_input_files/2015/.

For modifications in plots and analysis cuts go to settings.py
For minor modification of the analysis (stability, linearity) go first to lumianalysis.py

*********************** IMPLEMENTED OPTIONS ***********************

use the help: python run.py --help

*********************** using recorded or delivered lumi? ***********************

By default the code uses recorded luminosity from the detectors. If you want to use delivered please use: --lumi_type del

*********************** Analysis of the exclusion of bad zones (runs, fills) ***********************
For this it are needed two files for each detector: det.csv(selected data) and det_all.csv (all data)

For example: python run.py -a -y 2017 hfet pcc

Is also possible for 3 detectors: python run.py -a -y 2017 hfet pcc hfoc


*********************** Best/second best analysis ***********************

example for 2017: python run.py -p -y 2017 best second

** Where best, second should be in this path: csv_input_files/2017/ as best.csv and second.csv

*********************** Ramses(in principle can be used for any other cross calibration) ***********************

For this you need two kind of files: ramses_channels_raw_data.csv (containing time and dose) and csv file from a reference detector: hfet_all.csv

Here one example:  python run.py --ramses_crosscal csv_input_files/ramses_2017_raw/ramses_dataframe_1.csv csv_input_files/ramses_2017_raw/ramses_dataframe_2.csv csv_input_files/2017/hfet_all.csv