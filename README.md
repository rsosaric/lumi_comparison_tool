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

"-i": Path to input Dir" (default=None)

"-a": all data vs. selected analysis" (default=False)

"-y": Year (default=None) !This is only used for input folder, year is determined using the data in the .csv

"-l": linearity analysis" (default=False)