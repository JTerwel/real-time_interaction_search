#!/bin/bash

#This shell script is designed to run the entire pipeline in manual mode
#Last modified on 13-11-23

list_to_use="sourcelists/list_0.csv"

#Make sure the correct environment is activated
source ~/miniconda3/etc/profile.d/conda.sh
conda activate real-time_search

#Set the lc location
lc_loc='lcs/real_time/list_0'

#Run the lc_generator
python3 lc_generator.py $list_to_use $lc_loc 32 32

#Run the Binning_program
python3 Binning_program.py $lc_loc results config.csv

#Run the auto_check
python3 auto_check.py manual 4 4 28
#put in manual when running in manual mode and only emailing the 1st person on the list

#Remove temporary files
#./scrapper.sh
