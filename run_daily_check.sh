#!/bin/bash

#This shell script is designed to run the entire real-time pipeline on a daily basis
#Last modified on 24-10-23

#set number of sources lists
nr_lists=28

#Use current mjd to loop through the lists on a daily basis
cur_date=$(( $(date +%s) / 86400 + 40587))
list_today=$((cur_date - $(expr $cur_date / $nr_lists)*nr_lists))
list_to_use="sourcelists/list_${list_today}.csv"

#Make sure the correct environment is activated
source ~/miniconda3/etc/profile.d/conda.sh
conda activate real-time_search

#Set the lc location
lc_loc="lcs/real_time/list_${list_today}"

#Run the lc_generator
python3 lc_generator.py $list_to_use $lc_loc

#Run the Binning_program
python3 Binning_program.py $lc_loc results config.csv

#Run the auto_check
python3 auto_check.py daily
#put in manual when running in manual mode and only emailing the 1st person on the list

#Remove temporary files
./scrapper.sh
