#!/bin/bash

# This script will remove all the temporary files created from running the pipeline
# Putting this part in its own little script makes it easier to deactive if desired
# Last modified on 24-10-23
# WARNING: FILES THAT ARE DELETED THIS WAY CANNOT BE RECOVERED! EDIT THIS FILE AT YOUR OWN RISK!

# Removing all downloaded cutouts
rm -r ZTFdata/sci/*

# Removing all light curves
rm -r ZTFdata/forcephotometry/*.csv

# Removing all results
rm -r results/*
