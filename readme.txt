##############################################
# Real-time Late-time CSM interaction search #
##############################################

This folder contains the necessary code and files to automatically run a real-time search for late-time interaction signatures in Type Ia SNe in ZTF.
This is done by running lc_generator.py, Binning_program.py, auto_check.py, and scrapper.py in order.
This file is meant as a guide to get it working, containing an inventory, installation and use guide, mongo configuration details, and known issues with solutions.

Author: 	Jacco Terwel
email:		terwelj@tcd.ie
Last modified:	16-11-23


=====================
|| Folder contents ||
=====================
readme.txt		(text file):		This file
config.csv		(csv file):		Configuration file for Binning_program.py for automated use
manual_config.csv	(csv file):		Configuration file for Binning_program.py for manual use
blacklist.csv		(csv file):		Contains objects that should not be run for some specific reason
testlist.csv		(csv file):		List of 16 ZTF objects used to test the programs
ZTFdata			(folder): 		Contains the output generated by the fpbot and ztfquery packages (image cutouts and light curves)
sourceslists		(folder):		Contains the lists of objects, one of these should be used every day, circulating through them all
obs_files		(folder):		Contains the other csv files needed for Binning_program.py to work (ref_im_data, SN_host_data & locs, and zp_thresholds_quadID)
mongo			(folder):		Contains everything to run mongodb locally when installing it normally would lead to security concerns.
results			(folder):		Contains the output of the python scripts that isn't saved in ZTFdata.
lcs			(folder):		Contains the light curves generated with fpbot (Different subfolders for different projects)
backup			(folder):		Results are saved here if the email fails to prevent them from being erased
lc_generator.py		(python script):	Reads in a list of objects and checks for new images. If new images are found, cutouts are downloaded for forced photometry, and a light curve is created.
						Usage: python3 lc_generator.py <object list> <lc directory>
							or python3 lc_generator.py <object list> <lc directory> <nr. processes download> <nr. processes fit>
						(Tresholds default = 32, 32)
Binning_program.py	(python script):	Runs the binning program on a given folder of lcs, storing the results in a given location, and using a specific configuration.
						Usage: python3 Binning_program.py <lc loc> <save loc> <config file>
auto_check.py		(python script):	Checks the output of the binning program for currently active rebrightenings and sends the results by email. Only emails 1st person listed in manual mode.
						Usage: python3 auto_check.py <mode> or python3 auto_check.py <mode> <nr. of dets to pass> <nr. of final bin dets to pass> <nr. of days considered recent>
						(Tresholds default = 4, 4, 28)
scrapper.sh		(bash script):		Deletes all downloaded image cutouts, light curves, and other output that was created by the python scripts.
						Usage: ./scrapper.sh
run_daily_check.sh	(bash script):		Automated script to check one of the lists in sourcelists each day and run all python programs in order.
						Usage: ./run_dily_check.sh
run_manually.sh		(bash script):		Script to run the programs on a specific set of objects, to be used for manual checks. May be adapted to specific cases without breaking the automated one
						Usage: ./run_manually.sh
list_creator_updater.ipynb (jupyter notebook):	Notebook for creating / updatings the source lists in a consistent manner
check_lcs.ipynb		(jupyter notebook):	Simple notebook to check light curves and binning program results. Has some basic functions for reading the files and plotting light curves and bins


========================
|| Installation guide ||
========================
- Install python 3.10 or newer
- Install ztfquery: (pip install ztfquery)
- Create the $ZTFDATA global variable (e.g. in ~/.bash_profile add: export ZTFDATA="<path to this folder>/ZTFdata/"  NOTE: DO NOT FORGET THE FINAL / OR FPBOT WILL CRASH WHEN MAKING THE LIGHT CURVES)
- Install lmfit: (pip install lmfit)
- Install fpbot: (Follow steps 1-4 on https://github.com/simeonreusch/fpbot, make sure to have the latest version or fpbot will download entire images instead of cutouts)
- Get mongo in a container and set it up properly: (steps are dependent on containerization, example assumes apptainer)
	- apptainer pull docker://mongo
	- Add a custom config file for the mongodb (see https://www.mongodb.com/docs/manual/administration/configuration/) Define port, dbpath and path properly. The rest can stay as in the example
	  This circumvents an error which occurs when using the default port. See below for another example of the contents of the config file
	- apptainer run mongo_latest.sif --config <path to config file>
	- Check that mongodb is indeed running by using pgrep mongo (This should return one or more numbers. If not, it isn't running)
	- Make sure that fpbot uses the same port as is specified in the mongo config file (Change line 32 of fpbot/database.py)
- Make sure all folders listed above are present
- Start a python instance to set the final things and make sure mongodb is set up correctly:
	- import ztfquery and query a location:
		>> from ztfquery import query
		>> q = query.ZTFQuery()
		>> q.load_metadata(radec=[7.1178554, 48.8587525], size=0.01, sql_query='obsjd BETWEEN 2458000 AND 2462000')
		>> q.metatable
		- As this is the first time, you will be asked for your irsa credentials when running the 2nd line
		- Give the credentials, they will be saved for the next time
		- If everything works correctly, the last line should return a DataFrame
	- import ForcedPhotometryPipeline from fpbot.pipeline and get a light curve:
		>> from fpbot.pipeline import ForcedPhotometryPipeline
		>> pl = ForcedPhotometryPipeline(file_or_name="ZTF20abqetet", ra=7.1178554, dec=48.8587525, nprocess=24, ampel=False)
		>> pl.download()
		>> pl.psffit()
		>> import pandas as pd
		>> lc = pd.read_csv('<Whatever $ZTFDATA was>forcephotometry/ZTF20abqetet.csv', header=0, comment='#')
		>> lc
		- If mongodb isn't set up properly the 2nd line will fail
		- If everything works correctly, the last line should return a DataFrame of >=4029 rows and 47 columns
- Open auto_check.py and add the following details:
	- Recepients (line 38)
	- Bot email address (line 94)
	- Bot email password (line 95)
	Note: if Bot email is not gmail, mailserver & mailport (lines 90 & 91) should be updated as well


=========================
|| Mongo configuration ||
=========================
The following lines should be copied into mongo_config.conf in the mongo folder and serve to correctly configure the mongo database to run as a daemon:
processManagement:
   fork: true
net:
   bindIp: localhost
   port: 27053
storage:
   dbPath: <path to the mongo folder>/db
systemLog:
   destination: file
   path: "<path to the mongo folder>/mongod.log"
   logAppend: true


================
|| User guide ||
================
- Before running the program for the 1st time, a few final settings need to be set in order for things to work properly:
	- Create the directories that do not exist yet (GitHub doesn't like empty folders, making this step necessary. See the contents for what is needed)
	- Update config.csv to your requirements
	- Make sure fpbot and ztfquery work properly (see installation guide above)
	- Set the email address details from which to send the result email (auto_ceck.py, see installation guide) (or disable email completely)
		- NOTE: When using e.g. a Gmail account for this, 2-step verification needs to be enabled before a token can be created to allow 3rd party programs like python to use it
	- Set the recepient list (auto_check.py, line 38)
	- (Optionally) run the programs on testlist.csv to make sure everything works
- Make the lists of objects that need to be checked.
- (Optional) run each list manually to generate the bulk of the light curves immediately. This way the daily check only has to update the light curves making it faster on the first run
- Get run_daily_check.sh to automatically run each day:
	- 


======================================
|| Known issues and quick solutions ||
======================================
- According to the output generated by lc_generator.py, IPAC has more images than I can download. I can't seem to get all of them for some reason?
	A: Unfortunately it is not possible to check for permissions when querying IPAC for what images it has. All are returned, but you will only be able to download those you are allowed to.
	   This results in the query usually finding more, but then failing to download the missing images. lc_generator.py should be able to handle these cases without issues.
	   Should you feel like you are still getting less data than you should, there might be an permissions issue with your account. Please contact IPAC to resolve it if this is the case.
- When downloading images I get a connection timeout after a certain amount of objects have been downloaded. How can I keep the program from timing out?
	A: IPAC has (at times) limits in place for the amount of requests that can be made per account per unit of time (e.g. 100 requests per minute) in order to accomodate all user requests.
	   fpbot runs with parallel downloads, which might cause it to hit this limit. Try lowering the amount of parallel processes run during the download (lc_generator.py, line 59)
- When downloading images it crashes, giving the following errror: "NameError: name 'session' is not defined" How do I fix this?
	A: This is a bug in ztfquery that has been reported on its GitHub page, but may not be patched yet. It occurs after the download, when checking for corrupted files.
	   If one is found, a redownload is initiated without properly defining the session. This bug can be fixed by though these steps:
		- Open <path to the used cona env>/lib/python3.10/site-packages/ztfquery/io.py (assuming you use python 3.10 in this environment
		- Find the line that causes the crash (line 720), it part of download_url(...) in the test_files() function.
		- Add the following directly before the call to download_url(...): session = open_irsa_session()
		- Make sure that it has the same indentation as download_url(...) (Watch out for the difference between a space and a tab)
