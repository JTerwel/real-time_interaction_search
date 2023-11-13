'''
A program to download all available ZTF observations of a set of objects and perform forced photometry on them

Author: Jacco Terwel
Date: 13-11-23

Added:
	- Main structure and functions
	- Own version of the fpbot download function to avoid having to download everything all the time
	- Added something to deal with the rate caps placed by irsa at times
	- Make list of updated lcs for the binning program to run on
'''

import sys
import os
import traceback
import pandas as pd
import shutil
import numpy as np
import time
from pathlib import Path
from fpbot.pipeline import ForcedPhotometryPipeline
from fpbot.utils import calculate_magnitudes, get_wise_ra_dec, is_wise_name, is_ztf_name
from fpbot.clean_lc import clean_lc
from ztfquery import query
from ztfquery.io import test_files
from astropy.time import Time
from datetime import datetime
from fpbot import database
from ztflc import forcephotometry

FORCEPHOTODATA = Path(os.getenv("ZTFDATA"))/"forcephotometry"
MARSHALDATA = Path(os.getenv("ZTFDATA"))/"marshal"

def gen_lcs(listloc, lc_dir, nprocess_download, nprocess_fit):
	'''
	Main function of the program

	- listloc (str): Location of the list containing the objects to query
	- lc_dir (str): Location of the directory containing the light curves that are already present
	'''
	#Make sure the lc directory exists before starting
	lc_dir = Path(lc_dir)
	if not lc_dir.is_dir():
		lc_dir.mkdir()
	# Load the list
	listloc = Path(listloc)
	objs = load_list(listloc)
	# Keep notes throughout the process & set where they should be saved
	noteloc = Path('results/notes.txt')
	notes = f'Time started: {datetime.now()}\n'
	lcs_generated = 0
	#Load blacklisted objects that need to be skipped
	blacklist = pd.read_csv('blacklist.csv', comment='#')
	# For each object, check if new images have become available since last time
	# If so, run the fpbot pipeline on that object and create a light curve
	for i in range(len(objs)): # Possible to do this with multiprocessing as well or the download clog up?
		if objs.ztfname.loc[i] in blacklist.ztfname.values: #skip blacklisted objects
			continue
		try:
			# Check if an old lc of this object already exists
			if (lc_dir/f'{objs.ztfname.loc[i]}.csv').is_file():
				lc_loc = lc_dir/f'{objs.ztfname.loc[i]}.csv'
			else:
				lc_loc = None
			# Generate the light curve
			nr_files, download_time = run_fpbot(objs.ztfname.loc[i], objs.ra.loc[i], objs.dec.loc[i],
												objs.last_checked_nr.loc[i], lc_loc, lc_dir,
												nprocess_download, nprocess_fit)
			# Update last checked mjd and nr
			objs.loc[i, 'last_checked_mjd'] = download_time
			objs.loc[i, 'last_checked_nr'] = nr_files
			lcs_generated += 1
		except: #Make a note of the error and continue with the next object
			notes += f'\n\nSomething went wrong with {objs.ztfname.loc[i]}, this is the error traceback:\n'.join(traceback.format_stack())
			notes += '\n\n'
			continue
	# Save the notes
	notes += f'Out of {len(objs)} that were in {listloc.name}, {lcs_generated} had new data and had their new light curves generated\nTime finished: {datetime.now()}'
	ftext = open(noteloc, 'w')
	ftext.write(notes)
	ftext.close()
	# Save the updated version of objs if it has been altered
	if lcs_generated > 0:
		objs.to_csv(listloc)
	return

def load_list(path):
	'''
	Load the list of objects to check

	- path (Path): location of the list
	'''
	objs = pd.read_csv(path, header=0, index_col=0)
	# Check if all needed columns are present: Name, loc, when and nr. of images when last lc was made.
	for _ in ['ztfname', 'ra', 'dec', 'last_checked_mjd', 'last_checked_nr']:
		if _ not in objs.columns:
			raise ValueError(f'{_} not found in DataFrame columns, not all needed information is present')
	return objs

def check_available_ims(ra, dec):
	'''
	Use ztfquery to check what images are available for download
	Note that the query searches between 4-9-2017 and 17-08-2028. This should be enough for now

	- ra  (float): Object RA
	- dec (float): Object Dec
	'''
	q = query.ZTFQuery()
	q.load_metadata(radec=[ra, dec], size=0.01, sql_query='obsjd BETWEEN 2458000 AND 2462000')
	return q.metatable

def run_fpbot(name, ra, dec, last_checked_nr, lc_loc, lc_dir, nprocess_download, nprocess_fit=64):
	'''
	Download the cutouts and make a light curve at the given position

	- name  (str): Object name
	- ra  (float): Object RA
	- dec (float): Object Dec
	'''
	#Need to play with nprocesses		###
	pl = ForcedPhotometryPipeline(file_or_name=name, ra=ra, dec=dec, nprocess=nprocess_fit, ampel=False)
	#Add custom download function to use
	pl.new_download_function = new_download_function
	nr_files, download_time = pl.new_download_function(pl, last_checked_nr, lc_loc,
													   nprocess_download=nprocess_download)
	pl.new_psffit_function = new_psffit_function
	pl.new_psffit_function(pl, lc_loc)
	if lc_loc is not None:
		# Add old and new lc together (header is lost, but I don't really care for it)
		lc_old = pd.read_csv(lc_loc, header=0, comment='#')
		try:
			lc_new = pd.read_csv(FORCEPHOTODATA/(name+'.csv'), header=0, comment='#')
			lc = pd.concat([lc_old, lc_new], ignore_index=True)
			#Check and remove duplicated rows (consider rows with the same filename as duplicates)
			lc = lc.loc[~lc.duplicated(subset='filename', keep='first'), :]
			print(f'Old lc has {len(lc_old)} lines, new lc has {len(lc_new)} lines, combined lc has {len(lc)} lines, the rest were duplicates')
			lc.to_csv(lc_loc, index=False)
		except:
			print(f'No new light curve was made for {name}')
	else:
		# Move lc to storage location
		shutil.move(FORCEPHOTODATA/(name+'.csv'), lc_dir/(name+'.csv'))
	return nr_files, download_time

def new_download_function(self, last_checked_nr, lc_loc, no_local_check=False, nprocess_filecounts=16,
						  nprocess_download=32):
	'''
	Mimic the download function in the ForcedPhotometryPipeline but with a few changes
	(may have to import some things again as this function works from somewhere else, e.g. database)
	Changes made:
	- Added a proper return at the end of the function, returns ipac file nr & last download time
	- Need to pass lc_loc, the location of the lc we already have
	- Need to pass last_checked_nr, which is the last image count, can ingore with no_local_check=True
	- Removed outdir and dummy file lines as they aren't used anywhere
	- Added the option to change nprocess for the filecount & download
	- Might build option to disable logger later
	- Reformatted a bit to my style
	'''
	number_of_objects = len(self.object_list)
	download_requested = []
	query = database.read_database(self.object_list, ["ra", "dec", "last_download"])
	last_download = query["last_download"]
	# In case no_new_downloads option is passed (download_newest = False): Download only if it has never been downloaded before (useful for bulk downloads which repeatedly fail because IPAC is unstable) Else: try to download everything.
	if self.download_newest is False:
		for index, name in enumerate(self.object_list):
			if last_download[index] is None:
				download_requested.append(name)
	else:
		download_requested = self.object_list
	# Check with IRSA how many images are present for each object. Only if this number is bigger than the local number of images, download will start.
	download_needed = []
	query = database.read_database(download_requested, ["ra", "dec"])
	ras = query["ra"]
	decs = query["dec"]
	from fpbot.connectors import get_ipac_and_local_filecount
	self.logger.info(f'Trying to get the filecount using nprocess={nprocess_filecounts}')
	ipac_filecounts = get_ipac_and_local_filecount(ztf_names=download_requested, ras=ras, decs=decs,
												   jdmin=self.jdmin, jdmax=self.jdmax,
												   nprocess=nprocess_filecounts)
	for index, name in enumerate(download_requested): #Check if download is required
		if ipac_filecounts[name]["local"] < ipac_filecounts[name]["ipac"]: #Compare with local ims
			if last_checked_nr < ipac_filecounts[name]["ipac"] or no_local_check: #Compare with old lc size
				download_needed.append(name)
	if download_needed:
		self.logger.info(f"{len(download_needed)} of {len(self.object_list)} objects have additional images available at IRSA (IRSA: {ipac_filecounts[name]['ipac']}, local: {ipac_filecounts[name]['local']}).\nThese will be downloaded now.")
	else:
		self.logger.info("For none of the transients are new images available, so no download needed.")
	for i, name in enumerate(download_needed):
		query = database.read_database(name, ["ra", "dec", "jdmin", "jdmax", "local_filecount"])
		ra = query["ra"][0]
		dec = query["dec"][0]
		jdmin = self.jdmin
		jdmax = self.jdmax
		fp = forcephotometry.ForcePhotometry.from_coords(ra=ra, dec=dec, jdmin=jdmin, jdmax=jdmax, name=name)
		self.logger.info(f"{name} ({i+1} of {len(download_needed)}) Downloading data.")
		fp.load_metadata(clean=self.ztfquery_clean_metatable)
		#metatable has been loaded, time to remove the rows whose points we already have
		if lc_loc is not None:
			#load the lc
			lc = pd.read_csv(lc_loc, header=0, comment='#', usecols=['filename', 'flag'])
			#Get the names associated with each metatable line
			metatable_file_names = [Path(i).name.split('_q')[0]+'.fits' for i in fp.io.zquery.get_data_path(suffix='diffimgpsf.fits', source="None")]
			#Keep only those lines that aren't in lc.filename
			keepers = [i not in lc.filename.values for i in metatable_file_names]#Not sure if this evaluation works
			fp.io.zquery._metatable = fp.io.zquery.metatable[keepers]
		#Should now only have those lines whose data we don't have yet
		which = ["scimrefdiffimg.fits.fz", "diffimgpsf.fits"]
		#First get the download location
		download_url, download_location = fp.io.zquery.download_data(download_dir=self.get_outdir(name),
																	 cutouts=True, radec=[ra, dec],
																	 cutout_size=30,
																	 nprocess=nprocess_download, 
																	 overwrite=False, show_progress=True,
																	 ignore_warnings=True, nodl=True)
		self.logger.info(f'Trying to download using nprocess={nprocess_download}')
		fp.io.download_data(download_dir=self.get_outdir(name), cutouts=True, radec=[ra, dec],
							cutout_size=30, nprocess=nprocess_download, overwrite=False, show_progress=True,
							ignore_warnings=True, which=which)
		if self.sciimg:
			which = ["scimrefdiffimg.fits.fz", "diffimgpsf.fits", "sciimg.fits"]
			self.logger.info(f'Trying to download using nprocess={nprocess_download}')
			fp.io.download_data(download_dir=self.get_outdir(name), cutouts=True, radec=[ra, dec],
								cutout_size=30, nprocess=nprocess_download, overwrite=False,
								show_progress=True, ignore_warnings=True, which=which)
		last_download = Time(time.time(), format="unix", scale="utc").jd
		database.update_database(name, {"lastdownload": last_download})
	return ipac_filecounts[name]["ipac"], last_download

def new_psffit_function(self, lc_loc, nprocess=None, force_refit=False):
	"""
	Mimic the psffit function in the ForcedPhotometryPipeline but with a few changes
	(may have to import some things again as this function works from somewhere else, e.g. database)
	Changes made:
	- Added a proper return at the end of the function, returns ipac file nr & last download time
	- Added the lines selecting only the new datapoints to work on (so it doesn't use non-existing ones)
	- Removed the lines to make the header in the csv file, I'm not using it anyway.
	- Removed dummyfile lines
	- Reformatted a bit to my style
	"""
	if nprocess is None:
		nprocess = self.nprocess
	query = database.read_database(self.object_list, ["ra", "dec", "jdmin", "jdmax", "lastobs",
													  "lastdownload", "lastfit", "coords_per_filter",
													  "fitted_datapoints",])
	for i, name in enumerate(self.object_list):
		objects_total = len(self.object_list)
		ra = query["ra"][i]
		dec = query["dec"][i]
		jdmin = self.jdmin
		jdmax = self.jdmax
		lastobs = query["lastobs"][i]
		lastdownload = query["lastdownload"][i]
		lastfit = query["lastfit"][i]
		coords_per_filter = query["coords_per_filter"][i]
		fitted_datapoints = query["fitted_datapoints"][i]
		"""
		Automatically rerun fit if last fit was before
		March 24, 2022 (to ensure header and quality flags)
		"""
		if lastfit:
			if lastfit < 2459662.50000:
				force_refit = True
		"""
		Check if there are different centroids for the
		different filters
		If a filter is missing, replace with total (all filters) median RA/Dec.
		If the source is a WISE object, we
		only have one RA/Dec
		"""
		if is_wise_name(name):
			coords_per_filter = [ra, dec]
		coords_per_filter[0] = np.nan_to_num(x=coords_per_filter[0], nan=ra).tolist()
		coords_per_filter[1] = np.nan_to_num(x=coords_per_filter[1], nan=dec).tolist()
		fp = forcephotometry.ForcePhotometry.from_coords(ra=coords_per_filter[0], dec=coords_per_filter[1],
														 jdmin=jdmin, jdmax=jdmax, name=name)
		self.logger.info(f"{name} ({i+1} of {objects_total}) loading metadata.")
		fp.load_metadata(clean=self.ztfquery_clean_metatable)
		#metatable has been loaded, time to remove the rows whose points we already have
		if lc_loc is not None:
			#load the lc
			lc = pd.read_csv(lc_loc, header=0, comment='#', usecols=['filename', 'flag'])
			#Get the names associated with each metatable line
			metatable_file_names = [Path(i).name.split('_q')[0]+'.fits' for i in fp.io.zquery.get_data_path(suffix='diffimgpsf.fits', source="None")]
			#Keep only those lines that aren't in lc.filename
			keepers = [i not in lc.filename.values for i in metatable_file_names]#Not sure if this evaluation works
			fp.io.zquery._metatable = fp.io.zquery.metatable[keepers]
		#Should now only have those lines whose data we don't have yet
		self.logger.info(f"{name} ({i+1} of {objects_total}) metadata loaded.")
		self.logger.info(f"{name} ({i+1} of {objects_total}) loading paths to files.")
		fp.load_filepathes(filecheck=self.filecheck, download_dir=self.get_outdir(name))
		self.logger.info(f"{name} ({i+1} of {objects_total}) paths to files loaded.")
		"""
		Check how many forced photometry datapoints
		there SHOULD exist for this object
		"""
		number_of_fitted_datapoints_expected = len(fp.filepathes)
		if fitted_datapoints is None:
			fitted_datapoints = 0
		df_file = os.path.join(FORCEPHOTODATA, f"{name}.csv")
		if os.path.isfile(df_file):
			if os.stat(df_file).st_size > 1:
				_df = pd.read_csv(df_file, comment="#", index_col=0)
				if len(_df) == 0:
					force_refit = True
			else:
				os.remove(df_file)
				force_refit = True
		else:
			force_refit = True
		# Compare to number of fitted datapoints from database
		if number_of_fitted_datapoints_expected > fitted_datapoints or force_refit:
			self.logger.info(f"{name} ({i+1} of {objects_total}): Fitting PSF.")
			fp.run_forcefit(nprocess=nprocess, store=True, force_refit=force_refit, no_badsub=False)
			self.logger.info(f"{len(fp.data_forcefit)} points fitted")
			fp.store()
			lastfit = Time(time.time(), format="unix", scale="utc").jd
			df_file = os.path.join(FORCEPHOTODATA, f"{name}.csv")
			if len(fp.data_forcefit) > 0: #df_file only exists is it had been saved before
				df = pd.read_csv(df_file, comment="#", index_col=0)
				# Calculate cloudiness parameter, add 'pass' column (only rows with pass=1 should be used)
				if len(df) > 0:
					df = clean_lc(df, trim=False)
				os.remove(df_file)
				f = open(df_file, "a")
				df.to_csv(f, index=False)
				f.close()
				database.update_database(name, {"lastfit": lastfit,
												"fitted_datapoints": number_of_fitted_datapoints_expected})
		else:
			self.logger.info(f"{name} ({i+1} of {objects_total}) No new images to fit, skipping PSF fit.")
	return

if (__name__ == "__main__"):
	try:
		gen_lcs(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
	except:
		gen_lcs(sys.argv[1], sys.argv[2], 32, 32)