'''
A program to check the binning program output for interesting transients and send the results by email

Author: Jacco Terwel
Date: 24-10-23

Added:
	- Main structure and functions
'''

import sys
import pandas as pd
import smtplib
import shutil
from astropy.time import Time
from email.mime.application import MIMEApplication
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.utils import formatdate
from datetime import datetime, date
from pathlib import Path

def run_check(mode, tresh_dets, tresh_final_bins, recent):
	'''
	Main function of the program

	- mode (str): Mode this program is used in, currently only controls who the email is sent to
	- tresh_dets (int): Minimum nr. of binning attempts that found proper detections
	- tresh_final_bins (int): Minimum nr. of binning attempts where the last bin was a detection
	- recent (float): Up until how many days ago would still be considered currently active?
	'''
	if mode not in ['daily', 'manual']:
		print(f'ERROR: mode {mode} not recognised. Please use either daily or manual')
		print('Assuming manual mode and only sending an email to the first recepient only')
		mode = 'manual'
	# Check for interesting objects
	nr_old_dets, nr_active = check_fv(tresh_dets, tresh_final_bins, recent)
	# Send email or save in backup if the email fails
	try:
		for recepient in ['terwelj@tcd.ie']: # Might want to have recepients depending on what is being run (config file?)
			send_email(recepient, nr_old_dets, nr_active)
			# In manual mode, only send an email to the 1st person in the list
			if mode == 'manual':
				break
	except:
		make_backup()
	return

def check_fv(tresh_dets, tresh_final_bins, recent):
	'''
	Check the final verdicts for interesting objects and whether or not they are currently active

	- tresh_dets (int): Minimum nr. of binning attempts that found proper detections
	- tresh_final_bins (int): Minimum nr. of binning attempts where the last bin was a detection
	- recent (float): Up until how many days ago would still be considered currently active?
	'''
	# Read the final verdicts
	fv = pd.read_csv('results/final_verdicts.csv',
					 usecols=['name', 'suc_g', 'suc_r', 'suc_i', 'nr_suc_g', 'nr_suc_r', 'nr_suc_i',
					 		  'final_bin_det_g', 'last_obs_mjd_g', 'final_bin_det_r', 'last_obs_mjd_r',
					 		  'final_bin_det_i', 'last_obs_mjd_i'])
	# Select the objects that have robust late-time detections
	has_late_dets = fv[((fv.nr_suc_g>=tresh_dets)|(fv.nr_suc_r>=tresh_dets)|(fv.nr_suc_i>=tresh_dets))]
	# Select the objects that also have the last bin a detection
	final_bin_dets = has_late_dets[((has_late_dets.final_bin_det_g >= tresh_final_bins) |
									(has_late_dets.final_bin_det_r >= tresh_final_bins) |
									(has_late_dets.final_bin_det_i >= tresh_final_bins))]
	# Select the object where the last observation in the band whith the final bin a detection was recently
	mjd_tresh = Time(datetime.now()).mjd - recent
	currently_active = final_bin_dets[(((final_bin_dets.final_bin_det_g >= tresh_final_bins) &
										  	(final_bin_dets.last_obs_mjd_g >= mjd_tresh)) |
										((final_bin_dets.final_bin_det_r >= tresh_final_bins) &
										  	(final_bin_dets.last_obs_mjd_r >= mjd_tresh)) |
										((final_bin_dets.final_bin_det_i >= tresh_final_bins) &
										  	(final_bin_dets.last_obs_mjd_i >= mjd_tresh)))]
	# Make a DataFrame containing objects with dets that aren't in the currently_active DataFrame
	old_dets = has_late_dets[~has_late_dets.name.isin(currently_active.name.values)]
	# Save the currently active and old dets separately and return their sizes
	currently_active.reset_index(drop=True, inplace=True)
	old_dets.reset_index(drop=True, inplace=True)
	currently_active.to_csv('results/currently_active.csv')
	old_dets.to_csv('results/old_dets.csv')
	return len(old_dets), len(currently_active)

def send_email(recepient, nr_old_dets, nr_active): #Need to test if sending to multiple recepients at once is possible as well
	'''
	Send an email with the results

	- recepient (str): location to send the email to
	- nr_old_dets (int): length of the old_dets.csv DataFrame
	- nr_active (int): length of the currently_active.csv DataFrame
	'''
	# Define mailserver
	mailserver = "smtp.gmail.com"
	mailport = 587
	# Specify login details
	send_from = 'latetimesearcher@gmail.com'
	password = 'lotv hbcy youh dzoo'
	# Specify email
	subject = f'real-time late-time CSM interaction search results {date.today()}'
	text = f"The light curves have been generated, and the late-time observations have been binned.\nPlease check notes.txt to see the used list, number of objects queried, and any possible issues tht came up.\n{nr_active} objects that are currently active with late-time detections have been found, they are listed in currently_active.csv\n{nr_old_dets} objects with late-time detections that have already faded away again have been found, they are listed in old_dets.csv\nLet's see what we can find next time.\n\nCheers\n\nYour friendly neighbourhood robot"
	# Build the message
	msg = MIMEMultipart()
	msg["From"] = send_from
	msg["To"] = recepient
	msg["Date"] = formatdate(localtime=True)
	msg["Subject"] = subject
	msg.attach(MIMEText(text))
	# Attach files
	for loc in ['results/notes.txt', 'results/old_dets.csv', 'results/currently_active.csv']:
		file = Path(loc)
		# Check to make sure the file exists
		if file.is_file():
			with open(file, "rb") as csv:
			    part = MIMEApplication(csv.read(), Name="Dataframe")
			part["Content-Disposition"] = f'attachment; filename="{file.name}"'
			msg.attach(part)
	# Send email
	smtp = smtplib.SMTP(mailserver, mailport)
	smtp.starttls()
	smtp.ehlo()
	smtp.login(send_from, password)
	smtp.sendmail(send_from, recepient, msg.as_string())
	smtp.close()
	return

def make_backup():
	'''
	Move the files that were supposed to be emailed to the backup folder in a dated subfolder
	'''
	nu = datetime.date()
	van = Path('results')
	naar = Path(f'backup/{nu.day}-{nu.month}-{nu.year}')
	naar.mkdir(parents=True, exist_ok=True)
	for _ in ['notes.txt', 'old_dets.csv', 'currently_active.csv']:
		shutil.move(van/_, naar/_)
	return

if (__name__ == "__main__"):
	try:
		run_check(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
	except:
		run_check(sys.argv[1], 4, 4, 28)
