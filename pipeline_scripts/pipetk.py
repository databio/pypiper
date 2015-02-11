#!/usr/bin/python2.7
"""
Pipe ToolKit (pipetk) documentation
"""

import os
import errno
import re
import platform
import sys
from subprocess import call
import subprocess
from time import sleep, time, strftime


def make_sure_path_exists(path):
	try:
		os.makedirs(path)
	except OSError as exception:
		if exception.errno != errno.EEXIST:
			raise


def wait_for_lock(lock_file):
	sleeptime = 10
	while os.path.isfile(lock_file):
		print "Waiting for file lock (" + str(sleeptime) + " sec): " + lock_file
		sleep(sleeptime)
		sleeptime = min(sleeptime+5, 3600)


def create_file(file):
	with open(file, 'w') as fout:
		fout.write('')
	fout.close()


def timestamp(message):
	message += " (" + strftime("%Y-%m-%d %H:%M:%S") + ")"
	if re.match("^###", message):
		message = "\n" + message + "\n"
	print(message)


def callprint(cmd):
	print(cmd)
	call(cmd, shell=True)



def call_lock(cmd, lock_name, folder, output_file=None):
	# Create lock file:
	lock_file = os.path.join(folder,  lock_name)
	wait_for_lock(lock_file)
	if output_file is None or not os.path.isfile(output_file):
		create_file(lock_file)		# Create lock
		callprint(cmd)				# Run command
		os.remove(lock_file)		# Remove lock file
	else:
		print("File already exists: " + output_file)


def time_elapsed(time_since):
	return round(time() - time_since,2)

def start_pipeline(paths, args):
	"""Do some setup, like tee output, print some diagnostics, create temp files"""
	make_sure_path_exists(paths.pipeline_outfolder)

	# Mirror every operation on sys.stdout to log file
	sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0) # Unbuffer output
	# a for append to file
	tee = subprocess.Popen(["tee", "-a", paths.log_file], stdin=subprocess.PIPE)
	os.dup2(tee.stdin.fileno(), sys.stdout.fileno())
	os.dup2(tee.stdin.fileno(), sys.stderr.fileno())
	start_time = time()
	print("################################################################################")
	timestamp("Script start time: ")
	print "Run outfolder:\t\t" + paths.pipeline_outfolder
	print "Compute host:\t\t" + platform.node()
	print("Python version:\t\t" + platform.python_version())
	print("Project root:\t\t" + args.project_root)
	print("Paired end mode:\t\t" + str(args.paired_end))
	print("################################################################################")

	# Create a temporary file to indicate that this pipeline is currently running in this folder.
	pipeline_temp_marker = paths.pipeline_outfolder + "/" + "WGBS-running.temp"
	create_file(pipeline_temp_marker)
	return start_time

def stop_pipeline(paths, args, start_time):

	"""Remove temporary marker files to complete the pipeline"""
	pipeline_temp_marker = paths.pipeline_outfolder + "/" + "WGBS-running.temp"
	os.remove(os.path.join(pipeline_temp_marker))

	timestamp("### Script end time: ");
	print ("Total elapsed time: " + str(time_elapsed(start_time)))

