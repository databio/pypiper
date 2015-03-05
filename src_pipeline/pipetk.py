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
#import psutil

# Define global variables
# These just record some details about the pipeline

PIPELINE_NAME = ""
PEAKMEM = 0			# memory high water mark
STARTTIME = time()
LAST_TIMESTAMP = STARTTIME	# time of the last call to timestamp()


def make_sure_path_exists(path):
	try:
		os.makedirs(path)
	except OSError as exception:
		if exception.errno != errno.EEXIST:
			raise


def wait_for_lock(lock_file):
	sleeptime = 5
	while os.path.isfile(lock_file):
		print "Waiting for file lock (" + str(sleeptime) + " sec): " + lock_file
		sleep(sleeptime)
		sleeptime = min(sleeptime+5, 120)


def create_file(file):
	with open(file, 'w') as fout:
		fout.write('')
	fout.close()


def timestamp(message):
	global LAST_TIMESTAMP
	message += " (" + strftime("%m-%d %H:%M:%S") + ")"
	message += " elapsed:" + str(time_elapsed(LAST_TIMESTAMP))
	message += " _TIME_"
	if re.match("^###", message):
		message = "\n" + message + "\n"
	print(message)
	LAST_TIMESTAMP = time()


# split the command to use shell=False;
# leave it together to use shell=True; I use False so I can get the PID
# and poll memory use.
def callprint(cmd, shell=False):
	global PEAKMEM
	print(cmd)
	if not shell:
		if ("|" in cmd or ">" in cmd):
			print("Should this command run in a shell instead of directly in a subprocess?")
		cmd = cmd.split()
	#call(cmd, shell=shell) # old way (no memory profiling)

	p = subprocess.Popen(cmd, shell=shell)
	local_maxmem = -1
	sleeptime=5
	while p.poll() == None:
		if not shell:
			local_maxmem = max(local_maxmem, memory_usage(p.pid))
			#print ("int.maxmem (pid:" + str(p.pid) + ") " + str(local_maxmem))
		sleep(sleeptime)
		sleeptime = min(sleeptime+5, 60)

	# set global maxmem
	PEAKMEM = max(PEAKMEM, local_maxmem)

	info = "Process " + str(p.pid) + " returned: (" + str(p.returncode) + ")."
	info += " Peak memory: (Process: " + str(local_maxmem) + "b;"
	info += " Pipeline: " +  str(PEAKMEM) +"b)"
	print (info)
	if p.returncode != 0:
		raise Exception("Process returned nonzero result.")
	return [p.returncode, local_maxmem]



def report_result(key, value, paths):
	message = key + "\t " + str(value).strip()
	print(message + "\t" + "_RES_")
	with open(paths.pipe_stats, "a") as myfile:
		myfile.write(message + "\n")




def call_lock(cmd, lock_name, folder, output_file=None, shell=False):
	# Create lock file:
	lock_file = os.path.join(folder,  lock_name)
	wait_for_lock(lock_file)
	ret = 0
	if output_file is not None:
		print ("Looking for file: " + output_file)
	if output_file is None or not (os.path.exists(output_file)):
		create_file(lock_file)		# Create lock
		ret, local_maxmem = callprint(cmd, shell)				# Run command
		os.remove(lock_file)		# Remove lock file
	else:
		print("File already exists: " + output_file)

	return ret


def time_elapsed(time_since):
	return round(time() - time_since,2)


def start_pipeline(paths, args, pipeline_name):
	"""Do some setup, like tee output, print some diagnostics, create temp files"""
	# add variables for this pipeline
	global PIPELINE_NAME
	PIPELINE_NAME = pipeline_name
	paths.pipeline_outfolder = os.path.join(args.project_root + args.sample_name + "/")
	paths.pipe_stats = paths.pipeline_outfolder + "/" + "stats_" + pipeline_name
	paths.log_file = paths.pipeline_outfolder + pipeline_name  + ".log.md"
	make_sure_path_exists(paths.pipeline_outfolder)
	global STARTTIME
	STARTTIME = time()
	# Mirror every operation on sys.stdout to log file
	sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0) # Unbuffer output
	# a for append to file
	tee = subprocess.Popen(["tee", "-a", paths.log_file], stdin=subprocess.PIPE)
	os.dup2(tee.stdin.fileno(), sys.stdout.fileno())
	os.dup2(tee.stdin.fileno(), sys.stderr.fileno())
	start_time = time()
	print("################################################################################")
	timestamp("Script start time: ")
	print "Cmd: " + str(" ".join(sys.argv))
	print "Working dir : %s" % os.getcwd()
	print "Run outfolder:\t\t" + paths.pipeline_outfolder
	print "Compute host:\t\t" + platform.node()
	print("Python version:\t\t" + platform.python_version())
	print("Project root:\t\t" + args.project_root)
	print("Paired end mode:\t\t" + str(args.paired_end))
	print("################################################################################")

	# Create a temporary file to indicate that this pipeline is currently running in this folder.
	pipeline_temp_marker = paths.pipeline_outfolder + "/" + pipeline_name + "-running.temp"
	create_file(pipeline_temp_marker)


	return paths


def stop_pipeline(paths):
	global PEAKMEM
	global STARTTIME
	global PIPELINE_NAME
	"""Remove temporary marker files to complete the pipeline"""
	pipeline_temp_marker = paths.pipeline_outfolder + "/" + PIPELINE_NAME + "-running.temp"
	os.remove(os.path.join(pipeline_temp_marker))
	pipeline_done_marker = paths.pipeline_outfolder + "/" + PIPELINE_NAME + "-completed"
	create_file(pipeline_done_marker)
	timestamp("### Script end time: ");
	print ("Total elapsed time: " + str(time_elapsed(STARTTIME)))
	#print ("Peak memory used: " + str(memory_usage()["peak"]) + "kb")
	print ("Peak memory used: " + str(PEAKMEM/1e6) + " GB")


# Thanks Martin Geisler:
def memory_usage(pid='self', category="peak"):
	"""Memory usage of the current process in kilobytes."""
	status = None
	result = {'peak': 0, 'rss': 0}
	try:
		# This will only work on systems with a /proc file system
		# (like Linux).
		#status = open('/proc/self/status')
		proc_spot = '/proc/%s/status' % pid
		status = open(proc_spot)
		for line in status:
			parts = line.split()
			key = parts[0][2:-1].lower()
			if key in result:
				result[key] = int(parts[1])
	finally:
		if status is not None:
			status.close()
	return result[category]
