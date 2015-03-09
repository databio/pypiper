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

import atexit
#import psutil

# Define global variables
# These just record some details about the pipeline

PIPELINE_NAME = ""
PEAKMEM = 0			# memory high water mark
STARTTIME = time()
LAST_TIMESTAMP = STARTTIME	# time of the last call to timestamp()
STATUS = "initializing"
PIPELINE_OUTFOLDER = ""

def exit_handler():
	# Catch-all for uncaught exceptions...
	global STATUS
	if STATUS != "completed":
		set_status_flag("failed")


def set_status_flag(status):
	global PIPELINE_NAME
	global PIPELINE_OUTFOLDER
	global STATUS
	prev_status = STATUS
	# Remove previous status
	flag_file = PIPELINE_OUTFOLDER + "/" + PIPELINE_NAME + "_" + prev_status + ".flag"
	try:
		os.remove(os.path.join(flag_file))
	except:
		pass

	# Set new status
	STATUS = status
	flag_file = PIPELINE_OUTFOLDER + "/" + PIPELINE_NAME + "_" + status + ".flag"
	create_file(flag_file)

	print ("Change status from " + prev_status + " to " + status)


def start_pipeline(paths, args, pipeline_name):
	"""Do some setup, like tee output, print some diagnostics, create temp files"""
	# add variables for this pipeline
	atexit.register(exit_handler)

	global PIPELINE_NAME
	PIPELINE_NAME = pipeline_name
	paths.pipeline_outfolder = os.path.join(args.project_root + args.sample_name + "/")
	paths.pipe_stats = paths.pipeline_outfolder + "/" + pipeline_name + "_stats.txt"
	paths.log_file = paths.pipeline_outfolder + pipeline_name  + "_log.md"
	make_sure_path_exists(paths.pipeline_outfolder)
	global PIPELINE_OUTFOLDER
	PIPELINE_OUTFOLDER = paths.pipeline_outfolder

	global STARTTIME
	STARTTIME = time()
	# Mirror every operation on sys.stdout to log file
	sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0) # Unbuffer output
	# a for append to file
	tee = subprocess.Popen(["tee", "-a", paths.log_file], stdin=subprocess.PIPE)
	os.dup2(tee.stdin.fileno(), sys.stdout.fileno())
	os.dup2(tee.stdin.fileno(), sys.stderr.fileno())

	# Record the git version of the pipeline being run. This code gets:
	# hash: the commit id of the last commit in this repo
	# date: the date of the last commit in this repo
	# diff: a summary of any differences in the current (run) version vs. the committed version
	git_commit_hash = subprocess.check_output("cd " + os.path.dirname(os.path.realpath(__file__)) + "; git rev-parse --verify HEAD", shell=True)
	git_commit_date = subprocess.check_output("cd " + os.path.dirname(os.path.realpath(__file__)) + "; git show -s --format=%ai HEAD", shell=True)
	git_commit_diff = subprocess.check_output("cd " + os.path.dirname(os.path.realpath(__file__)) + "; git diff --shortstat HEAD", shell=True)
	if (git_commit_diff==""):
		git_commit_diff = "No uncommitted changes."
	start_time = time()
	print("################################################################################")
	timestamp("Pipeline started at: ")
	print "Cmd: " + str(" ".join(sys.argv))
	print "Working dir : %s" % os.getcwd()
	print "Run outfolder:\t\t" + paths.pipeline_outfolder
	print "Compute host:\t\t" + platform.node()
	print "Git pipeline version:\t\t" + git_commit_hash.strip()
	print "Git pipeline date:\t\t" + git_commit_date.strip()
	print "Git pipeline diff: \t\t" + git_commit_diff.strip()
	print("Python version:\t\t" + platform.python_version())
	print("Project root:\t\t" + args.project_root)
	print("Paired end mode:\t\t" + str(args.paired_end))
	print("################################################################################")

	set_status_flag("running")

	return paths

def fail_pipeline(reason):
	set_status_flag("failed")
	raise Exception(reason)

def make_sure_path_exists(path):
	try:
		os.makedirs(path)
	except OSError as exception:
		if exception.errno != errno.EEXIST:
			raise


def wait_for_lock(lock_file):
	sleeptime = 5
	first_message_flag = False
	dot_count=0
	while os.path.isfile(lock_file):
		if first_message_flag == False:
			timestamp("Waiting for file lock: " + lock_file)
		else:
			sys.stdout.write(".")
			dot_count = dot_count+1
			if dot_count == 60:
				print "" # linefeed
				dot_count = 0
		sleep(sleeptime)
		sleeptime = min(sleeptime+5, 120)
	timestamp("File unlocked.")


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
		if isinstance(cmd, list):
			for cmd_i in cmd:
				ret, local_maxmeme = callprint(cmd_i, shell)
			else:
				ret, local_maxmem = callprint(cmd, shell)	# Run command
		os.remove(lock_file)		# Remove lock file
	else:
		print("File already exists: " + output_file)
	return ret



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
		fail_pipeline("Process returned nonzero result.")
		#raise Exception("Process returned nonzero result.")
	return [p.returncode, local_maxmem]



def time_elapsed(time_since):
	return round(time() - time_since,2)



def stop_pipeline(paths):
	global PEAKMEM
	global STARTTIME
	"""Remove temporary marker files to complete the pipeline"""
	set_status_flag("completed")
	print ("Total elapsed time: " + str(time_elapsed(STARTTIME)))
	#print ("Peak memory used: " + str(memory_usage()["peak"]) + "kb")
	print ("Peak memory used: " + str(PEAKMEM/1e6) + " GB")
	timestamp("### Pipeline completed at: ");


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
