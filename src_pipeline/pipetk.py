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
	sleeptime = 15
	while os.path.isfile(lock_file):
		print "Waiting for file lock (" + str(sleeptime) + " sec): " + lock_file
		sleep(sleeptime)
		sleeptime = min(sleeptime+5, 3600)


def create_file(file):
	with open(file, 'w') as fout:
		fout.write('')
	fout.close()


def timestamp(message, time_since=None):
	message += " (" + strftime("%Y-%m-%d %H:%M:%S") + ")"
	if time_since!=None:
		message += " elapsed:" + time_elapsed(time_since)
	if re.match("^###", message):
		message = "\n" + message + "\n"
	print(message)


# split the command to use shell=False;
# leave it together to use shell=True; I use False so I can get the PID
# and poll memory use.
def callprint(cmd, shell=False):
	print(cmd)
	if not shell:
		cmd = cmd.split()
	call(cmd, shell=shell)



def report_result(key, value, paths):
	with open(paths.pipe_stats, "a") as myfile:
		myfile.write(key + ": " + str(value).strip() + "\n")




def call_lock(cmd, lock_name, folder, output_file=None):
	# Create lock file:
	lock_file = os.path.join(folder,  lock_name)
	wait_for_lock(lock_file)
	if output_file is None or not os.path.isfile(output_file):
		create_file(lock_file)		# Create lock
		callprint(cmd)				# Run command
		os.remove(lock_file)		# Remove lock file
		print ("Peak memory: " + "{:,}".format(round(memory_usage()["peak"]/1e6, 0)) + "gb")
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
	print "Cmd: " + str(" ".join(sys.argv))
	print "Working dir : %s" % os.getcwd()
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
	#print ("Peak memory used: " + str(memory_usage()["peak"]) + "kb")
	print ("Peak memory used: " + "{:,}".format(memory_usage()["peak"]) + "kb")





# Thanks Martin Geisler:
def memory_usage():
    """Memory usage of the current process in kilobytes."""
    status = None
    result = {'peak': 0, 'rss': 0}
    try:
        # This will only work on systems with a /proc file system
        # (like Linux).
        status = open('/proc/self/status')
        for line in status:
            parts = line.split()
            key = parts[0][2:-1].lower()
            if key in result:
                result[key] = int(parts[1])
    finally:
        if status is not None:
            status.close()
    return result
