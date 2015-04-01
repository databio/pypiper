#!/usr/env python

"""
Pypiper is a python module with two components: 1) the Pypiper class, and 2) other toolkits with functions for more
specific pipeline use-cases. The Pypiper class can be used to create a procedural pipeline in python.
"""

import os
import errno
import re
import platform
import sys
# from subprocess import call
import subprocess
from time import sleep, time, strftime

import atexit
import signal



class Pypiper:
	'Base class for instantiating a pypiper pipeline object'
	# Define pipeline-level variables to keep track of global state and some pipeline stats

	PIPELINE_NAME = ""
	PEAKMEM = 0         # memory high water mark
	STARTTIME = time()
	LAST_TIMESTAMP = STARTTIME  # time of the last call to timestamp()
	STATUS = "initializing"
	RUNNING_SUBPROCESS = None

	# Some file paths:
	PIPELINE_OUTFOLDER = ""
	PIPELINE_STATS = ""
	LOGFILE = ""  # required for termination signal code only.


	def __init__(self, name, ):
		self.PIPELINE_NAME = name


	def start_pipeline(self, paths, args):
		"""Do some setup, like tee output, print some diagnostics, create temp files"""
		# add variables for this pipeline
		atexit.register(self.exit_handler)

		paths.pipeline_outfolder = os.path.join(args.project_root + args.sample_name + "/")
		paths.pipe_stats = paths.pipeline_outfolder + "/" + self.PIPELINE_NAME + "_stats.txt"
		paths.log_file = paths.pipeline_outfolder + self.PIPELINE_NAME + "_log.md"

		self.LOGFILE = paths.log_file
		self.PIPELINE_OUTFOLDER = paths.pipeline_outfolder
		self.PIPELINE_STATS = paths.pipe_stats

		self.make_sure_path_exists(paths.pipeline_outfolder)

		self.STARTTIME = time()

		# Register handler functions to deal with interrupt and termination signals;
		# If received, we would then clean up properly (set pipeline status to FAIL, etc).
		signal.signal(signal.SIGINT, self.signal_int_handler)
		signal.signal(signal.SIGTERM, self.signal_term_handler)

		# Mirror every operation on sys.stdout to log file
		sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)  # Unbuffer output
		# a for append to file
		tee = subprocess.Popen(["tee", "-a", self.LOGFILE], stdin=subprocess.PIPE)
		# In case this pipeline process is terminated with SIGTERM, make sure we kill this spawned process as well.
		atexit.register(self.kill_child, tee.pid)
		os.dup2(tee.stdin.fileno(), sys.stdout.fileno())
		os.dup2(tee.stdin.fileno(), sys.stderr.fileno())

		# Record the git version of the pipeline being run. This code gets:
		# hash: the commit id of the last commit in this repo
		# date: the date of the last commit in this repo
		# diff: a summary of any differences in the current (run) version vs. the committed version
		git_pypiper_commit_hash = subprocess.check_output("cd " + os.path.dirname(os.path.realpath(__file__)) + "; git rev-parse --verify HEAD", shell=True)
		git_pypiper_commit_date = subprocess.check_output("cd " + os.path.dirname(os.path.realpath(__file__)) + "; git show -s --format=%ai HEAD", shell=True)
		git_pypiper_commit_diff = subprocess.check_output("cd " + os.path.dirname(os.path.realpath(__file__)) + "; git diff --shortstat HEAD", shell=True)
		git_pipe_commit_hash = subprocess.check_output("cd " + os.path.dirname(os.path.realpath(sys.argv[0])) + "; git rev-parse --verify HEAD", shell=True)
		git_pipe_commit_date = subprocess.check_output("cd " + os.path.dirname(os.path.realpath(sys.argv[0])) + "; git show -s --format=%ai HEAD", shell=True)
		git_pipe_commit_diff = subprocess.check_output("cd " + os.path.dirname(os.path.realpath(sys.argv[0])) + "; git diff --shortstat HEAD", shell=True)

		# Print out a header section in the pipeline log:
		print("################################################################################")
		self.timestamp("Pipeline started at: ")
		print("Compute host:\t\t" + platform.node())
		print("Working dir : %s" % os.getcwd())
		print("Git pypiper version:\t\t" + git_pypiper_commit_hash.strip())
		print("Git pypiper date:\t\t" + git_pypiper_commit_date.strip())
		if (git_pypiper_commit_diff != ""):
			print("Git pypiper diff: \t\t" + git_pypiper_commit_diff.strip())
		print("Git pipeline version:\t\t" + git_pipe_commit_hash.strip())
		print("Git pipeline date:\t\t" + git_pipe_commit_date.strip())
		if (git_pipe_commit_diff != ""):
			print("Git pipeline diff: \t\t" + git_pipe_commit_diff.strip())
		print("Python version:\t\t" + platform.python_version())
		print("Cmd: " + str(" ".join(sys.argv)))
		# Print all arguments
		argsDict = vars(args)
		for arg in argsDict:
			print(arg + ":\t\t" + str(argsDict[arg]))
		print("Run outfolder:\t\t" + self.PIPELINE_OUTFOLDER)
		print("################################################################################")

		self.set_status_flag("running")

		return paths


	def set_status_flag(self, status):
		prev_status = self.STATUS
		# Remove previous status flag file
		flag_file = self.PIPELINE_OUTFOLDER + "/" + self.PIPELINE_NAME + "_" + prev_status + ".flag"
		try:
			os.remove(os.path.join(flag_file))
		except:
			pass

		# Set new status
		self.STATUS = status
		flag_file = self.PIPELINE_OUTFOLDER + "/" + self.PIPELINE_NAME + "_" + status + ".flag"
		self.create_file(flag_file)

		print ("Change status from " + prev_status + " to " + status)


	def timestamp(self, message):
		message += " (" + strftime("%m-%d %H:%M:%S") + ")"
		message += " elapsed:" + str(self.time_elapsed(self.LAST_TIMESTAMP))
		message += " _TIME_"
		if re.match("^###", message):
			message = "\n" + message + "\n"
		print(message)
		LAST_TIMESTAMP = time()


	def report_result(self, key, value):
		message = key + "\t " + str(value).strip()
		print(message + "\t" + "_RES_")
		with open(self.PIPELINE_STATS, "a") as myfile:
			myfile.write(message + "\n")


	def make_sure_path_exists(self, path):
		try:
			os.makedirs(path)
		except OSError as exception:
			if exception.errno != errno.EEXIST:
				raise


	def wait_for_lock(self, lock_file):
		'''Just sleep until the lock_file does not exist'''
		sleeptime = 5
		first_message_flag = False
		dot_count = 0
		while os.path.isfile(lock_file):
			if first_message_flag is False:
				self.timestamp("Waiting for file lock: " + lock_file)
				first_message_flag = True
			else:
				sys.stdout.write(".")
				dot_count = dot_count + 1
				if dot_count % 60 == 0:
					print ""  # linefeed
			sleep(sleeptime)
			sleeptime = min(sleeptime + 5, 60)

		if first_message_flag:
			self.timestamp("File unlocked.")


	def create_file(self, file):
		''' An older function that could succumb to race conditions'''
		with open(file, 'w') as fout:
			fout.write('')


	def create_file_racefree(self, file):
		'''
		This function will only succeed if this process actually
		creates the file; if the file already exists, it will
		raise an OSError
		'''
		write_lock_flags = os.O_CREAT | os.O_EXCL | os.O_WRONLY
		os.open(file, write_lock_flags)


	def call_lock_nofail(self, cmd, lock_name, folder, output_file=None, shell=False):
		"""
		@pass_failure Should the pipeline bail on a nonzero return from a process? Default:  True
		Nofail can be used to implement non-essential parts of the pipeline; if these processes fail,
		they will not cause the pipeline to bail out.
		"""
		self.call_lock_internal(cmd, lock_name, folder, output_file, shell, pass_failure=False)


	def call_lock(self, cmd, lock_name, folder, output_file=None, shell=False):
		"""
		The primary workhorse function of pypiper. This is the command execution function, which enforces
		file-locking to enable restartability, and multiple pipelines using the same files. The function will
		wait for the file lock if it exists, and not produce new output (by default) if the target output file already
		exists. If the output is to be created, it will first create a lock file to prevent other calls to call_lock
		(for example, in parallel pipelines) from touching the file while it is being created.
		It also record the memory of the process and provides some logging output.
		"""
		self.call_lock_internal(cmd, lock_name, folder, output_file, shell, pass_failure=True)


	def call_lock_internal(self, cmd, lock_name, folder, output_file=None, shell=False, pass_failure=True):
		# Create lock file:
		lock_file = os.path.join(folder, lock_name)
		ret = 0
		local_maxmem = 0

		# The loop here is unlikely to be triggered, and is just a wrapper
		# to prevent race conditions; the lock_file must be created by
		# the current loop. If not, we wait again for it and then
		# re-do the tests.
		while True:
			self.wait_for_lock(lock_file)
			if output_file is None or not (os.path.exists(output_file)):
				try:
					self.create_file_racefree(lock_file)     # Create lock
				except OSError as e:
					if e.errno == errno.EEXIST:  # File already exists
						print ("Lock file created after test! Looping again.")
						continue  # Go back to start
				# If you make it past the try block, we successfully
				# created the lock file and should proceed.
				if output_file is not None:
					print ("File to produce: " + output_file)
				try:
					if isinstance(cmd, list):  # Handle command lists
						for cmd_i in cmd:
							list_ret, list_maxmem = self.callprint(cmd_i, shell)
							local_maxmem = max(local_maxmem, list_maxmem)
							ret = max(ret, list_ret)
					else:  # Single command (most common)
						ret, local_maxmem = self.callprint(cmd, shell)   # Run command
				except Exception as e:
					if (pass_failure):
						self.fail_pipeline(e)
					else:
						print(e)
						print("Process failed, but pipeline is continuing because pass_failure=False")
				os.remove(lock_file)        # Remove lock file
			else:
				if output_file is not None:
					print("File already exists: " + output_file)
			# If you make it to the end of the while loop, you're done
			break

		return ret


	def callprint(self, cmd, shell=False):
		"""
		split the command to use shell=False;
		leave it together to use shell=True; I use False so I can get the PID
		and poll memory use.
		"""

		print(cmd)
		if not shell:
			if ("|" in cmd or ">" in cmd):
				print("Should this command run in a shell instead of directly in a subprocess?")
			cmd = cmd.split()
		# call(cmd, shell=shell) # old way (no memory profiling)

		p = subprocess.Popen(cmd, shell=shell)
		self.RUNNING_SUBPROCESS = p.pid
		local_maxmem = -1
		sleeptime = 5
		while p.poll() == None:
			if not shell:
				local_maxmem = max(local_maxmem, self.memory_usage(p.pid))
				# print ("int.maxmem (pid:" + str(p.pid) + ") " + str(local_maxmem))
			sleep(sleeptime)
			sleeptime = min(sleeptime + 5, 60)

		# set self.maxmem
		self.PEAKMEM = max(self.PEAKMEM, local_maxmem)
		self.RUNNING_SUBPROCESS = None

		info = "Process " + str(p.pid) + " returned: (" + str(p.returncode) + ")."
		info += " Peak memory: (Process: " + str(local_maxmem) + "b;"
		info += " Pipeline: " + str(self.PEAKMEM) + "b)"
		print (info)
		if p.returncode != 0:
			raise Exception("Process returned nonzero result.")
		return [p.returncode, local_maxmem]

	####################################
	# Pipeline termination functions
	###################################


	def signal_term_handler(self, signal, frame):
		"""
		TERM signal handler function: this function is run if the process receives a termination signal (TERM).
		This may be invoked, for example, by SLURM if the job exceeds its memory or time limits.
		It will simply record a message in the log file, stating that the process was terminated, and then
		gracefully fail the pipeline. This is necessary to 1. set the status flag and 2. provide a meaningful
		error message in the tee'd output; if you do not handle this, then the tee process will be terminated
		before the TERM error message, leading to a confusing log file.
		"""
		message = "Got SIGTERM; Failing gracefully..."
		with open(self.LOGFILE, "a") as myfile:
			myfile.write(message + "\n")
		self.fail_pipeline(Exception("SIGTERM"))
		sys.exit(1)


	def signal_int_handler(self, signal, frame):
		"""
		For catching interrupt (Ctrl +C) signals. Fails gracefully.
		"""
		message = "Got SIGINT; Failing gracefully..."
		with open(self.LOGFILE, "a") as myfile:
			myfile.write(message + "\n")
		self.fail_pipeline(Exception("SIGINT"))
		sys.exit(1)


	def exit_handler(self):
		"""
		This function I register with atexit to run whenever the script is completing.
		A catch-all for uncaught exceptions...
		"""
		print("Exit handler")
		if self.RUNNING_SUBPROCESS is not None:
				self.kill_child(self.RUNNING_SUBPROCESS)
		if self.STATUS != "completed":
			self.set_status_flag("failed")


	def kill_child(self, child_pid):
		if child_pid is None:
			pass
		else:
			print("Killing child process " + str(child_pid))
			os.kill(child_pid, signal.SIGTERM)


	def fail_pipeline(self, e):
		"""
		Stops the pipeline gracefully: sets the status flag to failed and skips the
		normal pipeline completion procedure.
		"""
		self.set_status_flag("failed")
		raise e


	def time_elapsed(self, time_since):
		"""Returns the number of seconds that have elapsed since the time_since parameter"""
		return round(time() - time_since, 2)


	def stop_pipeline(self):
		"""Set status flag to completed and record some time and memory statistics to the log file."""
		self.set_status_flag("completed")
		print ("Total elapsed time: " + str(self.time_elapsed(self.STARTTIME)))
		# print ("Peak memory used: " + str(memory_usage()["peak"]) + "kb")
		print ("Peak memory used: " + str(self.PEAKMEM / 1e6) + " GB")
		self.report_result("Success", strftime("%m-%d %H:%M:%S"))
		self.timestamp("### Pipeline completed at: ")


	def memory_usage(self, pid='self', category="peak"):
		"""Memory usage of the current process in kilobytes."""
		# Thanks Martin Geisler:
		status = None
		result = {'peak': 0, 'rss': 0}
		try:
			# This will only work on systems with a /proc file system
			# (like Linux).
			# status = open('/proc/self/status')
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
