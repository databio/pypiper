#!/usr/env python

"""
Pypiper is a python module with two components:
1) the Pypiper class, and
2) other toolkits with functions for more specific pipeline use-cases.
The Pypiper class can be used to create a procedural pipeline in python.
"""

import os, sys, subprocess, errno
import re, glob
import time
import atexit, signal, platform

class Pypiper:
	"""
	Base class for instantiating a Pypiper object.

	:param name: Choose a name for your pipeline; it's used to name the output files, flags, etc.
	:param outfolder: Folder in which to store the results.
	:param args: Optional args object from ArgumentParser; Pypiper will simply record these arguments from your script
	:param overwrite_locks: Advanced debugging options; this will cause Pypiper to ignore locked files, enabling
	you to restart a failed pipeline where it left off. Be careful, though, this invalidates the typical file locking
	that enables multiple pipelines to run in parallel on the same intermediate files.
	:param fresh_start: NOT IMPLEMENTED
	:param multi: Enables running multiple pipelines in one script; or for interactive use. It simply disables the tee
	of the output, so you won't get output logged to a file.
	:param manual_clean: Overrides the pipeline's clean_add() manual parameters, to *never* clean up intermediate files
	automatically. Useful for debugging; all cleanup files are added to the manual cleanup script.
	"""
	def __init__(self, name, outfolder, args=None, multi=False,
					manual_clean=False, recover=False, fresh=False):
		# Params defines the set of options that could be updated via
		# command line args to a pipeline run, that can be forwarded
		# to Pypiper. If any pypiper arguments are passed
		# (via add_pypiper_args()), these will override the constructor 		# defaults for these arguments.

		params = { 'manual_clean':manual_clean, 'recover':recover, 
					'fresh': fresh }
		
		if args is not None:
			# vars(...) transforms ArgumentParser's Namespace into a dict
			params.update(vars(args)) 

		# Define pipeline-level variables to keep track of global state and some pipeline stats
		# Pipeline settings
		self.pipeline_name = name
		self.overwrite_locks = params['recover']
		self.fresh_start = params['fresh']
		self.manual_clean = params['manual_clean']

		# File paths:
		self.pipeline_outfolder = os.path.join(outfolder, '')
		self.pipeline_log_file = self.pipeline_outfolder + self.pipeline_name + "_log.md"
		self.pipeline_stats_file = self.pipeline_outfolder + self.pipeline_name + "_stats.txt"
		self.cleanup_file = self.pipeline_outfolder + self.pipeline_name + "_cleanup.sh"


		# Pipeline status variables
		self.peak_memory = 0         # memory high water mark
		self.starttime = time.time()
		self.last_timestamp = self.starttime  # time of the last call to timestamp()
		self.status = "initializing"
		self.running_subprocess = None
		self.wait = True # turn off for debugging

		# Register handler functions to deal with interrupt and termination signals;
		# If received, we would then clean up properly (set pipeline status to FAIL, etc).
		atexit.register(self.exit_handler)
		signal.signal(signal.SIGINT, self.signal_int_handler)
		signal.signal(signal.SIGTERM, self.signal_term_handler)

		# Pypiper can keep track of intermediate files to clean up at the end
		self.cleanup_list = []
		self.cleanup_list_conditional = []
		self.start_pipeline(args, multi)

	@staticmethod
	def add_pypiper_args(parser):
		"""
		Use this to take an ArgumentParser in your pipeline, and also parse default pypiper arguments
		:param parser: an ArgumentParser object from your pipeline
		:returns: A new ArgumentParser object, with default pypiper arguments added
		"""
		parser.add_argument('-R', '--recover', dest='recover', action='store_true',
						default=False, help='Recover mode, overwrite locks')
		parser.add_argument('-F', '--fresh-start', dest='fresh', action='store_true',
						default=False, help='Fresh start mode, overwrite all')
		parser.add_argument('-D', '--dirty', dest='dirty', action='store_true',
						default=False, help='Make all cleanups manual') #Useful for debugging
		parser.add_argument('-L', '--overwrite_locks', dest='overwrite_locks', action='store_true',
						default=False, help='Overwrite locks') #Useful for debugging
		return(parser)

	def start_pipeline(self, args=None, multi=False):
		"""
		Do some setup, like tee output, print some diagnostics, create temp files.
		You provide only the output directory (used for pipeline stats, log, and status flag files).
		"""
		# Perhaps this could all just be put into __init__, but I just kind of like the idea of a start function
		self.make_sure_path_exists(self.pipeline_outfolder)

		# Mirror every operation on sys.stdout to log file
		if not multi:
			sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)  # Unbuffer output
			# a for append to file
			tee = subprocess.Popen(["tee", "-a", self.pipeline_log_file], stdin=subprocess.PIPE)
			# In case this pipeline process is terminated with SIGTERM, make sure we kill this spawned process as well.
			atexit.register(self.kill_child, tee.pid)
			os.dup2(tee.stdin.fileno(), sys.stdout.fileno())
			os.dup2(tee.stdin.fileno(), sys.stderr.fileno())

		# Record the git version of the pipeline and pypiper used. This gets (if it is in a git repo):
		# dir: the directory where the code is stored
		# hash: the commit id of the last commit in this repo
		# date: the date of the last commit in this repo
		# diff: a summary of any differences in the current (run) version vs. the committed version

		# Wrapped in try blocks so that the code will not fail if the pipeline or pypiper are not git repositories
		gitvars = {}
		try:
			gitvars['pypiper_dir'] = os.path.dirname(os.path.realpath(__file__))
			gitvars['pypiper_hash'] = subprocess.check_output("cd " + os.path.dirname(os.path.realpath(__file__)) + "; git rev-parse --verify HEAD", shell=True)
			gitvars['pypiper_date'] = subprocess.check_output("cd " + os.path.dirname(os.path.realpath(__file__)) + "; git show -s --format=%ai HEAD", shell=True)
			gitvars['pypiper_diff'] = subprocess.check_output("cd " + os.path.dirname(os.path.realpath(__file__)) + "; git diff --shortstat HEAD", shell=True)
		except Exception:
			pass
		try:
			gitvars['pipe_dir'] = os.path.dirname(os.path.realpath(sys.argv[0]))
			gitvars['pipe_hash'] = subprocess.check_output("cd " + os.path.dirname(os.path.realpath(sys.argv[0])) + "; git rev-parse --verify HEAD", shell=True)
			gitvars['pipe_date'] = subprocess.check_output("cd " + os.path.dirname(os.path.realpath(sys.argv[0])) + "; git show -s --format=%ai HEAD", shell=True)
			gitvars['pipe_diff'] = subprocess.check_output("cd " + os.path.dirname(os.path.realpath(sys.argv[0])) + "; git diff --shortstat HEAD", shell=True)
		except Exception:
			pass

		# Print out a header section in the pipeline log:
		print("################################################################################")
		print("##### [Pipeline run code and environment:]")
		print("Command".rjust(20) + ":  " + str(" ".join(sys.argv)))
		print("Compute host".rjust(20) + ":  " + platform.node())
		print("Working dir".rjust(20) + ":  " + os.getcwd())
		print("Outfolder".rjust(20) + ":  " + self.pipeline_outfolder)

		self.timestamp("Pipeline started at".rjust(20) + ":  ")

		print("\n##### [Version log:]")
		print("Python version".rjust(20) + ":  " + platform.python_version())
		try:
			print("Pypiper dir".rjust(20) + ":  " + gitvars['pypiper_dir'].strip())
			print("Pypiper version".rjust(20) + ":  " + gitvars['pypiper_hash'].strip())
			print("Pypiper date".rjust(20) + ":  " + gitvars['pypiper_date'].strip())
			if (gitvars['pypiper_diff'] != ""):
				print("Pypiper diff".rjust(20) + ":  " + gitvars['pypiper_diff'].strip())
		except KeyError:
			# If any of the keys aren't set, that's OK. It just means pypiper isn't being run from a git repo.
			pass

		try:
			print("Pipeline dir".rjust(20) + ":  " + gitvars['pipe_dir'].strip())
			print("Pipeline version".rjust(20) + ":  " + gitvars['pipe_hash'].strip())
			print("Pipeline date:".rjust(20) + ":  " + gitvars['pipe_date'].strip())
			if (gitvars['pipe_diff'] != ""):
				print("Pipeline diff".rjust(20) + ":  " + gitvars['pipe_diff'].strip())
		except KeyError:
			# If any of the keys aren't set, that's OK. It just means the pipeline isn' a git repo.
			pass


		# Print all arguments (if any)
		print("\n##### [Arguments passed to pipeline:]")
		if args is not None:
			argsDict = vars(args)
			for arg in argsDict:
				print("* " + arg.rjust(20) + ":  " + str(argsDict[arg]))
		print("################################################################################")
		self.set_status_flag("running")

	def set_status_flag(self, status):
		prev_status = self.status
		# Remove previous status flag file
		flag_file = self.pipeline_outfolder + "/" + self.pipeline_name + "_" + prev_status + ".flag"
		try:
			os.remove(os.path.join(flag_file))
		except:
			pass

		# Set new status
		self.status = status
		flag_file = self.pipeline_outfolder + "/" + self.pipeline_name + "_" + status + ".flag"
		self.create_file(flag_file)

		print("Change status from " + prev_status + " to " + status)

	###################################
	# Process calling functions
	###################################

	def call_lock(self, cmd, target=None, lock_name=None, shell=False, nofail=False, clean=False):
		"""
		The primary workhorse function of pypiper. This is the command execution function, which enforces
		file-locking, enables restartability, and multiple pipelines can produce/use the same files. The function will
		wait for the file lock if it exists, and not produce new output (by default) if the target output file already
		exists. If the output is to be created, it will first create a lock file to prevent other calls to call_lock
		(for example, in parallel pipelines) from touching the file while it is being created.
		It also records the memory of the process and provides some logging output.

		:param cmd: Bash command(s) to be run.
		:type cmd: str or list
		:param target: Output file to be produced. Optional.
		:type target: str or None
		:param lock_name: Name of lock file. Optional.
		:type lock_name: str or None
		:param shell: If command requires should be run in its own shell. Optional. Default: False.
		:type shell: bool
		:param nofail: Should the pipeline bail on a nonzero return from a process? Default: False
		Nofail can be used to implement non-essential parts of the pipeline; if these processes fail,
		they will not cause the pipeline to bail out.
		:type nofail: bool
		:param clean: True means the files will be automatically added to a auto cleanup list. Optional.
		:type clean: bool
		:returns: Return code of process. If a list of commands is passed, this is the maximum of all return codes for all commands.
		:rtype: int
		"""

		# The default lock name is based on the target name. Therefore, a targetless command that you want
		# to lock must specify a lock_name manually.
		if target is None and lock_name is None:
				raise Exception("You must provide either a target or a lock name.")

		# Create lock file:
		# Default lock_name (if not provided) is based on the target file name,
		# but placed in the parent pipeline outfolder, and not in a subfolder, if any.
		if lock_name is None:
			lock_name = target.replace(self.pipeline_outfolder, "").replace("/", "__")

		# Prepend "lock." to make it easy to find the lock files.
		lock_name = "lock." + lock_name
		lock_file = os.path.join(self.pipeline_outfolder, lock_name)
		process_return_code = 0
		local_maxmem = 0

		# The loop here is unlikely to be triggered, and is just a wrapper
		# to prevent race conditions; the lock_file must be created by
		# the current loop. If not, we wait again for it and then
		# re-do the tests.
        # TODO: maybe output a message if when repeatedly going through the loop
		
		while True:
			##### Tests block
			# Base case: Target exists.	# Scenario 3: Target exists (and we don't overwrite); break loop, don't run process.
			if target is not None and os.path.isfile(target) and not os.path.isfile(lock_file):
				print("Target exists: " + target)
				break # Do not run command

			# Scenario 1: Lock file exists, but we're supposed to overwrite target; Run process.
			if os.path.isfile(lock_file):
				if self.overwrite_locks:
					print("Found lock file; overwriting this target...")
				else: # don't overwite locks
					self.wait_for_lock(lock_file)
					# when it's done loop through again to try one more time (to see if the target exists now)
					continue

			# if you get to this point, the target doesn't exist, and the lock_file doesn't exist (or we should overwrite).
			# create the lock (if you can)
			if not self.overwrite_locks:
				try:
					self.create_file_racefree(lock_file)  # Create lock
				except OSError as e:
					if e.errno == errno.EEXIST:  # File already exists
						print ("Lock file created after test! Looping again.")
						continue  # Go back to start
			else:
				self.create_file(lock_file)

			##### End tests block
			# If you make it past theses tests, we should proceed to run the process.

			if target is not None:
				print("Target to produce: " + target)
			else:
				print("Targetless command, running...")

			if isinstance(cmd, list):  # Handle command lists
				for cmd_i in cmd:
					list_ret, list_maxmem = self.callprint(cmd_i, shell, nofail)
					local_maxmem = max(local_maxmem, list_maxmem)
					process_return_code = max(process_return_code, list_ret)

			else:  # Single command (most common)
				process_return_code, local_maxmem = self.callprint(cmd, shell, nofail)  # Run command

			# For temporary files, you can specify a clean option to automatically add them to the clean list,
			# saving you a manual call to clean_add
			if target is not None and clean:
				self.clean_add(target)

			os.remove(lock_file)  # Remove lock file

			# If you make it to the end of the while loop, you're done
			break

		return process_return_code

	def checkprint(self, cmd, shell=False, nofail=False):
		"""
		A wrapper around subprocess.check_output() that also prints the command,
		and can understand the nofail parameter. Just like callprint, but checks output. This is equivalent to
		running subprocess.check_output() instead of subprocess.call().

		:param cmd: Bash command(s) to be run.
		:type cmd: str or list
		:param shell: If command requires should be run in its own shell. Optional. Default: False.
		:type shell: bool
		:param nofail: Should the pipeline bail on a nonzero return from a process? Default: False
		Nofail can be used to implement non-essential parts of the pipeline; if these processes fail,
		they will not cause the pipeline to bail out.
		:type nofail: bool
		"""
		print(cmd)
		if not shell:
			if ("|" in cmd or ">" in cmd):
				print("Should this command run in a shell instead of directly in a subprocess?")
			cmd = cmd.split()

		try:
			returnvalue = subprocess.check_output(cmd, shell=shell)

		except Exception as e:
			if not nofail:
				self.fail_pipeline(e)
			else:
				print(e)
				print("ERROR: Subprocess returned nonzero result, but pipeline is continuing because nofail=True")

		return returnvalue

	def callprint(self, cmd, shell=False, nofail=False):
		"""
		Prints the command, and then executes it, then prints the memory use and return code of the command.

		Uses python's subprocess.Popen() to execute the given command. The shell argument is simply
		passed along to Popen(). You should use shell=False (default) where possible, because this enables memory
		profiling. You should use shell=True if you require shell functions like redirects (>) or pipes (|), but this
		will prevent the script from monitoring memory use.

		cmd can also be a series (a dict object) of multiple commands, which will be run in succession.

		:param cmd: Bash command(s) to be run.
		:type cmd: str or list
		:param shell: If command requires should be run in its own shell. Optional. Default: False.
		:type shell: bool
		:param nofail: Should the pipeline bail on a nonzero return from a process? Default: False
		Nofail can be used to implement non-essential parts of the pipeline; if these processes fail,
		they will not cause the pipeline to bail out.
		:type nofail: bool
		"""
		# The Popen shell argument works like this:
		# if shell=False, then we format the command (with split()) to be a list of command and its arguments.
		# Split the command to use shell=False;
		# leave it together to use shell=True;
		print(cmd)
		if not shell:
			if ("|" in cmd or ">" in cmd):
				print("Should this command run in a shell instead of directly in a subprocess?")
			cmd = cmd.split()
		# call(cmd, shell=shell) # old way (no memory profiling)

		# Try to execute the command:
		# Putting it in a try block enables us to catch exceptions from bad subprocess
		# commands, as well as from valid commands that just fail

		returncode = -1  # set default return values for failed command
		local_maxmem = -1
		try:
			p = subprocess.Popen(cmd, shell=shell)

			# Keep track of the running process ID in case we need to kill it when the pipeline is interrupted.
			self.running_subprocess = p

			sleeptime = .25

			if not self.wait:
				print ("Not waiting for subprocess: " + str(p.pid))
				return [0, -1]

			while p.poll() is None:
				if not shell:
					local_maxmem = max(local_maxmem, self.memory_usage(p.pid))
					# print("int.maxmem (pid:" + str(p.pid) + ") " + str(local_maxmem))
				time.sleep(sleeptime)
				sleeptime = min(sleeptime + 5, 60)

			returncode = p.returncode
			info = "Process " + str(p.pid) + " returned: (" + str(p.returncode) + ")."
			if not shell:
				info += " Peak memory: (Process: " + str(local_maxmem) + "b;"
				info += " Pipeline: " + str(self.peak_memory) + "b)"

			print(info)
			# set self.maxmem
			self.peak_memory = max(self.peak_memory, local_maxmem)
			self.running_subprocess = None

			if p.returncode != 0:
				raise Exception("Subprocess returned nonzero result.")

		except Exception as e:
			if not nofail:
				self.fail_pipeline(e)
			else:
				print(e)
				print("ERROR: Subprocess returned nonzero result, but pipeline is continuing because nofail=True")

		return [returncode, local_maxmem]

	def wait_for_process(self, p, shell=False):
		"""
		Debug function used in unit tests.

		:param p: A subprocess.Popen process.
		:param shell: If command requires should be run in its own shell. Optional. Default: False.
		:type shell: bool
		"""
		local_maxmem = -1
		sleeptime = .5
		while p.poll() is None:
			if not shell:
				local_maxmem = max(local_maxmem, self.memory_usage(p.pid))
				# print("int.maxmem (pid:" + str(p.pid) + ") " + str(local_maxmem))
			time.sleep(sleeptime)
			sleeptime = min(sleeptime + 5, 60)

		self.peak_memory = max(self.peak_memory, local_maxmem)
		self.running_subprocess = None

		info = "Process " + str(p.pid) + " returned: (" + str(p.returncode) + ")."
		if not shell:
			info += " Peak memory: (Process: " + str(local_maxmem) + "b;"
			info += " Pipeline: " + str(self.peak_memory) + "b)"

		print(info)
		if p.returncode != 0:
			raise Exception("Process returned nonzero result.")
		return [p.returncode, local_maxmem]

	def wait_for_lock(self, lock_file):
		"""
		Just sleep until the lock_file does not exist.

		:param lock_file: Lock file to wait upon.
		:type lock_file: str
		"""
		sleeptime = .5
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
					print("")  # linefeed
			time.sleep(sleeptime)
			sleeptime = min(sleeptime + 2.5, 60)

		if first_message_flag:
			self.timestamp("File unlocked.")

	###################################
	# Logging functions
	###################################

	def timestamp(self, message=""):
		"""
		Prints your given message, along with the current time, and time elapsed since the previous timestamp() call.
		If you specify a HEADING by beginning the message with "###", it surrounds the message with newlines
		for easier readability in the log file.

		:param message: Message to timestamp.
		:type message: str
		"""
		message += " (" + time.strftime("%m-%d %H:%M:%S") + ")"
		message += " elapsed:" + str(self.time_elapsed(self.last_timestamp))
		message += " _TIME_"
		if re.match("^###", message):
			message = "\n" + message + "\n"
		print(message)
		self.last_timestamp = time.time()

	def time_elapsed(self, time_since):
		"""
		Returns the number of seconds that have elapsed since the time_since parameter.

		:param time_since: Time as a float given by time.time().
		:type time_since: float
		"""
		return round(time.time() - time_since, 2)

	def report_result(self, key, value):
		"""
		Writes a string to self.pipeline_stats_file.

		:type key: str
		"""
		message = key + "\t " + str(value).strip()
		print(message + "\t" + "_RES_")
		with open(self.pipeline_stats_file, "a") as myfile:
			myfile.write(message + "\n")

	def create_file(self, file):
		"""
		Creates a file, but will not fail if the file already exists. (Vulnerable to race conditions).
		Use this for cases where it doesn't matter if this process is the one that created the file.

		:param file: File to create.
		:type file: str
		"""
		with open(file, 'w') as fout:
			fout.write('')

	def create_file_racefree(self, file):
		"""
		Creates a file, but fails if the file already exists.
		This function will thus only succeed if this process actually creates the file;
		if the file already exists, it will	raise an OSError, solving race conditions.

		:param file: File to create.
		:type file: str
		"""
		write_lock_flags = os.O_CREAT | os.O_EXCL | os.O_WRONLY
		os.open(file, write_lock_flags)

	def make_sure_path_exists(self, path):
		"""
		Creates all directories in a path if it does not exist. Raise exception otherwise.

		:param path: Path to create.
		:type path: str
		"""
		try:
			os.makedirs(path)
		except OSError as exception:
			if exception.errno != errno.EEXIST:
				raise


	###################################
	# Pipeline termination functions
	###################################

	def stop_pipeline(self):
		"""
		The normal pipeline completion function, to be run by the user at the end of the script.
		It sets status flag to completed and records some time and memory statistics to the log file.
		"""
		self.set_status_flag("completed")
		self.cleanup()
		print("Total elapsed time: " + str(self.time_elapsed(self.starttime)))
		# print("Peak memory used: " + str(memory_usage()["peak"]) + "kb")
		print("Peak memory used: " + str(self.peak_memory / 1e6) + " GB")
		self.report_result("Success", time.strftime("%m-%d %H:%M:%S"))
		self.timestamp("### Pipeline completed at: ")

	def fail_pipeline(self, e):
		"""
		If the pipeline does not complete, this function will stop the pipeline gracefully.
		It sets the status flag to failed and skips the	normal success completion procedure.

		:param e: Exception to raise.
		:type e: Exception
		"""
		self.set_status_flag("failed")
		self.timestamp("### Pipeline failed at: ")
		raise e

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
		with open(self.pipeline_log_file, "a") as myfile:
			myfile.write(message + "\n")
		self.fail_pipeline(Exception("SIGTERM"))
		sys.exit(1)

	def signal_int_handler(self, signal, frame):
		"""
		For catching interrupt (Ctrl +C) signals. Fails gracefully.
		"""
		message = "Got SIGINT (Ctrl +C); Failing gracefully..."
		with open(self.pipeline_log_file, "a") as myfile:
			myfile.write(message + "\n")
		self.fail_pipeline(Exception("SIGINT"))
		sys.exit(1)

	def exit_handler(self):
		"""
		This function I register with atexit to run whenever the script is completing.
		A catch-all for uncaught exceptions, setting status flag file to failed.
		"""
		#print("Exit handler")
		# Make the cleanup file executable if it exists
		if os.path.isfile(self.cleanup_file):
			# Make the cleanup file self destruct.
			with open(self.cleanup_file, "a") as myfile:
				myfile.write("rm " + self.cleanup_file + "\n")
			os.chmod(self.cleanup_file, 0755)

		if self.running_subprocess is not None:
				self.kill_child(self.running_subprocess.pid)
		if self.status != "completed":
			self.set_status_flag("failed")

	def kill_child(self, child_pid):
		"""
		Pypiper spawns subprocesses. We need to kill them to exit gracefully,
		in the event of a pipeline termination or interrupt signal.
		By default, child processes are not automatically killed when python
		terminates, so Pypiper must clean these up manually.
		Given a process ID, this function just kills it.

		:param child_pid: Child process id.
		:type child_pid: int
		"""
		if child_pid is None:
			pass
		else:
			print("Pypiper terminating spawned child process " + str(child_pid))
			os.kill(child_pid, signal.SIGTERM)

	def clean_add(self, regex, conditional=False, manual=False):
		"""
		Add files (or regexs) to a cleanup list, to delete when this pipeline completes successfully.
		When making a call with call_lock that produces intermediate files that should be
		deleted after the pipeline completes, you flag these files for deletion with this command.
		Files added with clean_add will only be deleted upon success of the pipeline.

		:param regex:  A unix-style regular expression that matches files to delete
		(can also be a file name).
		:type regex: str
		:param conditional: True means the files will only be deleted if no other
		pipelines are currently running; otherwise they are added to a manual cleanup script
		called {pipeline_name}_cleanup.sh
		:type conditional: bool
		:param manual: True means the files will just be added to a manual cleanup script.
		:type manual: bool
		"""
		if self.manual_clean:
			# Override the user-provided option and force manual cleanup.
			manual = True

		if manual:
			try:
				files = glob.glob(regex)
				for file in files:
					with open(self.cleanup_file, "a") as myfile:
						myfile.write("rm " + file + "\n")
			except:
				pass
		elif conditional:
			self.cleanup_list_conditional.append(regex)
		else:
			self.cleanup_list.append(regex)
			# Remove it from the conditional list if added to the absolute list
			while regex in self.cleanup_list_conditional: self.cleanup_list_conditional.remove(regex)

	def cleanup(self):
		"""
		Cleans up (removes) intermediate files. You can register intermediate files, which will
		be deleted automatically when the pipeline completes. This function deletes them, either
		absolutely or conditionally.
		"""
		if len(self.cleanup_list) > 0:
			print("Cleaning up flagged intermediate files...")
			for expr in self.cleanup_list:
				print("Removing glob: " + expr)
				try:
					# Expand regular expression
					files = glob.glob(expr)
					# Remove entry from cleanup list
					while files in self.cleanup_list: self.cleanup_list.remove(files)
					# and delete the files
					for file in files:
						print("rm " + file)
						os.remove(os.path.join(file))
				except:
					pass

		if len(self.cleanup_list_conditional) > 0:
			# flag_files = glob.glob(self.pipeline_outfolder + "*.flag")
			flag_files = [fn for fn in glob.glob(self.pipeline_outfolder + "*.flag") if not "completed" in os.path.basename(fn) and not self.pipeline_name + "_running.flag" == os.path.basename(fn)]
			if (len(flag_files) == 0):
				print("Cleaning up conditional list...")
				for expr in self.cleanup_list_conditional:
					print("Removing glob: " + expr)
					try:
						files = glob.glob(expr)
						while files in self.cleanup_list_conditional: self.cleanup_list_conditional.remove(files)
						for file in files:
							print("rm " + file)
							os.remove(os.path.join(file))
					except:
						pass
			else:
				print("Conditional flag found: " + str([os.path.basename(i) for i in flag_files]))
				print("These conditional files were left in place:" + str(self.cleanup_list_conditional))
				# Produce a cleanup script
				for expr in self.cleanup_list_conditional:
					try:
						files = glob.glob(expr)
						for file in files:
							with open(self.cleanup_file, "a") as myfile:
								myfile.write("rm " + file + "\n")
					except:
						pass

	def memory_usage(self, pid='self', category="peak"):
		"""
		Memory usage of the current process in kilobytes.

		:param pid:
		:type pid: str
		:param category:
		:type category: str
		"""
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
