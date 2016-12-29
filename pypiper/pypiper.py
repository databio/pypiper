#!/usr/env python

"""
Pypiper is a python module with two components:
1) the PipelineManager class, and
2) other toolkits (currently just NGSTk) with functions for more specific pipeline use-cases.
The PipelineManager class can be used to create a procedural pipeline in python.
"""

import os, sys, subprocess, errno
import re, glob
import time
import atexit, signal, platform
import __main__
from AttributeDict import AttributeDict

import shlex # for splitting commands like a shell does
import datetime

class PipelineManager(object):
	"""
	Base class for instantiating a PipelineManager object, the main class of Pypiper

	:param name: Choose a name for your pipeline; it's used to name the output files, flags, etc.
	:param outfolder: Folder in which to store the results.
	:param args: Optional args object from ArgumentParser; Pypiper will simply record these arguments from your script
	:param overwrite_locks: Advanced debugging options; this will cause Pypiper to ignore locked files, enabling
		you to restart a failed pipeline where it left off. Be careful, though, this invalidates the typical file locking
		that enables multiple pipelines to run in parallel on the same intermediate files.
	:param fresh: NOT IMPLEMENTED
	:param follow: Force run all follow functions, even if the preceeding command is not run.
		By default, follow functions are only run if the preceeding command is run.
	:param recover: Specify recover mode, to overwrite lock files. If pypiper encounters a locked
		target, it will ignore the lock, and recompute this step. Useful to restart a failed pipeline.
	:param multi: Enables running multiple pipelines in one script; or for interactive use. It simply disables the tee
		of the output, so you won't get output logged to a file.
	:param manual_clean: Overrides the pipeline's clean_add() manual parameters, to *never* clean up intermediate files
		automatically. Useful for debugging; all cleanup files are added to the manual cleanup script.
	"""
	def __init__(
		self, name, outfolder, args=None, multi=False,
		manual_clean=False, recover=False, fresh=False,
		force_follow=False,
		cores=1, mem="1000",
		config_file=None, output_parent=None):
		# Params defines the set of options that could be updated via
		# command line args to a pipeline run, that can be forwarded
		# to Pypiper. If any pypiper arguments are passed
		# (via add_pypiper_args()), these will override the constructor
		# defaults for these arguments.

		# Establish default params
		params = {
			'manual_clean': manual_clean,
			'recover': recover,
			'fresh': fresh,
			'force_follow': force_follow,
			'config_file': config_file,
			'output_parent': output_parent,
			'cores': cores,
			'mem': mem}

		# Update them with any passed via 'args'
		if args is not None:
			# vars(...) transforms ArgumentParser's Namespace into a dict
			params.update(vars(args))

		# Define pipeline-level variables to keep track of global state and some pipeline stats
		# Pipeline settings
		self.pipeline_name = name
		self.overwrite_locks = params['recover']
		self.fresh_start = params['fresh']
		self.force_follow = params['force_follow']
		self.manual_clean = params['manual_clean']
		self.cores = params['cores']
		self.output_parent = params['output_parent']
		# For now, we assume the memory is in megabytes.
		# this could become customizable if necessary
		self.mem = params['mem'] + "m"

		# We use this memory to pass a memory limit to processes like java that
		# can take a memory limit, so they don't get killed by a SLURM (or other
		# cluster manager) overage. However, with java, the -Xmx argument can only
		# limit the *heap* space, not total memory use; so occasionally SLURM will
		# still kill these processes because total memory goes over the limit.
		# As a kind of hack, we'll set the java processes heap limit to 95% of the
		# total memory limit provided.
		# This will give a little breathing room for non-heap java memory use.
		self.javamem = str(int(int(params['mem']) * 0.95)) + "m"

		# Set relative output_parent directory to absolute
		# not necessary after all...
		#if self.output_parent and not os.path.isabs(self.output_parent):
		#	self.output_parent = os.path.join(os.getcwd(), self.output_parent)

		# File paths:
		self.pipeline_outfolder = os.path.join(outfolder, '')
		self.pipeline_log_file = self.pipeline_outfolder + self.pipeline_name + "_log.md"
		self.pipeline_profile_file = self.pipeline_outfolder + self.pipeline_name + "_profile.tsv"
		#self.pipeline_stats_file = self.pipeline_outfolder + self.pipeline_name + "_stats.tsv"
		self.pipeline_stats_file = self.pipeline_outfolder + "stats.tsv"
		self.pipeline_commands_file = self.pipeline_outfolder + self.pipeline_name + "_commands.sh"
		self.cleanup_file = self.pipeline_outfolder + self.pipeline_name + "_cleanup.sh"

		# Pipeline status variables
		self.peak_memory = 0  # memory high water mark
		self.starttime = time.time()
		self.last_timestamp = self.starttime  # time of the last call to timestamp()
		self.status = "initializing"
		self.proc_name = None
		self.proc_lock_name = None
		self.proc_start_time = None
		self.running_subprocess = None
		self.wait = True  # turn off for debugging

		# In-memory holder for report_results
		self.stats_dict = {}

		# Register handler functions to deal with interrupt and termination signals;
		# If received, we would then clean up properly (set pipeline status to FAIL, etc).
		atexit.register(self._exit_handler)
		signal.signal(signal.SIGINT, self._signal_int_handler)
		signal.signal(signal.SIGTERM, self._signal_term_handler)

		# Pypiper can keep track of intermediate files to clean up at the end
		self.cleanup_list = []
		self.cleanup_list_conditional = []
		self.start_pipeline(args, multi)

		# Handle config file if it exists

		# Read YAML config file
		# TODO: This section should become a function, so toolkits can use it
		# to locate a config file.
		config_to_load = None  # start with nothing

		if args and args.config_file is not None:
			if os.path.isabs(args.config_file):
				# Absolute custom config file specified
				if os.path.isfile(args.config_file):
					config_to_load = args.config_file
				else:
					print("Can't find custom config file: " + args.config_file)
			else: 
				# Relative custom config file specified
				# Set path to be relative to pipeline script
				abs_config = os.path.join(os.path.dirname(sys.argv[0]), args.config_file)
				if os.path.isfile(abs_config):
					config_to_load = abs_config
				else:
					print("Can't find custom config file: " + abs_config)
		else:
			# No custom config file specified. Check for default
			default_config = os.path.splitext(os.path.basename(sys.argv[0]))[0] + ".yaml"
			if os.path.isfile(default_config):
				config_to_load = default_config

		# Finally load the config we found.
		if config_to_load is not None:
			with open(config_to_load, 'r') as config_file:
				# Set the args to the new config file, so it can be used
				# later to pass to, for example, toolkits
				args.config_file = config_to_load
				import yaml
				config = yaml.load(config_file)
				self.config = AttributeDict(config, default=True)
		else:
			self.config = None


	def _ignore_interrupts(self):
		"""
		Ignore interrupt and termination signals. Used as a pre-execution
		function (preexec_fn) for subprocess.Popen calls that Pyper will retain
		control over (meaning I will clean these processes up manually).
		"""
		signal.signal(signal.SIGINT, signal.SIG_IGN)
		signal.signal(signal.SIGTERM, signal.SIG_IGN)

	def start_pipeline(self, args = None, multi = False):
		"""
		Initialize setup. Do some setup, like tee output, print some diagnostics, create temp files.
		You provide only the output directory (used for pipeline stats, log, and status flag files).
		"""
		# Perhaps this could all just be put into __init__, but I just kind of like the idea of a start function
		self.make_sure_path_exists(self.pipeline_outfolder)

		# By default, Pypiper will mirror every operation so it is displayed both
		# on sys.stdout **and** to a log file. Unfortunately, interactive python sessions
		# ruin this by interfering with stdout. So, for interactive mode, we do not enable 
		# the tee subprocess, sending all output to screen only.
		# Starting multiple PipelineManagers in the same script has the same problem, and
		# must therefore be run in interactive_mode.
		
		interactive_mode = False  # Defaults to off

		if multi or not hasattr(__main__, "__file__"):
			print("Warning: You're running an interactive python session. This works, but pypiper cannot tee\
					the output, so results are only logged to screen.")

			interactive_mode = True
		
		if not interactive_mode:
			sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)  # Unbuffer output

			# The tee subprocess must be instructed to ignore TERM and INT signals;
			# Instead, I will clean up this process in the signal handler functions.
			# This is required because otherwise, if pypiper recieves a TERM or INT,
			# the tee will be automatically terminated by python before I have a chance to
			# print some final output (for example, about when the process stopped),
			# and so those things don't end up in the log files because the tee
			# subprocess is dead. Instead, I will handle the killing of the tee process
			# manually (in the exit handler).

			# a for append to file
			
			tee = subprocess.Popen(
				["tee", "-a", self.pipeline_log_file], stdin=subprocess.PIPE,
				preexec_fn=self._ignore_interrupts)

			# If the pipeline is terminated with SIGTERM/SIGINT,
			# make sure we kill this spawned tee subprocess as well.
			atexit.register(self._kill_child_process, tee.pid)
			os.dup2(tee.stdin.fileno(), sys.stdout.fileno())
			os.dup2(tee.stdin.fileno(), sys.stderr.fileno())

		# A future possibility to avoid this tee, is to use a Tee class; this works for anything printed here
		# by pypiper, but can't tee the subprocess output. For this, it would require using threading to
		# simultaneously capture and display subprocess output. I shelve this for now and stick with the tee option.
		# sys.stdout = Tee(self.pipeline_log_file)

		# Record the git version of the pipeline and pypiper used. This gets (if it is in a git repo):
		# dir: the directory where the code is stored
		# hash: the commit id of the last commit in this repo
		# date: the date of the last commit in this repo
		# diff: a summary of any differences in the current (run) version vs. the committed version

		# Wrapped in try blocks so that the code will not fail if the pipeline or pypiper are not git repositories
		gitvars = {}
		try:
			gitvars['pypiper_dir'] = os.path.dirname(os.path.realpath(__file__))
			gitvars['pypiper_hash'] = subprocess.check_output("cd " + os.path.dirname(os.path.realpath(__file__)) + "; git rev-parse --verify HEAD 2>/dev/null", shell=True)
			gitvars['pypiper_date'] = subprocess.check_output("cd " + os.path.dirname(os.path.realpath(__file__)) + "; git show -s --format=%ai HEAD 2>/dev/null", shell=True)
			gitvars['pypiper_diff'] = subprocess.check_output("cd " + os.path.dirname(os.path.realpath(__file__)) + "; git diff --shortstat HEAD 2>/dev/null", shell=True)
			gitvars['pypiper_branch'] = subprocess.check_output("cd " + os.path.dirname(os.path.realpath(__file__)) + "; git branch | grep '*' 2>/dev/null", shell=True)
		except Exception:
			pass
		try:
			gitvars['pipe_dir'] = os.path.dirname(os.path.realpath(sys.argv[0]))
			gitvars['pipe_hash'] = subprocess.check_output("cd " + os.path.dirname(os.path.realpath(sys.argv[0])) + "; git rev-parse --verify HEAD 2>/dev/null", shell=True)
			gitvars['pipe_date'] = subprocess.check_output("cd " + os.path.dirname(os.path.realpath(sys.argv[0])) + "; git show -s --format=%ai HEAD 2>/dev/null", shell=True)
			gitvars['pipe_diff'] = subprocess.check_output("cd " + os.path.dirname(os.path.realpath(sys.argv[0])) + "; git diff --shortstat HEAD 2>/dev/null", shell=True)
			gitvars['pipe_branch'] = subprocess.check_output("cd " + os.path.dirname(os.path.realpath(sys.argv[0])) + "; git branch | grep '*' 2>/dev/null", shell=True)
		except Exception:
			pass
		
		# Print out a header section in the pipeline log:
		# Wrap things in backticks to prevent markdown from interpreting underscores as emphasis.
		print("----------------------------------------")
		print("##### [Pipeline run code and environment:]")
		print("* " + "Command".rjust(20) + ":  " + "`" + str(" ".join(sys.argv))) + "`"
		print("* " + "Compute host".rjust(20) + ":  " + platform.node())
		print("* " + "Working dir".rjust(20) + ":  " + os.getcwd())
		print("* " + "Outfolder".rjust(20) + ":  " + self.pipeline_outfolder)

		self.timestamp("* " + "Pipeline started at".rjust(20) + ":  ")

		print("\n##### [Version log:]")
		print("* " + "Python version".rjust(20) + ":  " + platform.python_version())
		try:
			print("* " + "Pypiper dir".rjust(20) + ":  " + "`" + gitvars['pypiper_dir'].strip() + "`")
			print("* " + "Pypiper version".rjust(20) + ":  " + gitvars['pypiper_hash'].strip())
			print("* " + "Pypiper branch".rjust(20) + ":  " + gitvars['pypiper_branch'].strip())
			print("* " + "Pypiper date".rjust(20) + ":  " + gitvars['pypiper_date'].strip())
			if (gitvars['pypiper_diff'] != ""):
				print("* " + "Pypiper diff".rjust(20) + ":  " + gitvars['pypiper_diff'].strip())
		except KeyError:
			# If any of the keys aren't set, that's OK. It just means pypiper isn't being run from a git repo.
			pass

		try:
			print("* " + "Pipeline dir".rjust(20) + ":  " + "`" + gitvars['pipe_dir'].strip() + "`")
			print("* " + "Pipeline version".rjust(20) + ":  " + gitvars['pipe_hash'].strip())
			print("* " + "Pipeline branch".rjust(20) + ":  " + gitvars['pipe_branch'].strip())
			print("* " + "Pipeline date".rjust(20) + ":  " + gitvars['pipe_date'].strip())
			if (gitvars['pipe_diff'] != ""):
				print("* " + "Pipeline diff".rjust(20) + ":  " + gitvars['pipe_diff'].strip())
		except KeyError:
			# It is ok if keys aren't set, it means the pipeline isn't a git repo.
			pass

		# Print all arguments (if any)
		print("\n##### [Arguments passed to pipeline:]")
		if args is not None:
			argsDict = vars(args)
			for arg in argsDict:
				arg_print = "`" + arg + "`"
				print("* " + arg_print.rjust(20) + ":  " + "`" + str(argsDict[arg]) + "`")
		print("\n----------------------------------------\n")
		self.set_status_flag("running")

		# Record the start in PIPE_profile and PIPE_commands output files so we
		# can trace which run they belong to

		with open(self.pipeline_commands_file, "a") as myfile:
			myfile.write("# Pipeline started at " + time.strftime("%m-%d %H:%M:%S", time.localtime(self.starttime)) + "\n\n")

		with open(self.pipeline_profile_file, "a") as myfile:
			myfile.write("# Pipeline started at " + time.strftime("%m-%d %H:%M:%S", time.localtime(self.starttime)) + "\n\n")

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
		self._create_file(flag_file)

		print("\nChange status from " + prev_status + " to " + status)

	###################################
	# Process calling functions
	###################################

	def run(self, cmd, target=None, lock_name=None, shell="guess", nofail=False, clean=False, follow=None):
		"""
		Runs a shell command. This is the primary workhorse function of PipelineManager. This is the command  execution
		function, which enforces racefree file-locking, enables restartability, and multiple pipelines can produce/use
		the same files. The function will wait for the file lock if it exists, and not produce new output (by default)
		if the target output file already exists. If the output is to be created, it will first create a lock file to
		prevent other calls to run (for example, in parallel pipelines) from touching the file while it is being
		created. It also records the memory of the process and provides some logging output

		:param cmd: Shell command(s) to be run.
		:type cmd: str or list
		:param target: Output file to be produced. Optional.
		:type target: str or None
		:param lock_name: Name of lock file. Optional.
		:type lock_name: str or None
		:param shell: If command requires should be run in its own shell. Optional. Default: "guess" -- run will try to determine if the command requires a shell.
		:type shell: bool
		:param nofail: Should the pipeline bail on a nonzero return from a process? Default: False
			Nofail can be used to implement non-essential parts of the pipeline; if these processes fail,
			they will not cause the pipeline to bail out.
		:type nofail: bool
		:param clean: True means the target file will be automatically added to a auto cleanup list. Optional.
		:type clean: bool
		:returns: Return code of process. If a list of commands is passed, this is the maximum of all return codes for all commands.
		:rtype: int
		"""

		# The default lock name is based on the target name. Therefore, a targetless command that you want
		# to lock must specify a lock_name manually.
		if target is None and lock_name is None:
				raise Exception("You must provide either a target or a lock_name.")

		# If the target is a list, for now let's just strip it to the first target.
		# Really, it should just check for all of them.
		if type(target) == list:
			target = target[0]
			#primary_target = target[0]
		# Create lock file:
		# Default lock_name (if not provided) is based on the target file name,
		# but placed in the parent pipeline outfolder, and not in a subfolder, if any.
		if lock_name is None:
			lock_name = target.replace(self.pipeline_outfolder, "").replace("/", "__")

		# Prepend "lock." to make it easy to find the lock files.
		self.proc_lock_name = lock_name
		lock_name = "lock." + lock_name
		lock_file = os.path.join(self.pipeline_outfolder, lock_name)
		process_return_code = 0
		local_maxmem = 0

		# The loop here is unlikely to be triggered, and is just a wrapper
		# to prevent race conditions; the lock_file must be created by
		# the current loop. If not, we wait again for it and then
		# re-do the tests.
		# TODO: maybe output a message if when repeatedly going through the loop

		follow_result = None

		while True:
			##### Tests block
			# Base case: Target exists (and we don't overwrite); break loop, don't run process.
			# os.path.exists allows the target to be either a file or directory; .isfile is file-only
			if target is not None and os.path.exists(target) and not os.path.isfile(lock_file):
				print("\nTarget exists: `" + target + "`")

				# Normally we don't run the follow, but if you want to force...
				if self.force_follow and hasattr(follow, '__call__'):
					# Run the follow function
					print("Follow:")
					follow_result = follow()
					
				break  # Do not run command

			# Scenario 1: Lock file exists, but we're supposed to overwrite target; Run process.
			if os.path.isfile(lock_file):
				if self.overwrite_locks:
					print("Found lock file; overwriting this target...")
				else:  # don't overwite locks
					self._wait_for_lock(lock_file)
					# when it's done loop through again to try one more time (to see if the target exists now)
					continue

			# If you get to this point, the target doesn't exist, and the lock_file doesn't exist 
			# (or we should overwrite). create the lock (if you can)
			if not self.overwrite_locks:
				try:
					self._create_file_racefree(lock_file)  # Create lock
				except OSError as e:
					if e.errno == errno.EEXIST:  # File already exists
						print ("Lock file created after test! Looping again.")
						continue  # Go back to start
			else:
				self._create_file(lock_file)

			##### End tests block
			# If you make it past these tests, we should proceed to run the process.

			if target is not None:
				print("\nTarget to produce: `" + target + "`")
			else:
				print("\nTargetless command, running...")

			if isinstance(cmd, list):  # Handle command lists
				for cmd_i in cmd:
					list_ret, list_maxmem = self.callprint(cmd_i, shell, nofail)
					local_maxmem = max(local_maxmem, list_maxmem)
					process_return_code = max(process_return_code, list_ret)

			else:  # Single command (most common)
				process_return_code, local_maxmem = self.callprint(cmd, shell, nofail)  # Run command

			# For temporary files, you can specify a clean option to automatically
			# add them to the clean list, saving you a manual call to clean_add
			if target is not None and clean:
				self.clean_add(target)

			if hasattr(follow, '__call__'):
				# Run the follow function
				print("Follow:")
				follow_result = follow()

			os.remove(lock_file)  # Remove lock file

			# If you make it to the end of the while loop, you're done
			break

		# Bad idea: don't return follow_result; it seems nice but nothing else
		# in your pipeline can depend on this since it won't be run if that command 		# isn't required because target exists.
		return process_return_code

	def checkprint(self, cmd, shell="guess", nofail=False):
		"""
		Just like callprint, but checks output -- so you can get a variable
		in python corresponding to the return value of the command you call.
		This is equivalent to running subprocess.check_output() 
		instead of subprocess.call().

		:param cmd: Bash command(s) to be run.
		:type cmd: str or list
		:param shell: If command requires should be run in its own shell. Optional. 		Default: "guess" -- `run()` will 
		try to guess if the command should be run in a shell (based on the 
		presence of a pipe (|) or redirect (>), To force a process to run as
		a direct subprocess, set `shell` to False; to force a shell, set True.
		:type shell: bool
		:param nofail: Should the pipeline bail on a nonzero return from a process? Default: False
			Nofail can be used to implement non-essential parts of the pipeline; if these processes fail,
			they will not cause the pipeline to bail out.
		:type nofail: bool
		"""
		self._report_command(cmd)


		if shell == "guess":
			if ("|" in cmd or ">" in cmd or r"*" in cmd):
				shell = True
			else:
				shell = False

		if not shell:
			if ("|" in cmd or ">" in cmd or r"*" in cmd):
				print("Should this command run in a shell instead of directly in a subprocess?")
		
			#cmd = cmd.split()
			cmd = shlex.split(cmd)
		# else: # if shell: # do nothing (cmd is not split)
			

		try:
			returnvalue = subprocess.check_output(cmd, shell=shell)

		except Exception as e:
			if not nofail:
				self.fail_pipeline(e)
			else:
				print(e)
				print("ERROR: Subprocess returned nonzero result, but pipeline is continuing because nofail=True")

		return returnvalue

	def callprint(self, cmd, shell="guess", nofail=False):
		"""
		Prints the command, and then executes it, then prints the memory use and
		return code of the command.

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
		self._report_command(cmd)
		# self.proc_name = cmd[0] + " " + cmd[1]
		self.proc_name = "".join(cmd).split()[0]

		if shell == "guess":
			if ("|" in cmd or ">" in cmd or "*" in cmd):
				shell = True
			else:
				shell = False

		if not shell:
			if ("|" in cmd or ">" in cmd or "*" in cmd):
				print("Should this command run in a shell instead of directly in a subprocess?")
			#cmd = cmd.split()
			cmd = shlex.split(cmd)
		# call(cmd, shell=shell) # old way (no memory profiling)

		# Try to execute the command:
		# Putting it in a try block enables us to catch exceptions from bad subprocess
		# commands, as well as from valid commands that just fail

		returncode = -1  # set default return values for failed command
		local_maxmem = -1

		try:
			# Capture the subprocess output in <pre> tags to make it format nicely
			# if the markdown log file is displayed as HTML.
			print("<pre>")
			self.proc_start_time = time.time()
			p = subprocess.Popen(cmd, shell=shell)

			# Keep track of the running process ID in case we need to kill it when the pipeline is interrupted.
			self.running_subprocess = p

			sleeptime = .25

			if not self.wait:
				print("</pre>")
				print ("Not waiting for subprocess: " + str(p.pid))
				return [0, -1]

			while p.poll() is None:
				if not shell:
					local_maxmem = max(local_maxmem, self._memory_usage(p.pid)/1e6)
					# print("int.maxmem (pid:" + str(p.pid) + ") " + str(local_maxmem))
				time.sleep(sleeptime)
				sleeptime = min(sleeptime + 5, 60)

			returncode = p.returncode
			info = "Process " + str(p.pid) + " returned: (" + str(p.returncode) + ")."
			info += " Elapsed: " + str(datetime.timedelta(seconds=self.time_elapsed(self.proc_start_time))) + "."
			if not shell:
				info += " Peak memory: (Process: " + str(round(local_maxmem, 3)) + "GB;"
				info += " Pipeline: " + str(round(self.peak_memory, 3)) + "GB)"
			# Close the preformat tag for markdown output
			print("</pre>")
			print(info)
			# set self.maxmem
			self.peak_memory = max(self.peak_memory, local_maxmem)
			# report process profile
			self._report_profile(self.proc_name, self.proc_lock_name, time.time() - self.proc_start_time, local_maxmem)
			self.running_subprocess = None
			self.proc_lock_name = None
			self.proc_start_time = None

			if p.returncode != 0:
				raise Exception("Subprocess returned nonzero result.")

		except Exception as e:
			if not nofail:
				self.fail_pipeline(e)
			else:
				print(e)
				print("ERROR: Subprocess returned nonzero result, but pipeline is continuing because nofail=True")

		return [returncode, local_maxmem]

	def _wait_for_process(self, p, shell=False):
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
				local_maxmem = max(local_maxmem, self._memory_usage(p.pid) / 1e6)
				# print("int.maxmem (pid:" + str(p.pid) + ") " + str(local_maxmem))
			time.sleep(sleeptime)
			sleeptime = min(sleeptime + 5, 60)

		self.peak_memory = max(self.peak_memory, local_maxmem)
		self.running_subprocess = None

		info = "Process " + str(p.pid) + " returned: (" + str(p.returncode) + ")."
		if not shell:
			info += " Peak memory: (Process: " + str(round(local_maxmem,3)) + "GB;"
			info += " Pipeline: " + str(round(self.peak_memory,3)) + "GB)"

		print(info + "\n")
		if p.returncode != 0:
			raise Exception("Process returned nonzero result.")
		return [p.returncode, local_maxmem]

	def _wait_for_lock(self, lock_file):
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

	def _wait_for_file(self, file_name, lock_name = None):
		"""
		Just sleep until the file_name DOES exist.

		:param file_name: File to wait for.
		:type file_name: str
		"""
		# Build default lock name:
		if lock_name is None:
			lock_name = file_name.replace(self.pipeline_outfolder, "").replace("/", "__")

		lock_name = "lock." + lock_name
		lock_file = os.path.join(self.pipeline_outfolder, lock_name)


		sleeptime = .5
		first_message_flag = False
		dot_count = 0
		while not os.path.isfile(file_name):
			if first_message_flag is False:
				self.timestamp("Waiting for file: " + file_name)
				first_message_flag = True
			else:
				sys.stdout.write(".")
				dot_count = dot_count + 1
				if dot_count % 60 == 0:
					print("")  # linefeed
			time.sleep(sleeptime)
			sleeptime = min(sleeptime + 2.5, 60)

		if first_message_flag:
			self.timestamp("File exists.")

		self._wait_for_lock(lock_file)

	###################################
	# Logging functions
	###################################

	def timestamp(self, message=""):
		"""
		Prints your given message, along with the current time, and time elapsed since the 
		previous timestamp() call.	If you specify a HEADING by beginning the message with "###",
		it surrounds the message with newlines for easier readability in the log file.

		:param message: Message to timestamp.
		:type message: str
		"""
		message += " (" + time.strftime("%m-%d %H:%M:%S") + ")"
		message += " elapsed:" + str(datetime.timedelta(seconds=self.time_elapsed(self.last_timestamp)))
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
		return round(time.time() - time_since, 0)

	def _report_profile(self, command, lock_name, elapsed_time, memory):
		"""
		Writes a string to self.pipeline_profile_file.

		:type key: str
		"""
		messageRaw = str(command) + "\t " + \
			str(lock_name) + "\t" + \
			str(datetime.timedelta(seconds = round(elapsed_time, 2))) + "\t " + \
			str(memory)
		# messageMarkdown = "> `" + command + "`\t" + str(elapsed_time).strip() + "\t " + str(memory).strip() + "\t" + "_PROF_"
		# print(messageMarkdown)
		with open(self.pipeline_profile_file, "a") as myfile:
			myfile.write(messageRaw + "\n")

	def report_result(self, key, value, annotation=None):
		"""
		Writes a string to self.pipeline_stats_file.

		:type key: str
		:param annotation: By default, the stats will be annotated with the pipeline
		name, so you can tell which pipeline records which stats. If you want, you can
		change this; use annotation='shared' if you need the stat to be used by
		another pipeline (using get_stat()).
		:type annotation: str
		"""
		# Default annotation is current pipeline name
		if not annotation:
			annotation = self.pipeline_name

		# keep the value in memory:
		self.stats_dict[key] = str(value).strip()
		messageRaw = key + "\t " + str(value).strip() + "\t" + str(annotation)
		messageMarkdown = "> `" + key + "`\t" + str(value).strip()\
		 + "\t" + str(annotation) + "\t" + "_RES_"
		print(messageMarkdown)

		# Just to be extra careful, let's lock the file while we we write
		# in case multiple pipelines write to the same file.
		target = self.pipeline_stats_file
		lock_name = target.replace(self.pipeline_outfolder, "").replace("/", "__")
		lock_name = "lock." + lock_name
		lock_file = os.path.join(self.pipeline_outfolder, lock_name)

		while True:
			if os.path.isfile(lock_file):
				self._wait_for_lock(lock_file)
			else:
				try:
					self._create_file_racefree(lock_file)  # Create lock
				except OSError as e:
					if e.errno == errno.EEXIST:  # File already exists
						print ("Lock file created after test! Looping again.")
						continue  # Go back to start

				# Proceed with file writing
				with open(self.pipeline_stats_file, "a") as myfile:
					myfile.write(messageRaw + "\n")

				os.remove(lock_file)  # Remove lock file

				# If you make it to the end of the while loop, you're done
				break

	def _report_command(self, cmd):
		"""
		Writes a command to self.pipeline_commands_file.

		:type cmd: str
		"""
		print("> `" + cmd + "`\n")
		with open(self.pipeline_commands_file, "a") as myfile:
			myfile.write(cmd + "\n\n")

	def _create_file(self, file):
		"""
		Creates a file, but will not fail if the file already exists. (Vulnerable to race conditions).
		Use this for cases where it doesn't matter if this process is the one that created the file.

		:param file: File to create.
		:type file: str
		"""
		with open(file, 'w') as fout:
			fout.write('')

	def _create_file_racefree(self, file):
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
	# Pipeline stats calculation helpers
	###################################


	def _refresh_stats(self):
		"""
		Loads up the stats sheet created for this pipeline run and reads
		those stats into memory
		"""

		# regex identifies all possible stats files.
		#regex = self.pipeline_outfolder +  "*_stats.tsv"		
		#stats_files = glob.glob(regex)
		#stats_files.insert(self.pipeline_stats_file) # last one is the current pipeline
		#for stats_file in stats_files:

		stats_file = self.pipeline_stats_file
		if os.path.isfile(self.pipeline_stats_file):
			with open(stats_file, "rb") as stat_file:
				for line in stat_file:
					try:
						# Someone may have put something that's not 3 columns in the stats file
						# if so, shame on him, but we can just ignore it.
						key, value, annotation  = line.split('\t')
					except ValueError:
						print("WARNING: Each row in a stats file is expected to have 3 columns")

					if annotation.rstrip() == self.pipeline_name or annotation.rstrip() == "shared":
						self.stats_dict[key] = value.strip()


		#if os.path.isfile(self.pipeline_stats_file):

	def get_stat(self, key):
		"""
		Returns a stat that was previously reported. This is necessary for reporting new stats that are
		derived from two stats, one of which may have been reported by an earlier run. For example, 
		if you first use report_result to report (number of trimmed reads), and then in a later stage
		want to report alignment rate, then this second stat (alignment rate) will require knowing the 
		first stat (number of trimmed reads); however, that may not have been calculated in the current
		pipeline run, so we must retrieve it from the stats.tsv output file. This command will retrieve
		such previously reported stats if they were not already calculated in the current pipeline run.
		:param key: key of stat to retrieve		
		"""

		if self.stats_dict.has_key(key):
			return (self.stats_dict[key])
		else:
			self._refresh_stats()
			if self.stats_dict.has_key(key):
				return (self.stats_dict[key])
			else:
				print("Error: Missing stat: " + key + ".")
				return None


	###################################
	# Pipeline termination functions
	###################################

	def stop_pipeline(self):
		"""
		Terminates the pipeline.

		The normal pipeline completion function, to be run by the pipeline
		at the end of the script. It sets status flag to completed and records 
		some time and memory statistics to the log file.
		"""
		self.set_status_flag("completed")
		self._cleanup()
		self.report_result("Time", str(datetime.timedelta(seconds = self.time_elapsed(self.starttime))))
		self.report_result("Success", time.strftime("%m-%d %H:%M:%S"))
		print("\n##### [Epilogue:]")
		print("* " + "Total elapsed time".rjust(20) + ":  " + str(datetime.timedelta(seconds = self.time_elapsed(self.starttime))))
		# print("Peak memory used: " + str(memory_usage()["peak"]) + "kb")
		print("* " + "Peak memory used".rjust(20) + ":  " + str(round(self.peak_memory, 2)) + " GB")
		self.timestamp("* Pipeline completed at: ".rjust(20))

	def fail_pipeline(self, e):
		"""
		If the pipeline does not complete, this function will stop the pipeline gracefully.
		It sets the status flag to failed and skips the	normal success completion procedure.

		:param e: Exception to raise.
		:type e: Exception
		"""
		# Take care of any active running subprocess
		if self.running_subprocess is not None:
			pid_to_kill = self.running_subprocess.pid
			self.running_subprocess = None
			# Close the preformat tag that we opened when the process was spawned.
			print("</pre>")
			# record profile of any running processes before killing
			elapsed_time = time.time() - self.proc_start_time
			self._report_profile(self.proc_name, self.proc_lock_name, elapsed_time, self.peak_memory)

			self._kill_child_process(pid_to_kill)
		if self.status != "failed":  # and self.status != "completed":
			self.set_status_flag("failed")
			self.timestamp("### Pipeline failed at: ")
			print("Total time: ", str(datetime.timedelta(seconds = self.time_elapsed(self.starttime))))

		raise e

	def _signal_term_handler(self, signal, frame):
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

	def _signal_int_handler(self, signal, frame):
		"""
		For catching interrupt (Ctrl +C) signals. Fails gracefully.
		"""
		message = "Got SIGINT (Ctrl +C); Failing gracefully..."
		with open(self.pipeline_log_file, "a") as myfile:
			myfile.write(message + "\n")
		self.fail_pipeline(Exception("SIGINT"))
		sys.exit(1)

	def _exit_handler(self):
		"""
		This function I register with atexit to run whenever the script is completing.
		A catch-all for uncaught exceptions, setting status flag file to failed.
		"""
		# print("Exit handler")
		# Make the cleanup file executable if it exists
		if os.path.isfile(self.cleanup_file):
			# Make the cleanup file self destruct.
			with open(self.cleanup_file, "a") as myfile:
				myfile.write("rm " + self.cleanup_file + "\n")
			os.chmod(self.cleanup_file, 0o755)

		# If the pipeline hasn't completed successfully, or already been marked
		# as failed, then mark it as failed now.

		if self.status != "completed" and self.status != "failed":
			self.fail_pipeline(Exception("Unknown exit failure"))

	def _kill_child_process(self, child_pid):
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
			print("\nPypiper terminating spawned child process " + str(child_pid))
			# sys.stdout.flush()
			os.kill(child_pid, signal.SIGTERM)

	def clean_add(self, regex, conditional=False, manual=False):
		"""
		Add files (or regexs) to a cleanup list, to delete when this pipeline completes successfully.
		When making a call with run that produces intermediate files that should be
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
						if os.path.isfile(file): myfile.write("rm " + file + "\n")
						elif os.path.isdir(file): myfile.write("rmdir " + file + "\n")
			except:
				pass
		elif conditional:
			self.cleanup_list_conditional.append(regex)
		else:
			self.cleanup_list.append(regex)
			# Remove it from the conditional list if added to the absolute list
			while regex in self.cleanup_list_conditional: self.cleanup_list_conditional.remove(regex)

	def _cleanup(self):
		"""
		Cleans up (removes) intermediate files.

		You can register intermediate files, which will be deleted automatically
		when the pipeline completes. This function deletes them,
		either absolutely or conditionally. It is run automatically when the
		pipeline succeeds, so you shouldn't need to call it from a pipeline.
		"""
		if len(self.cleanup_list) > 0:
			print("\nCleaning up flagged intermediate files...")
			for expr in self.cleanup_list:
				print("\nRemoving glob: " + expr)
				try:
					# Expand regular expression
					files = glob.glob(expr)
					# Remove entry from cleanup list
					while files in self.cleanup_list: self.cleanup_list.remove(files)
					# and delete the files
					for file in files:
						if os.path.isfile(file):
							print("`rm " + file + "`")
							os.remove(os.path.join(file))
						elif os.path.isdir(file):
							print("`rmdir " + file + "`")
							os.rmdir(os.path.join(file))
				except:
					pass

		if len(self.cleanup_list_conditional) > 0:
			# flag_files = glob.glob(self.pipeline_outfolder + "*.flag")
			flag_files = [fn for fn in glob.glob(self.pipeline_outfolder + "*.flag") if "completed" not in os.path.basename(fn) and not self.pipeline_name + "_running.flag" == os.path.basename(fn)]
			if (len(flag_files) == 0):
				print("\nCleaning up conditional list...")
				for expr in self.cleanup_list_conditional:
					print("\nRemoving glob: " + expr)
					try:
						files = glob.glob(expr)
						while files in self.cleanup_list_conditional: self.cleanup_list_conditional.remove(files)
						for file in files:
							if os.path.isfile(file):
								print("`rm " + file + "`")
								os.remove(os.path.join(file))
							elif os.path.isdir(file):
								print("`rmdir " + file + "`")
								os.rmdir(os.path.join(file))
					except:
						pass
			else:
				print("\nConditional flag found: " + str([os.path.basename(i) for i in flag_files]))
				print("\nThese conditional files were left in place:" + str(self.cleanup_list_conditional))
				# Produce a cleanup script
				for expr in self.cleanup_list_conditional:
					try:
						files = glob.glob(expr)
						for file in files:
							with open(self.cleanup_file, "a") as myfile:
								if os.path.isfile(file): myfile.write("rm " + file + "\n")
								elif os.path.isdir(file): myfile.write("rmdir " + file + "\n")
					except:
						pass

	def _memory_usage(self, pid='self', category="hwm"):
		"""
		Memory usage of the process in kilobytes.

		:param pid: Process ID of process to check
		:type pid: str
		:param category: Memory type to check. 'hwm' for high water mark.
		:type category: str
		"""
		# Thanks Martin Geisler:
		status = None
		result = {'peak': 0, 'rss': 0, 'hwm': 0}
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
		# print(result[category])
		return result[category]


class Tee(object):
	def __init__(self, log_file):
		self.file = open(log_file, "a")
		self.stdout = sys.stdout
		sys.stdout = self

	def __del__(self):
		sys.stdout = self.stdout
		self.file.close()

	def write(self, data):
		self.file.write(data)
		self.stdout.write(data)
		self.file.flush()
		self.stdout.flush()

	def fileno(self):
		return(self.stdout.fileno())


# @staticmethod
def add_pypiper_args(parser, groups = ["pypiper"], args = [None], all_args = False):
	"""
	Adds default automatic args to an ArgumentParser. Use this to add standardized 
	pypiper arguments to your python pipeline.

	There are two ways to use `add_pypiper_args`: by specifying argument groups,
	or by specifying individual arguments. Specifying argument groups will add
	multiple arguments to your parser; these convenient argumenet groupings make it
	easy to add arguments to certain types of pipeline. For example, to make a
	looper-compatible pipeline, use `groups = ["pypiper", "looper"]`.

	:param parser: an ArgumentParser object from your pipeline
	:param groups: Adds arguments belong to specified group of args.
		 Options are: pypiper, config, looper, resources, common, ngs, all.
	:type groups: list
	:param args: You may specify a list of specific arguments one by one.
	:type args: list

	:returns: A new ArgumentParser object, with selected pypiper arguments added
	"""

	# A Dict that defines groups of arguments. If the keys in this dict
	# are specified under "groups", then all of the arguments for any specified
	# group will be added.

	args_list = {"pypiper" : ["recover", "new-start", "dirty", "follow"],
				"config" : ["config"],
				"resource" : ["mem", "cores"],
				"looper" : ["config", "output-parent", "mem", "cores"],
				"common" : ["input", "sample-name"],
				"ngs" : ["input", "sample-name", "input2", "genome", "single-or-paired"]
				}



	# any argument in any group is added to a new group called "all"
	args_list["all"] = list(set(sum(args_list.values(), [])))

	args_to_add = set()
	# accumulate args to add based on parameters passed to this method

	#if basic:
	#	args_to_add.update(args_list["basic"])

	if groups:
		if type(groups) != list:
			groups = [groups]

		# For backwards-compatilibity
		if all_args:
			groups.append("all")

		for group in groups:
			args_to_add.update(args_list[group])

	if args:
		if type(args) != list:
			args = [args]
		args_to_add.update(args)

	#print(args_to_add)

	for arg in args_to_add:	
		# Basic pypiper arguments actually used by pypiper
		if arg == "recover":
			parser.add_argument(
				'-R', '--recover', dest='recover', action='store_true',
				default=False, help='Recover mode, overwrite locks')
		if arg == "new-start":
			parser.add_argument(
				'-N', '--new-start', dest='fresh', action='store_true',
				default=False, help='Fresh start mode, overwrite all')
		if arg == "dirty":
			parser.add_argument(
				'-D', '--dirty', dest='manual_clean', action='store_true',
				default=False, help='Make all cleanups manual')  # Useful for debugging
		if arg == "follow":
			parser.add_argument(
				'-F', '--follow', dest='force_follow', action='store_true',
				default=False, help='Run all follow commands, even if command is not run')  # Recalculating stats

		if arg == "config":
			default_config = os.path.splitext(os.path.basename(sys.argv[0]))[0] + ".yaml"
			# Arguments to optimize the interface to looper
			parser.add_argument(
				"-C", "--config", dest="config_file", type=str,
				help="pipeline config file in YAML format; relative paths are \
				considered relative to the pipeline script. \
				defaults to " + default_config,
				required=False, default=default_config, metavar="CONFIG_FILE")
		if arg == "output-parent":	
			parser.add_argument(
				"-O", "--output-parent", dest="output_parent", type=str,
				help="parent output directory of the project (required).",
				required=False, metavar="PARENT_OUTPUT_FOLDER")
			# output_parent was previously called project_root

		if arg == "cores":
			parser.add_argument(
				"-P", "--cores", dest="cores", type=str,
				help="number of cores to use for parallel processes",
				required=False, default=1, metavar="NUMBER_OF_CORES")
		if arg == "mem":	
			parser.add_argument(
				"-M", "--mem", dest="mem", type=str,
				help="Memory string for processes that accept memory limits (like java)",
				required=False, default="4000", metavar="MEMORY_LIMIT")

		if arg == "input":
			# Arguments typically used in every pipeline
			parser.add_argument(
				"-I", "--input", dest="input", type=str, nargs="+",
				help="One or more primary input files (required)",
				required=False, metavar="INPUT_FILES")
			# input was previously called unmapped_bam
		if arg == "sample-name":
			parser.add_argument(
				"-S", "--sample-name", dest="sample_name", type=str,
				help="unique name for output subfolder and files (required)",
				required=False, metavar="SAMPLE_NAME")

		if arg == "input2":
			# Common arguments specific to NGS pipelines
			parser.add_argument(
				"-I2", "--input2", dest="input2", type=str, nargs="+",
				help="One or more secondary input files (if they exists); \
				for example, second read in pair.",
				required=False, default=None, metavar="INPUT_FILES2")
		if arg == "genome":
			parser.add_argument(
				"-G", "--genome", dest="genome_assembly", type=str,
				help="identifier for genome assempbly (required)",
				required=False)
		if arg == "single-or-paired":
			parser.add_argument(
				"-Q", "--single-or-paired", dest="single_or_paired", type=str,
				help="single or paired end? default: single",
				required=False, default="single")

	return(parser)
