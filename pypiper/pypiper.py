#!/usr/env python

"""
Pypiper documentation
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


# Define global variables
# These just record some details about the pipeline

PIPELINE_NAME = ""
PEAKMEM = 0         # memory high water mark
STARTTIME = time()
LAST_TIMESTAMP = STARTTIME  # time of the last call to timestamp()
STATUS = "initializing"
PIPELINE_OUTFOLDER = ""
LOGFILE = ""  # required for termination signal code only.


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
    paths.log_file = paths.pipeline_outfolder + pipeline_name + "_log.md"
    global LOGFILE
    LOGFILE = paths.log_file
    make_sure_path_exists(paths.pipeline_outfolder)
    global PIPELINE_OUTFOLDER
    PIPELINE_OUTFOLDER = paths.pipeline_outfolder

    global STARTTIME
    STARTTIME = time()

    # Register handler functions to deal with interrupt and termination signals;
    # If recieved, we would then clean up properly (set pipeline status to FAIL, etc).
    signal.signal(signal.SIGINT, signal_int_handler)
    signal.signal(signal.SIGTERM, signal_term_handler)

    # Mirror every operation on sys.stdout to log file
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)  # Unbuffer output
    # a for append to file
    tee = subprocess.Popen(["tee", "-a", paths.log_file], stdin=subprocess.PIPE)
    # In case this pipeline process is terminated with SIGTERM, make sure we kill this spawned process as well.
    atexit.register(kill_child, tee.pid)
    os.dup2(tee.stdin.fileno(), sys.stdout.fileno())
    os.dup2(tee.stdin.fileno(), sys.stderr.fileno())

    # Record the git version of the pipeline being run. This code gets:
    # hash: the commit id of the last commit in this repo
    # date: the date of the last commit in this repo
    # diff: a summary of any differences in the current (run) version vs. the committed version
    git_commit_hash = subprocess.check_output("cd " + os.path.dirname(os.path.realpath(__file__)) + "; git rev-parse --verify HEAD", shell=True)
    git_commit_date = subprocess.check_output("cd " + os.path.dirname(os.path.realpath(__file__)) + "; git show -s --format=%ai HEAD", shell=True)
    git_commit_diff = subprocess.check_output("cd " + os.path.dirname(os.path.realpath(__file__)) + "; git diff --shortstat HEAD", shell=True)
    if (git_commit_diff == ""):
        git_commit_diff = "No uncommitted changes."
    start_time = time()  # this variable is never used, consider removing
    print("################################################################################")
    timestamp("Pipeline started at: ")
    print("Compute host:\t\t" + platform.node())
    print("Working dir : %s" % os.getcwd())
    print("Git pipeline version:\t\t" + git_commit_hash.strip())
    print("Git pipeline date:\t\t" + git_commit_date.strip())
    print("Git pipeline diff: \t\t" + git_commit_diff.strip())
    print("Python version:\t\t" + platform.python_version())
    print("Cmd: " + str(" ".join(sys.argv)))
    # Print all arguments
    argsDict = vars(args)
    for arg in argsDict:
        print(arg + ":\t\t" + str(argsDict[arg]))
    print("Run outfolder:\t\t" + paths.pipeline_outfolder)
    print("################################################################################")

    set_status_flag("running")

    return paths


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


def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def wait_for_lock(lock_file):
    sleeptime = 5
    first_message_flag = False
    dot_count = 0
    while os.path.isfile(lock_file):
        if first_message_flag is False:
            timestamp("Waiting for file lock: " + lock_file)
            first_message_flag = True
        else:
            sys.stdout.write(".")
            dot_count = dot_count + 1
            if dot_count % 60 == 0:
                print ""  # linefeed
        sleep(sleeptime)
        sleeptime = min(sleeptime + 5, 60)

    if first_message_flag:
        timestamp("File unlocked.")


def create_file(file):
    with open(file, 'w') as fout:
        fout.write('')

# This version will only succeed if this process actually
# creates the file; if the file already exists, it will
# raise an OSError


def create_file_racefree(file):
    write_lock_flags = os.O_CREAT | os.O_EXCL | os.O_WRONLY
    os.open(file, write_lock_flags)


# @pass_failure Should the pipeline bail on a nonzero return from a process? Default:  True
#               This can be used to implement non-essential parts of the pipeline.


def call_lock_nofail(cmd, lock_name, folder, output_file=None, shell=False):
    call_lock_internal(cmd, lock_name, folder, output_file, shell, pass_failure=False)


def call_lock(cmd, lock_name, folder, output_file=None, shell=False):
    call_lock_internal(cmd, lock_name, folder, output_file, shell, pass_failure=True)


def call_lock_internal(cmd, lock_name, folder, output_file=None, shell=False, pass_failure=True):
    # Create lock file:
    lock_file = os.path.join(folder, lock_name)
    ret = 0
    local_maxmem = 0

    # The loop here is unlikely to be triggered, and is just a wrapper
    # to prevent race conditions; the lock_file must be created by
    # the current loop. If not, we wait again for it and then
    # re-do the tests.
    while True:
        wait_for_lock(lock_file)
        if output_file is None or not (os.path.exists(output_file)):
            try:
                create_file_racefree(lock_file)     # Create lock
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
                        list_ret, list_maxmem = callprint(cmd_i, shell)
                        local_maxmem = max(local_maxmem, list_maxmem)
                        ret = max(ret, list_ret)
                else:  # Single command (most common)
                    ret, local_maxmem = callprint(cmd, shell)   # Run command
            except Exception as e:
                if (pass_failure):
                    fail_pipeline(e)
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


def callprint(cmd, shell=False):
    """
    split the command to use shell=False;
    leave it together to use shell=True; I use False so I can get the PID
    and poll memory use.
    """
    print "new callprint"
    global PEAKMEM
    global RUNNING_SUBPROCESS
    print(cmd)
    if not shell:
        if ("|" in cmd or ">" in cmd):
            print("Should this command run in a shell instead of directly in a subprocess?")
        cmd = cmd.split()
    # call(cmd, shell=shell) # old way (no memory profiling)

    p = subprocess.Popen(cmd, shell=shell)
    RUNNING_SUBPROCESS = p.pid
    local_maxmem = -1
    sleeptime = 5
    while p.poll() == None:
        if not shell:
            local_maxmem = max(local_maxmem, memory_usage(p.pid))
            # print ("int.maxmem (pid:" + str(p.pid) + ") " + str(local_maxmem))
        sleep(sleeptime)
        sleeptime = min(sleeptime + 5, 60)

    # set global maxmem
    PEAKMEM = max(PEAKMEM, local_maxmem)
    RUNNING_SUBPROCESS = None

    info = "Process " + str(p.pid) + " returned: (" + str(p.returncode) + ")."
    info += " Peak memory: (Process: " + str(local_maxmem) + "b;"
    info += " Pipeline: " + str(PEAKMEM) + "b)"
    print (info)
    if p.returncode != 0:
        raise Exception("Process returned nonzero result.")
    return [p.returncode, local_maxmem]

####################################
# Pipeline termination functions
###################################


def signal_term_handler(signal, frame):
    """
    Attached TERM signal handler function
    This may be invoked, for example, by SLURM if the job exceeds its memory or time limits.
    """
    message = "Got SIGTERM; Failing gracefully..."
    global LOGFILE
    with open(LOGFILE, "a") as myfile:
        myfile.write(message + "\n")
    fail_pipeline(Exception("SIGTERM"))
    sys.exit(1)


def signal_int_handler(signal, frame):
    """
    For catching interrupt (Ctrl +C) signals
    """
    message = "Got SIGINT; Failing gracefully..."
    global LOGFILE
    with open(LOGFILE, "a") as myfile:
        myfile.write(message + "\n")
    fail_pipeline(Exception("SIGINT"))
    sys.exit(1)


def exit_handler():
    '''
    This function I register with atexit to run whenever the script is completing
    '''
    # Catch-all for uncaught exceptions...
    global STATUS
    print("Exit handler")
    if RUNNING_SUBPROCESS is not None:
            kill_child(RUNNING_SUBPROCESS)
    if STATUS != "completed":
        set_status_flag("failed")


def kill_child(child_pid):
    if child_pid is None:
        pass
    else:
        print("Killing child process " + str(child_pid))
        os.kill(child_pid, signal.SIGTERM)


def fail_pipeline(e):
    set_status_flag("failed")
    raise e


def time_elapsed(time_since):
    return round(time() - time_since, 2)


def stop_pipeline(paths):
    """Remove temporary marker files to complete the pipeline"""
    global PEAKMEM
    global STARTTIME
    set_status_flag("completed")
    print ("Total elapsed time: " + str(time_elapsed(STARTTIME)))
    # print ("Peak memory used: " + str(memory_usage()["peak"]) + "kb")
    print ("Peak memory used: " + str(PEAKMEM / 1e6) + " GB")
    report_result("Success", strftime("%m-%d %H:%M:%S"), paths)
    timestamp("### Pipeline completed at: ")


def memory_usage(pid='self', category="peak"):
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
