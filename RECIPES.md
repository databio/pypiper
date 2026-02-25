# Pypiper Cookbook

Copy-paste recipes for common pypiper patterns. Each recipe is a complete, runnable script. For API details, see the full documentation. For bioinformatics-specific examples, see `example_pipelines/`.

## Table of Contents

1. [Run a single command with restart support](#recipe-1-run-a-single-command-with-restart-support)
2. [Run a multi-step pipeline with checkpoints](#recipe-2-run-a-multi-step-pipeline-with-checkpoints)
3. [Capture command output into a Python variable](#recipe-3-capture-command-output-into-a-python-variable)
4. [Multiple PipelineManagers in one script](#recipe-4-multiple-pipelinemanagers-in-one-script)
5. [Report results to a structured output file](#recipe-5-report-results-to-a-structured-output-file)
6. [Skip expensive steps when output already exists](#recipe-6-skip-expensive-steps-when-output-already-exists)
7. [Clean up intermediate files](#recipe-7-clean-up-intermediate-files)
8. [Handle a command that might fail (nofail)](#recipe-8-handle-a-command-that-might-fail-nofail)
9. [Use follow functions for post-processing](#recipe-9-use-follow-functions-for-post-processing)
10. [Parallelism patterns](#recipe-10-parallelism-patterns)
11. [Interactive / notebook usage](#recipe-11-interactive--notebook-usage)

---

## Recipe 1: Run a single command with restart support

Run a shell command, skip it on re-run if the output file already exists, and handle stale lock files from crashed runs.

```python
#!/usr/bin/env python
"""Recipe 1: Run a single command with automatic restart support."""
import pypiper

pm = pypiper.PipelineManager(name="file_processor", outfolder="output/")
pm.run("wc -l /etc/passwd > output/line_count.txt", target="output/line_count.txt")
pm.stop_pipeline()
# Run again: the command is skipped because output/line_count.txt exists.
# Run with --recover: ignores stale lock files from a crashed run.
# Run with --new-start: forces re-run even if target exists.
```

**Alternative: context manager style** (calls `stop_pipeline()` automatically on exit):

```python
#!/usr/bin/env python
"""Recipe 1 (context manager variant)."""
import pypiper

with pypiper.PipelineManager(name="file_processor", outfolder="output/") as pm:
    pm.run("wc -l /etc/passwd > output/line_count.txt", target="output/line_count.txt")
# stop_pipeline() is called automatically when the `with` block exits.
```

**What happens:** The command runs and creates `output/line_count.txt`. On re-run, the command is skipped because the target file already exists. Pass `--new-start` to force re-run.

---

## Recipe 2: Run a multi-step pipeline with checkpoints

Use checkpoints to enable restarting a pipeline from any named step.

```python
#!/usr/bin/env python
"""Recipe 2: Multi-step pipeline with checkpoints for restart."""
import pypiper

pm = pypiper.PipelineManager(name="sorter", outfolder="output/")

pm.timestamp("### Step 1: Generate data", checkpoint="generate")
pm.run("shuf -i 1-1000 -n 100 > output/numbers.txt", target="output/numbers.txt")

pm.timestamp("### Step 2: Sort data", checkpoint="sort")
pm.run("sort -n output/numbers.txt > output/sorted.txt", target="output/sorted.txt")

pm.timestamp("### Step 3: Extract top 10", checkpoint="top10")
pm.run("head -n 10 output/sorted.txt > output/top10.txt", target="output/top10.txt")

pm.stop_pipeline()
# Run with: python recipe2.py --start-point sort
# This skips Step 1 and starts from Step 2.
```

**What happens:** Three steps run in sequence, each with a named checkpoint. Passing `--start-point sort` on re-run skips Step 1 and begins from Step 2.

---

## Recipe 3: Capture command output into a Python variable

Use `checkprint()` to capture stdout from a shell command and `report_result()` to log key-value metrics.

```python
#!/usr/bin/env python
"""Recipe 3: Capture command output and report structured results."""
import pypiper

pm = pypiper.PipelineManager(name="system_info", outfolder="output/")

hostname = pm.checkprint("hostname")
pm.report_result("hostname", hostname)

timestamp = pm.checkprint("date +%s")
pm.report_result("unix_timestamp", timestamp)

line_count = pm.checkprint("wc -l < /etc/passwd")
pm.report_result("passwd_lines", int(line_count))

pm.stop_pipeline()
# Results are written to output/stats.yaml via pipestat.
```

**What happens:** Each `checkprint()` runs a command and returns its stdout as a string. `report_result()` persists the key-value pairs to the stats file in the output folder.

---

## Recipe 4: Multiple PipelineManagers in one script

Use `multi=True` to run multiple independent PipelineManagers without log file conflicts.

```python
#!/usr/bin/env python
"""Recipe 4: Multiple PipelineManagers in one script."""
import pypiper

# multi=True disables output tee-ing so managers don't conflict on log files.
pm_a = pypiper.PipelineManager(name="task_a", outfolder="output/task_a/", multi=True)
pm_b = pypiper.PipelineManager(name="task_b", outfolder="output/task_b/", multi=True)

pm_a.run("echo 'A done' > output/task_a/result.txt", target="output/task_a/result.txt")
pm_b.run("echo 'B done' > output/task_b/result.txt", target="output/task_b/result.txt")

pm_a.stop_pipeline()
pm_b.stop_pipeline()
# Each manager has its own output folder, locks, and status tracking.
# multi=True is essential when running multiple PMs in one script.
# Commands here run sequentially; for parallelism patterns see Recipe 10.
```

**What happens:** Two independent pipeline managers run in the same script, each with its own output folder and state. The `multi=True` flag prevents log file conflicts.

---

## Recipe 5: Report results to a structured output file

Use `report_result()` to persist metrics extracted from command output.

```python
#!/usr/bin/env python
"""Recipe 5: Report structured results to a YAML output file."""
import pypiper

pm = pypiper.PipelineManager(name="analyzer", outfolder="output/")

# Generate a file and measure it
pm.run("seq 1 10000 > output/data.txt", target="output/data.txt")

line_count = pm.checkprint("wc -l < output/data.txt")
pm.report_result("line_count", int(line_count))

file_size = pm.checkprint("stat --format=%s output/data.txt")
pm.report_result("file_size_bytes", int(file_size))

word_count = pm.checkprint("wc -w < output/data.txt")
pm.report_result("word_count", int(word_count))

pm.stop_pipeline()
# Results are persisted via pipestat. Check output/ for the results file.
# report_result() also logs each result to the pipeline log.
```

**What happens:** The pipeline generates data, extracts metrics with `checkprint()`, and stores them via `report_result()`. Results are persisted to the pipestat-backed output file (YAML by default).

---

## Recipe 6: Skip expensive steps when output already exists

The `target` parameter is the key mechanism -- pypiper checks for file existence before running.

```python
#!/usr/bin/env python
"""Recipe 6: Skip expensive steps when output files already exist."""
import pypiper

pm = pypiper.PipelineManager(name="expensive", outfolder="output/")

# This simulates an expensive computation (3 seconds).
# On re-run, it's skipped instantly because the target file exists.
pm.run(
    "sleep 3 && seq 1 1000000 > output/big_data.txt",
    target="output/big_data.txt",
)

# This step also skips if its target exists.
pm.run(
    "wc -l output/big_data.txt > output/summary.txt",
    target="output/summary.txt",
)

pm.stop_pipeline()
# First run: takes ~3 seconds. Second run: completes instantly.
# Use --new-start to force re-computation of all steps.
```

**What happens:** The first run takes about 3 seconds. The second run completes instantly because both target files already exist. Pass `--new-start` to force re-computation.

---

## Recipe 7: Clean up intermediate files

Use `clean=True` or `clean_add()` to register intermediate files for automatic deletion on success.

```python
#!/usr/bin/env python
"""Recipe 7: Clean up intermediate files after pipeline success."""
import pypiper

pm = pypiper.PipelineManager(name="cleaner", outfolder="output/")

# Method 1: clean=True in run() auto-registers target for cleanup
pm.run(
    "seq 1 100 > output/intermediate.txt",
    target="output/intermediate.txt",
    clean=True,  # This file will be deleted on success
)

# Method 2: Explicit clean_add() for manual registration
pm.run("cp output/intermediate.txt output/temp_copy.txt", target="output/temp_copy.txt")
pm.clean_add("output/temp_copy.txt")

# Final output (not cleaned up)
pm.run(
    "cat output/intermediate.txt output/temp_copy.txt > output/final.txt",
    target="output/final.txt",
)

pm.stop_pipeline()
# On success: intermediate.txt and temp_copy.txt are deleted; final.txt remains.
# Run with --dirty to keep all intermediate files (useful for debugging).
```

**What happens:** On successful completion, intermediate files are deleted automatically. Only `output/final.txt` remains. Pass `--dirty` to preserve all intermediate files for debugging.

---

## Recipe 8: Handle a command that might fail (nofail)

Use `nofail=True` to let the pipeline continue past a command that returns a non-zero exit code.

```python
#!/usr/bin/env python
"""Recipe 8: Continue past a command that might fail."""
import pypiper

pm = pypiper.PipelineManager(name="resilient", outfolder="output/")

# This command will fail (nonexistent path), but nofail=True lets the pipeline continue.
pm.run("ls /nonexistent_path_12345 > output/maybe.txt", target="output/maybe.txt", nofail=True)

# This command runs regardless of the previous failure.
pm.run("echo 'Pipeline continued' > output/success.txt", target="output/success.txt")

pm.stop_pipeline()
# Without nofail=True, the failed command would halt the entire pipeline.
# With nofail=True, the failure is logged and execution continues.
```

**What happens:** The first command fails but the pipeline continues because of `nofail=True`. Without it, the pipeline would halt on the first failure.

---

## Recipe 9: Use follow functions for post-processing

The `follow` parameter runs a Python function after the command, only when the command actually executes.

```python
#!/usr/bin/env python
"""Recipe 9: Post-process command output with follow functions."""
import pypiper

pm = pypiper.PipelineManager(name="follower", outfolder="output/")

def report_line_count():
    """This runs after the command, only when the command actually executes."""
    count = pm.checkprint("wc -l < output/data.txt")
    pm.report_result("line_count", int(count))

pm.run(
    "seq 1 500 > output/data.txt",
    target="output/data.txt",
    follow=report_line_count,
)

pm.stop_pipeline()
# The follow function runs only when the command runs (target doesn't exist).
# On re-run, both the command and follow are skipped.
# Use --force-follow to run follow functions even when commands are skipped.
```

**What happens:** The follow function runs after the command completes. On re-run, both the command and follow function are skipped because the target exists. Pass `--force-follow` to run follow functions even when commands are skipped.

---

## Recipe 10: Parallelism patterns

Pypiper runs commands sequentially by design. Here are the two ways to get parallelism.

```python
#!/usr/bin/env python
"""Recipe 10: Parallelism patterns with pypiper.

Pypiper runs commands sequentially. Passing a list to pm.run() does NOT
parallelize -- it runs commands one after another under a single lock.

To get parallelism, use one of these approaches:
"""
import pypiper

# APPROACH 1: Parallel pipelines via separate PipelineManager instances.
# Each PM has its own output folder, locks, and status tracking.
# Use multi=True to prevent log file conflicts.
# For large-scale parallelism, use looper to submit each as a separate job.
pm_a = pypiper.PipelineManager(name="task_a", outfolder="output/task_a/", multi=True)
pm_b = pypiper.PipelineManager(name="task_b", outfolder="output/task_b/", multi=True)

pm_a.run("echo 'A done' > output/task_a/result.txt", target="output/task_a/result.txt")
pm_b.run("echo 'B done' > output/task_b/result.txt", target="output/task_b/result.txt")

pm_a.stop_pipeline()
pm_b.stop_pipeline()

# APPROACH 2: Tool-level parallelism -- let the tool use multiple threads.
# Pypiper manages the lifecycle; the tool manages internal parallelism.
# pm.run("samtools sort -@ 8 input.bam", target="sorted.bam")
# pm.run("pigz -p 4 big_file.txt", target="big_file.txt.gz")

# ANTI-PATTERN: Don't call pm.run() from multiple threads on the same PM.
# The shared state (locks, running_procs) is not thread-safe, and the
# check-target-then-run logic would race. Use separate PMs instead.
```

**What happens:** Each PipelineManager runs independently with its own locks and output folder. For true parallelism across samples, use looper to submit pipeline scripts as separate processes/jobs.

---

## Recipe 11: Interactive / notebook usage

Use `multi=True` in Jupyter notebooks or interactive sessions to avoid log tee-ing conflicts.

```python
"""Recipe 11: Interactive / notebook usage."""
import pypiper

# multi=True is essential for notebooks -- it disables log tee-ing
# which would interfere with notebook output.
pm = pypiper.PipelineManager(name="notebook", outfolder="output/", multi=True)

# Run commands interactively
pm.run("echo 'hello' > output/greeting.txt", target="output/greeting.txt")

# Capture output directly into Python variables
result = pm.checkprint("cat output/greeting.txt")
print(f"Got: {result}")

# Report results as you go
pm.report_result("greeting", result)

# When done exploring, finalize:
pm.stop_pipeline()
# In a notebook, you can also just let the PM go out of scope.
# multi=True is the key -- without it, pypiper tries to tee stdout
# to a log file, which conflicts with notebook output handling.
```

**What happens:** The pipeline manager works in interactive mode without interfering with notebook output. `multi=True` is the key setting that makes this work.

---

## Quick Reference

| Method | Purpose |
|--------|---------|
| `PipelineManager(name, outfolder)` | Create a pipeline manager |
| `pm.run(cmd, target=)` | Run a shell command, skip if target exists |
| `pm.checkprint(cmd)` | Run a command, return its stdout as a string |
| `pm.report_result(key, value)` | Log a key-value result to the stats file |
| `pm.timestamp(msg, checkpoint=)` | Log a timestamp with optional checkpoint name |
| `pm.clean_add(path)` | Register a file for cleanup on success |
| `pm.stop_pipeline()` | Finalize the pipeline (write stats, clean up) |

| CLI Flag | Purpose |
|----------|---------|
| `--recover` | Ignore stale lock files from crashed runs |
| `--new-start` | Force re-run even if target files exist |
| `--start-point NAME` | Restart from a named checkpoint |
| `--dirty` | Keep intermediate files (skip cleanup) |
| `--force-follow` | Run follow functions even when commands are skipped |
