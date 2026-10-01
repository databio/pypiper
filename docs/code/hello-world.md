# Hello world

This brief tutorial will run your first basic pypiper pipeline to ensure you have everything set up correctly. 

Just run these 3 lines of code and you're running your first pypiper pipeline!

### Install the latest version of pypiper

```{console}
pip install --user piper
```


### Download hello_pypiper.py
```{console}
wget https://raw.githubusercontent.com/databio/pypiper/master/example_pipelines/hello_pypiper.py
```

This is a basic pipeline. Here are the contents:


```bash
cat ../example_pipelines/hello_pypiper.py
```

    #!/usr/bin/env python
    
    import pypiper
    outfolder = "hello_pypiper_results" # Choose a folder for your results
    
    # Create a PipelineManager, the workhorse of pypiper
    pm = pypiper.PipelineManager(name="hello_pypiper", outfolder=outfolder)
    
    # Timestamps to delineate pipeline sections are easy:
    pm.timestamp("Hello!")
    
    # Now build a command-line command however you like, and pass it to pm.run()
    target_file = "hello_pypiper_results/output.txt"
    cmd = "echo 'Hello, Pypiper!' > " + target_file
    pm.run(cmd, target_file)
    
    pm.stop_pipeline()


### Run it!


```bash
python3 ../example_pipelines/hello_pypiper.py
```

    ### [Pipeline run code and environment:]
    
    *              Command:  `../example_pipelines/hello_pypiper.py`
    *         Compute host:  nox
    *          Working dir:  /home/sheffien/code/pypiper/docs_jupyter
    *            Outfolder:  hello_pypiper_results/
    *  Pipeline started at:   (03-16 23:47:37) elapsed: 0.0 _TIME_
    
    ### [Version log:]
    
    *       Python version:  3.6.7
    *          Pypiper dir:  `/home/sheffien/.local/lib/python3.6/site-packages/pypiper`
    *      Pypiper version:  0.9.5dev
    *         Pipeline dir:  `/home/sheffien/code/pypiper/example_pipelines`
    *     Pipeline version:  None
    *        Pipeline hash:  b'134e8c8f723da66697ab4f5b204315979b4e1042\n'
    *      Pipeline branch:  b'* dev\n'
    *        Pipeline date:  b'2019-03-16 11:41:56 -0400\n'
    *        Pipeline diff:  b' 1 file changed, 16 insertions(+), 1 deletion(-)\n'
    
    ### [Arguments passed to pipeline:]
    
    
    ----------------------------------------
    
    
    Changed status from initializing to running.
    No config file
    Hello! (03-16 23:47:37) elapsed: 0.0 _TIME_
    
    Target to produce: `hello_pypiper_results/output.txt`
    
    
    > `echo 'Hello, Pypiper!' > hello_pypiper_results/output.txt`
    
    <pre>
    </pre>
    Process 128 returned: (0). Elapsed: 0:00:00. Peak memory: (Process: None; Pipeline: 0GB)
    
    Changed status from running to completed.
    
    > `Time`	0:00:00	hello_pypiper	_RES_
    
    > `Success`	03-16-23:47:37	hello_pypiper	_RES_
    
    ##### [Epilogue:]
    *   Total elapsed time:  0:00:00
    *     Peak memory used:  0 GB
    * Pipeline completed at:  (03-16 23:47:37) elapsed: 0.0 _TIME_
    
    Pypiper terminating spawned child process 114...(tee)


This output is printed to your screen and also recorded in a log file (called ``hello_pypiper_log.md``). There are a few other outputs from the pipeline as well. All results are placed in a folder called ``hello_pypiper_results``. Navigate to that folder to observe the output of the pipeline, which will include these files:

 * hello_pypiper_commands.sh
 * hello_pypiper_completed.flag
 * hello_pypiper_log.md
 * hello_pypiper_profile.tsv
 * output.txt
 * stats.tsv

These files are explained in more detail in the reference section [outputs explained](../outputs.md). 

What's next? That depends on if you're interested in just *running* pypiper pipelines, or if you want to *develop* pypiper pipelines. The next sections are a series of HOW-TO articles that address each of these scenarios.

