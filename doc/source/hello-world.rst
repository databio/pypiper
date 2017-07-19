Installing and Hello World
==============================

Release versions are posted on the GitHub `releases page <https://github.com/epigen/pypiper/releases>`_. You can install the latest version directly from GitHub using pip:

.. code-block:: bash

	pip install --user https://github.com/epigen/pypiper/zipball/master


Update with:

.. code-block:: bash

	pip install --user --upgrade https://github.com/epigen/pypiper/zipball/master

Now, to test pypiper, follow the commands in the ``Hello, Pypiper!`` tutorial: just run these 3 lines of code and you're running your first pypiper pipeline!

.. code:: bash

	# Install the latest version of pypiper:
	pip install --user https://github.com/epigen/pypiper/zipball/master

	# download hello_pypiper.py
	wget https://raw.githubusercontent.com/epigen/pypiper/master/example_pipelines/hello_pypiper.py
	
	# Run it:
	python hello_pypiper.py


You should see printed to screen some output like this:

.. code::

	----------------------------------------
	##### [Pipeline run code and environment:]
	*              Command:  `hello_pypiper.py`
	*         Compute host:  puma
	*          Working dir:  /home/nsheff/code/pypiper/example_pipelines
	*            Outfolder:  hello_pypiper_results/
	*  Pipeline started at:   (03-02 14:58:20) elapsed:0:00:00 _TIME_

	##### [Version log:]
	*       Python version:  2.7.12
	*          Pypiper dir:  `/home/nsheff/.local/lib/python2.7/site-packages/pypiper`
	*         Pipeline dir:  `/home/nsheff/code/pypiper/example_pipelines`
	*     Pipeline version:  fa2c12d8fba0a9b50a678e7682f39a7a756e0ead
	*      Pipeline branch:  * dev
	*        Pipeline date:  2017-01-27 14:05:15 -0500

	##### [Arguments passed to pipeline:]

	----------------------------------------


	Change status from initializing to running
	Hello! (03-02 14:58:20) elapsed:0:00:00 _TIME_

	Target to produce: `hello_pypiper_results/output.txt`
	> `echo 'Hello, Pypiper!' > hello_pypiper_results/output.txt`

	<pre>
	</pre>
	Process 22693 returned: (0). Elapsed: 0:00:00.

	Change status from running to completed
	> `Time`	0:00:00	hello_pypiper	_RES_
	> `Success`	03-02 14:58:20	hello_pypiper	_RES_

	##### [Epilogue:]
	*   Total elapsed time:  0:00:00
	*     Peak memory used:  0.0 GB
	* Pipeline completed at:  (03-02 14:58:20) elapsed:0:00:00 _TIME_

	Pypiper terminating spawned child process 22679

Now observe your results in the folder ``hello_pypiper_results``:

 * hello_pypiper_commands.sh
 * hello_pypiper_completed.flag
 * hello_pypiper_log.md
 * hello_pypiper_profile.tsv
 * output.txt
 * stats.tsv

These files are explained in more detail in the :doc:`Outputs <outputs>` section.
