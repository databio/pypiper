
Basic functions
=========================

Pypiper is simple, but powerful. You really only need to know about 3 functions to get started. ``PipelineManager`` can do:

.. currentmodule:: pypiper.manager.PipelineManager
.. autosummary:: 
	start_pipeline
	run
	stop_pipeline


With those 3 functions, you can create a simple pipeline. Click on each function to view its in-depth documentation. There are quite a few optional parameters to the ``run`` function, which is where most of Pypiper's power comes from.

When you've mastered the basics and are ready to get more powerful, add in a few new (optional) commands that make debugging and development easier:

.. autosummary:: 
	timestamp
	report_result
	clean_add
	get_stat


The complete documentation for these functions can be found in the :doc:`API <api>`.
