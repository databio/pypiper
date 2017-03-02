
FAQ
=========================

-   **How can I run my pipeline on more than 1 sample?**
	Pypiper only handles individual-sample pipelines. To run it on multiple samples, write a loop, or use `Looper <http://looper.readthedocs.io/>`_. Dividing multi-sample handling from individual sample handling is a conceptual advantage that allows us to write a nice, universal, generic sample-handler that you only have to learn once.

-   **What cluster resources can pypiper use?** 
	PyPiper is compute-agnostic. You run it wherever you want; If you want a nice way to submit pipelines for samples any cluster manager, check out `Looper <http://looper.readthedocs.io/>`_.
