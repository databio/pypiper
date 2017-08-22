
Introduction
=========================

Pypiper is a lightweight python toolkit for gluing together restartable command
line pipelines. With Pypiper, **simplicity is paramount**. It should take less
than 15 minutes to build your first pipeline. Learning all the
:doc:`features and benefits <features>` takes just an hour or two. At
the same time, Pypiper provides immediate advantages over a
simple shell script.

Who should use Pypiper?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you have a shell script that would benefit from a layer of "handling code",
Pypiper helps you convert that set of shell commands into a production-scale
workflow, automatically handling the annoying details (restartablilty, file
integrity, logging) to make your pipeline robust and restartable.

If you need a full-blown, datacenter-scale environment that can do everything,
look elsewhere. Pypiper's strength is its simplicity. If all you want is a
shell-like script, but now with the power of python, and restartability, then
Pypiper is for you.

This emphasis on simplicity provides a few advantages:

- Write your pipeline in pure python (no new language to learn).
- Pypiper is easy to learn (3 or 4 functions will be all you need for simple
  stuff)
- Pypiper does not assume you want a complex dependency structure. You write a
  simple **ordered sequence of commands**, just like a shell script.


What Pypiper does NOT do:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Pypiper tries to exploit the `Pareto principle
<https://en.wikipedia.org/wiki/Pareto_principle>`_ -- you'll get 80% of the
features with only 20% of the work of other pipeline management systems. So,
there are a few things Pypiper deliberately doesn't do:


- Job dependencies. Pypiper runs sequential pipelines. If you want to implement
  a pipeline with complex task dependencies, there are better options. We view
  this as an advantage because it makes the pipeline easier to write, easier to
  understand, and easier to debug. For developmental pipelines, the complexity
  cost is not worth the minimal benefit -- read this `post on parallelism in bioinformatics <http://databio.org/posts/paralellism_in_bioinformatics.html>`_ 
  for an explanation.

- Cluster submission. Pypiper does not handle any sort of cluster job submission
  or  resource requesting. Instead, we have divided this into a
  separate project called `looper <http://looper.readthedocs.io/>`_. This makes
  a modular system: you can use whatever system you want for cluster management.
  `Pypiper <http://pypiper.readthedocs.io/>`_ builds individual, single-sample
  pipelines that can be run one sample at a time. `Looper
  <http://looper.readthedocs.io/>`_ then processes groups of samples, submitting
  appropriate pipelines to a cluster or server. The two projects are independent
  and can be used separately, but they are most powerful when combined.



Does the world need another pipeline management system?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As I began to put together production-scale pipelines, I found a
lot of relevant pipelining systems, but was universally disappointed. For my
needs, they were all overly complex. I wanted something simple enough to
**quickly write and maintain** a pipeline without having to learn a lot of new
functions and conventions, but robust enough to handle requirements like
restartability and memory usage monitoring. Everything related was either a pre-
packaged pipeline for a defined purpose, or a heavy-duty development environment
that was overkill for a simple pipeline. Both of these seemed to be targeted
toward ultra-efficient uses, and neither fit my needs: I had a set of commands
already in mind -- I just needed a wrapper that could take that code and make it
automatically restartable, logged, robust to crashing, easy to debug, and so
forth.

Pypiper fills a niche: the individual graduate student who is just putting
together a simple series of commands and wants to run that on a few hundred or
thousand samples. It's not worth learning a workflow development language and a
complex data-center-class tool for your simple research pipeline.