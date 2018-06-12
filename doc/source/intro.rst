
Introduction
=========================

Pypiper is a lightweight python toolkit for gluing together restartable command
line pipelines. With Pypiper, **simplicity is paramount**. It should take less
than 15 minutes to build your first pipeline. Learning all the
:doc:`features and benefits <features>` takes just an hour or two. At
the same time, Pypiper provides immediate advantages over a
simple shell script.

Pypiper is an example of a simple  `bioinformatics pipeline framework
<http://databio.org/pipeline_frameworks/>`_. It differs from existing frameworks in its focus on **simplicity** and **sequential pipelines**. 
To employ pypiper, you will just take your bash script and pass those commands through ``PipelineManager.run()``. This will give you automatic restartability, process monitoring for memory use and compute time, pipeline status monitoring, copious log output, robust error handling, easy debugging tools, guaranteed file output integrity, and a bunch of useful pipeline development helper functions.

A simple example pipeline
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To demonstrate the simplicity, take a look at a very simple but complete pipeline:

.. literalinclude:: ../../example_pipelines/hello_pypiper.py

There's nothing complex here: we are choosing an output folder, and running a single command, much like you may do in a shell script. Building pypiper pipelines is as simple as stringing together shell commands. That's it. We'll actually run this example pipeline in the ``Hello World`` section.


Who should use Pypiper?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The target audience for pypiper is an individual student or researcher or lab
that wants to build a basic pipeline, but to do it better job than just writing
a shell script. Many bioinformatics pipelines are written by students or
technicians who don't have time to learn a full-scale pipelining framework, so
they just end up using simple bash scripts to piece together commands because
that seems the most accessible. Pypiper tries to give 80% of the benefits of a
professional-scale pipelining system while requiring very little additional
effort.

If you have a shell script that would benefit from a layer of "handling code",
Pypiper helps you convert that set of shell commands into a production-scale
workflow, automatically handling the annoying details (restartablilty, file
integrity, logging) to make your pipeline robust and restartable.

If you need a full-blown, datacenter-scale environment that can do everything,
look elsewhere. Pypiper's strength is its simplicity. If all you want is a
shell-like script, but now with the power of python, some built-in benefits, and
syntactic sugar, then Pypiper is for you.

This emphasis on simplicity provides a few advantages:

- Write your pipeline in pure python (no new language to learn).
- Pypiper is easy to learn (3 or 4 functions will be all you need for simple
  stuff)
- Pypiper does not assume you want a complex dependency structure. You write a
  simple **ordered sequence of commands**, just like a shell script.


What Pypiper does NOT do
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Pypiper tries to exploit the `Pareto principle
<https://en.wikipedia.org/wiki/Pareto_principle>`_ -- you'll get 80% of the
features with only 20% of the work of other pipeline management systems. So,
there are a few things Pypiper deliberately doesn't do:


- Task dependencies. Pypiper runs sequential pipelines. If you want to implement
  a pipeline with complex task dependencies, there are better options. We view
  this as an advantage because it makes the pipeline easier to write, easier to
  understand, and easier to debug -- critical things for pipelines that are
  still under active development (which is, really, *all* pipelines). For
  developmental pipelines, the complexity cost of encoding task dependencies is
  not worth the minimal benefit -- read this `post on parallelism in
  bioinformatics <http://databio.org/posts/paralellism_in_bioinformatics.html>`_
  for an explanation.

- Cluster submission. Pypiper does not handle any sort of cluster job submission
  or resource requesting. Instead, we have divided this into a separate project
  called `looper <http://looper.readthedocs.io/>`_. This makes a modular system:
  you can use whatever system you want for cluster management. `Pypiper
  <http://pypiper.readthedocs.io/>`_ builds individual, single-sample pipelines
  that can be run one sample at a time. `Looper
  <http://looper.readthedocs.io/>`_ then processes groups of samples, submitting
  appropriate pipelines to a cluster or server. The two projects are independent
  and can be used separately, but they are most powerful when combined. This
  keeps things simple and modular.


Yet another pipeline system?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As I began to put together production-scale pipelines, I found a lot of relevant
pipelining systems, but was universally disappointed. For my needs, they were
all overly complex. I wanted something **simple enough to quickly write and
maintain** a pipeline without having to learn a lot of new functions and
conventions, but robust enough to handle requirements like restartability and
memory usage monitoring. Everything related was either a pre-packaged pipeline
for a defined purpose, or a heavy-duty development environment that was overkill
for a simple pipeline. Both of these seemed to be targeted toward ultra-
efficient uses, and neither fit my needs: I had a set of commands already in
mind -- I just needed a wrapper that could take that code and make it
automatically restartable, logged, robust to crashing, easy to debug, and so
forth.
