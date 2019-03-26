# Pypiper's development philosophy

## Who should use Pypiper?

The target audience for pypiper is an individual who wants to build a basic
pipeline, but **wants to do better job than just writing a shell script, without
requiring hours of energy to learn a new language or system**. Many
bioinformatics pipelines are written by students or technicians who don't have
time to learn a full-scale pipelining framework, so
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

## What Pypiper does NOT do

Pypiper tries to exploit the [Pareto principle](https://en.wikipedia.org/wiki/Pareto_principle) -- you'll get 80% of the
features with only 20% of the work of other pipeline management systems. So,
there are a few things Pypiper deliberately doesn't do:


- Task dependencies. Pypiper runs sequential pipelines. We view this as an
  advantage because it makes the pipeline easier to write, easier to understand,
  and easier to debug -- critical things for pipelines that are still under
  active development (which is, really, *all* pipelines). For developmental
  pipelines, the complexity introduced by task dependencies is not worth the
  minimal benefit -- read this [post on parallelism in
  bioinformatics](http://databio.org/posts/paralellism_in_bioinformatics.html)
  for an explanation.

- Cluster submission. Pypiper pipelines are scripts. You can use whatever system
  you want to run them on whatever computing resources you have. We have divided
  cluster resource management into a separate project called
  [looper](http://looper.readthedocs.io/).Pypiper  builds individual, single-sample
  pipelines that can be run one sample at a time.
  [Looper](http://looper.readthedocs.io/) then processes groups of samples,
  submitting appropriate pipelines to a cluster or server. The two projects are
  independent and can be used separately, keeping things simple and modular.


## Yet another pipeline system?

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
