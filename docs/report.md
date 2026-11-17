# Reporting statistics

One of the most useful features of pypiper is the `report_result` function. This function provides a way to record small-scale results, like summary statistics. It standardizes the output so that universal tools can be built to process all the pipeline results from any pipeline, because the results are all reported in the same way.

When you call `pm.report_result(key, value)`, pypiper simply writes the key-value pair to a `tsv` file (`stats.tsv`) in the pipeline output folder. These `stats.tsv` files can then later be read and aggregated systematically by other tools, such as `looper summarize`.

## Reporting objects

**Note**: Reporting objects will be deprecated in a future release. It is recommended to use `report_result`.

Starting in version 0.8, pypiper now implements a second reporting function, `report_object`. This is analogous to the `report_result` function, but instead of reporting simple key-value pairs, it lets you record any produced file as an output. Most commonly, this is used to record figures (PDFs, PNGs, etc.) produced by the pipeline. It can also be used to report other files, like HTML files.

Pypiper writes results to `objects.tsv`, which can then be aggregated for project-level summaries of plots and other pipeline result files.


## Re-using previously reported results

We frequently want to use the `report_result` capability in `follow` functions. It's a convenient place to do something like count or assess the result of a long-running command, and then report some summary statistic on it. One potential hangup with this strategy is dealing with secondary results after a pipeline is interrupted and restarted. By secondary result, I mean one that requires knowing the value of an earlier result. For example, if you want to compute the **percentage of reads that aligned**, you need to first know the **total reads** -- but what if your pipeline got interrupted and calculation of **total reads** happened in an earlier pipeline run?

To solve this issue, Pypiper has a neat function called `get_stat` that lets you retrieve any value you've reported with `report_result` so you could use it to calculate statistics elsewhere in the pipeline. It will retrieve this either from memory, if the calculation of that result happened during the current pipeline run, or from the `stats.tsv` file, if the result was reported by an earlier run (or even another pipeline). So you could, in theory, calculate statistics based on results across pipelines.

An example for how to use this is how we handle calculating the alignment rate in an NGS pipeline:

```{python}
x = myngstk.count_mapped_reads(bamfile, args.paired_end)
pm.report_result("Aligned_reads", x)
rr = float(pm.get_stat("Raw_reads"))
pm.report_result("Alignment_rate", round((rr * 100 / float(x), 3))
```

Here, we use `get_stat` to grab a result that we reported previously (with `report_result`), when we counted the number of `Raw_reads` (earlier in the pipeline). We need this after the alignment to calculate the alignment rate. Later, now that we've reported `Alignment_rate`, you could harvest this stat again for use with `pm.get_stat("Alignment_rate")`. This is useful because you could put this block of code in a `follow` statement so it may not be executed, but you can still grab a reported result like this even if the execution happened outside of the current pipeline run; you'd only have to do the calculation once.
