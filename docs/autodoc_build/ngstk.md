# Package pypiper Documentation

## Class NGSTk
Class to hold functions to build command strings used during pipeline runs. Object can be instantiated with a string of a path to a yaml `pipeline config file`. Since NGSTk inherits from `AttMapEcho`, the passed config file and its elements will be accessible through the NGSTk object as attributes under `config` (e.g. `NGSTk.tools.java`). In case no `config_file` argument is passed, all commands will be returned assuming the tool is in the user's $PATH.

**Parameters:**

- `config_file` -- `str`:  Path to pipeline yaml config file (optional).
- `pm` -- `pypiper.PipelineManager`:  A PipelineManager with which to associate this toolkit instance;that is, essentially a source from which to grab paths to tools, resources, etc.


**Example(s):**

```console
    from pypiper.ngstk import NGSTk as tk
    tk = NGSTk()
    tk.samtools_index("sample.bam")
    # returns: samtools index sample.bam

    # Using a configuration file (custom executable location):
    from pypiper.ngstk import NGSTk
    tk = NGSTk("pipeline_config_file.yaml")
    tk.samtools_index("sample.bam")
    # returns: /home/.local/samtools/bin/samtools index sample.bam
```

