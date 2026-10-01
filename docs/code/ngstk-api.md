# NGSTk API Documentation

## Overview

NGSTk (Next-Generation Sequencing Toolkit) is a toolkit class that provides helper functions for building command strings used in NGS pipelines. It can be configured with a YAML configuration file to specify custom tool paths, or it will use tools from the system PATH.

### Key Features

- **Command Building**: Generate command strings for common NGS tools
- **Configuration Management**: Use custom tool paths via YAML config
- **Tool Integration**: Built-in support for common tools like samtools, bedtools, etc.
- **Pipeline Integration**: Works seamlessly with PipelineManager

### Installation

NGSTk is included with pypiper:

```bash
pip install pypiper
```

### Quick Example

```python
from pypiper.ngstk import NGSTk

# Initialize NGSTk
tk = NGSTk()

# Generate a command
cmd = tk.samtools_index("sample.bam")
# Returns: "samtools index sample.bam"
```

## API Reference

### NGSTk Class

::: pypiper.ngstk.NGSTk
    options:
      docstring_style: google
      show_source: true
      show_signature: true
      merge_init_into_class: true
