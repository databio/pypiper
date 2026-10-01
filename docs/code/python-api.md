# Package `pypiper` Documentation

## Package Overview

The `pypiper` package provides a framework for building robust, restartable bioinformatics pipelines. It handles common pipeline tasks like checkpointing, logging, and resource monitoring.

### Key Features

- **Automatic Checkpointing**: Resume pipelines from where they left off
- **Resource Monitoring**: Track memory and CPU usage
- **Result Reporting**: Integrate with pipestat for standardized results
- **Container Support**: Run commands in Docker containers
- **Pipeline Management**: Built-in logging and status tracking

### Installation

```bash
pip install pypiper
```

### Quick Example

```python
from pypiper import PipelineManager

# Initialize a pipeline
pm = PipelineManager(
    name="my_pipeline",
    outfolder="results/"
)

# Run a command
pm.run("echo 'Hello, world!'")

# Stop the pipeline
pm.stop_pipeline()
```

## API Reference

### PipelineManager Class

The main class for building and managing pipelines:

::: pypiper.PipelineManager
    options:
      docstring_style: google
      show_source: true
      show_signature: true
      merge_init_into_class: true
