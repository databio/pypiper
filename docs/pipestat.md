# Pipestat

Starting with pypiper v0.13.0 [pipestat](http://pipestat.databio.org) is the recommended way of reporting pipeline statistics.
You can browse the pipestat documentation to learn more about it, but briefly pipestat is a tool that standardizes reporting of pipeline results. It provides 1) a standard specification for how pipeline outputs should be stored; and 2) an implementation to easily write results to that format from within Python or from the command line.

## Advancements

There are a multiple advantages of using piestat instead of the current pieline results reporiting system:

1. **Database results storage:** the results can be stored either in a database or a YAML-formatted results file. This way a pypiper pipeline running in an emphemeral compute environment can report the results to the database and exit. No need to sync the results with a central results storage.
2. **Strict and clear results definition:** all the results that can be reported by a pipeline run *must* be pre-defined in a [pipestat results schema](http://pipestat.databio.org/en/latest/pipestat_specification/#pipestat-schema-format) that in a simplest case just indicates the result's type. This presents piepstat clients with the possibility to *reliably* gather all the possible results and related metadata.
3. **On-the-fly results validation:** the schema is used to validate and/or convert the reported result to a strictly determined type, which makes the connection of pypiper with downstream pipeline results processing software seamless.
4. **Unified, pipeline-agnostic results interface:** other pipelines, possibly created with different pipeline frameworks, can read and write results via Python API or command line interface. This feature significantly incerases your pipeline interoperability.

## Setup

In order to start reporting results with pipestat in your pipeline all you need to do is:  

### Define a [pipestat resuts schema](http://pipestat.databio.org/en/latest/pipestat_specification/#pipestat-schema-format)

```yaml
my_int_result:
  type: integer
  description: "This is my first result"
my_str_result:
  type: string
```

### Pass the pipestat results schema to the `PipelineManager` object constructor

```python
pm = pypiper.PipelineManager(
  name="hello_pypiper",
  outfolfer="$HOME/hello_pypiper",
  pipestat_schema="pipestat_results_schema.yaml",
) 
```

If `pipestat_schema` argument is not provided, by default pypiper will look for a `pipestat_results_schema.yaml` file next to the pipeline Python script.

### Use `pipestat` property of the `PipelineManager` object to report/retrieve results. See usage for more details

And in the simplest case... that's it! Pypiper *by default* will use a YAML-formated file to store the reported results in the selected `outfolder`.

### Advanced features

Pypiper-pipestat integration really shines when more advanced features are used. Here's how to set them up.

#### Use a database to store reported results*

In order to establish a database connection pipestat requires few pieces of information, which *must* be provided in a [pipestat configuration file](http://pipestat.databio.org/en/latest/config/) passed to the `PipelineManager` constructor.

This is an example of such a file:

```yaml
database:
  name: pypiper # database name
  user: pypiper # database user name
  password: pypiper # database password
  host: localhost # database host address
  port: 5433 # port the database is running on
  dialect: postgresql # type of the databse 
  driver: psycopg2 # driver to use to communicate
```

For reference, here is a Docker command that would run a PostgreSQL instance that could be used to store the pipeline results when configured with with the configuration file above:

```console
docker volume create postgres-data

docker run -d --name pypiper-postgres \
-p 5432:5433 -e POSTGRES_PASSWORD=pypiper \
-e POSTGRES_USER=pypiper -e POSTGRES_DB=pypiper \
-v postgres-data:/var/lib/postgresql/data postgres
```

#### Highlight results

The pipestat results schema can include any number of additional attributes for results. An example of that is *results highlighting*. 

When a `highlight: true` attribute is included attribute under result identifier in the schema file the highlighted results can be later retrieved by pipestat clients via `PipelineManager.pipestat.highlighted_results` property, which simply returns a list of result identifiers. to be presented in a special way.

### Usage

Since a pipeline run-specific `PipestatManager` instance is attached to the `PipelineManager` object all the public pipestat API can be used. Please refer to the [pipestat API documentation](http://pipestat.databio.org/en/latest/autodoc_build/pipestat/) to read about all the currently available features.

Here we present the most commonly used features:

- results reporting

*report a result, convert to schema-defined type and overwrite previously reported result*

```python
results = {
  "my_int_result": 10,
  "my_str_result": "test"
}
pm.pipestat.report(
  values=results,
  strict_type=True,
  force_overwrite=True
)
```

- results retrieval

```python
pm.pipestat.retrieve(result_identifier="my_int_result")
```

- results schema exploration

```python
pm.pipestat.schema
```

- exploration of canonical [jsonschema](https://json-schema.org/) representation of result schemas

```python
pm.pipestat.result_schemas
```
