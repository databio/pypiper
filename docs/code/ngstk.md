<script>
document.addEventListener('DOMContentLoaded', (event) => {
  document.querySelectorAll('h3 code').forEach((block) => {
    hljs.highlightBlock(block);
  });
});
</script>

<style>
h3 .content { 
    padding-left: 22px;
    text-indent: -15px;
 }
h3 .hljs .content {
    padding-left: 20px;
    margin-left: 0px;
    text-indent: -15px;
    martin-bottom: 0px;
}
h4 .content, table .content, p .content, li .content { margin-left: 30px; }
h4 .content { 
    font-style: italic;
    font-size: 1em;
    margin-bottom: 0px;
}

</style>


# Package `ngstk` Documentation

 Package exports 
## <a name="UnsupportedFiletypeException"></a> Class `UnsupportedFiletypeException`
Restrict domain of file types.


```python
def count_fail_reads(file_name, paired_end, prog_path)
```

Count the number of reads that failed platform/vendor quality checks.
#### Parameters:

- `file_name` (`str`):  name/path to file to examine
- `paired_end` (``):  this parameter is ignored; samtools automaticallycorrectly responds depending on the data in the bamfile; we leave the option here just for consistency, since all the other counting functions require the parameter; this makes it easier to swap counting functions during pipeline development.
- `prog_path` (`str`):  path to main program/tool to use for the counting


#### Returns:

- `int`:  count of failed reads




```python
def count_flag_reads(file_name, flag, paired_end, prog_path)
```

Counts the number of reads with the specified flag.
#### Parameters:

- `file_name` (`str`):  name/path to file to examine
- `flag` (`str int |`):  SAM flag value to be read
- `paired_end` (``):  this parameter is ignored; samtools automaticallycorrectly responds depending on the data in the bamfile; we leave the option here just for consistency, since all the other counting functions require the parameter; this makes it easier to swap counting functions during pipeline development.
- `prog_path` (`str`):  path to main program/tool to use for the counting


#### Returns:

- `str`:  terminal-like text output




```python
def count_lines(file_name)
```

Uses the command-line utility wc to count the number of lines in a file.

For MacOS, must strip leading whitespace from wc.
#### Parameters:

- `file_name` (`str`):  name of file whose lines are to be counted


#### Returns:

- `str`:  terminal-like text output




```python
def count_lines_zip(file_name)
```

Count number of lines in a zipped file.

This function eses the command-line utility wc to count the number of lines
in a file. For MacOS, strip leading whitespace from wc.
#### Parameters:

- `file_name` (`str`):  path to file in which to count lines


#### Returns:

- `str`:  terminal-like text output




```python
def count_reads(file_name, paired_end, prog_path)
```

Count reads in a file.

Paired-end reads count as 2 in this function.
For paired-end reads, this function assumes that the reads are split
into 2 files, so it divides line count by 2 instead of 4.
This will thus give an incorrect result if your paired-end fastq files
are in only a single file (you must divide by 2 again).
#### Parameters:

- `file_name` (`str`):  Name/path of file whose reads are to be counted.
- `paired_end` (`bool`):  Whether the file contains paired-end reads.
- `prog_path` (`str`):  path to main program/tool to use for the counting


#### Returns:

- `str`:  terminal-like text output (if input is SAM/BAM), or actualcount value (if input isn't SAM/BAM)




```python
def get_file_size(filename)
```

Get size of all files in gigabytes (Gb).
#### Parameters:

- `filename` (`str | collections.Iterable[str]`):  A space-separatedstring or list of space-separated strings of absolute file paths.


#### Returns:

- `float`:  size of file(s), in gigabytes.




```python
def get_input_ext(input_file)
```

Get the extension of the input_file.

This function assumes you're using .bam, .fastq/.fq, or .fastq.gz/.fq.gz.
#### Parameters:

- `input_file` (`str`):  name/path of file for which to get extension


#### Returns:

- `str`:  standardized extension


#### Raises:

- `ubiquerg.ngs.UnsupportedFiletypeException`:  if the given file nameor path has an extension that's not supported




```python
def is_fastq(file_name)
```

Determine whether indicated file appears to be in FASTQ format.
#### Parameters:

- `file_name` (`str`):  Name/path of file to check as FASTQ.


#### Returns:

- `bool`:  Whether indicated file appears to be in FASTQ format, zippedor unzipped.




```python
def is_gzipped_fastq(file_name)
```

Determine whether indicated file appears to be a gzipped FASTQ.
#### Parameters:

- `file_name` (`str`):  Name/path of file to check as gzipped FASTQ.


#### Returns:

- `bool`:  Whether indicated file appears to be in gzipped FASTQ format.




```python
def is_sam_or_bam(file_name)
```

Determine whether a file appears to be in a SAM format.
#### Parameters:

- `file_name` (`str`):  Name/path of file to check as SAM-formatted.


#### Returns:

- `bool`:  Whether file appears to be SAM-formatted




```python
def is_unzipped_fastq(file_name)
```

Determine whether indicated file appears to be an unzipped FASTQ.
#### Parameters:

- `file_name` (`str`):  Name/path of file to check as unzipped FASTQ.


#### Returns:

- `bool`:  Whether indicated file appears to be in unzipped FASTQ format.




```python
def parse_ftype(input_file)
```

Checks determine filetype from extension.
#### Parameters:

- `input_file` (`str`):  String to check.


#### Returns:

- `str`:  filetype (extension without dot prefix)


#### Raises:

- `UnsupportedFiletypeException`:  if file does not appear of asupported type




```python
def peek_read_lengths_and_paired_counts_from_bam(bam, sample_size)
```

Counting read lengths and paired reads in a sample from a BAM.
#### Parameters:

- `bam` (`str`):  path to BAM file to examine
- `sample_size` (`int`):  number of reads to look at for estimation


#### Returns:

- `defaultdict[int, int], int`:  read length observation counts, andnumber of paired reads observed


#### Raises:

- `OSError`: 




```python
def samtools_view(file_name, param, prog_path, postpend='')
```

Run samtools view, with flexible parameters and post-processing.

This is used to implement the various read counting functions.
#### Parameters:

- `file_name` (`str`):  name/path of reads tile to use
- `param` (`str`):  String of parameters to pass to samtools view
- `prog_path` (`str`):  path to the samtools program
- `postpend` (`str`):  String to append to the samtools command;useful to add cut, sort, wc operations to the samtools view output.


#### Returns:

- `str`:  terminal-like text output







*Version Information: `ngstk` v0.0.1pre2, generated by `lucidoc` v0.4.3*