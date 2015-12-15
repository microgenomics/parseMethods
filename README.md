# parseMethods Module

Welcome to the parseMethods module! This module takes a directory with output from several metagenomics software and creates a standard table that can be easily imported into R for downstream analyses.
parseMethods is part of SEPA, a multi pipe for testing metagenomic profiling software using simulations. currently, parseMethods can process result files from pathoscope2, metamix, sigma, metaphlan (under development), and constrains (metaphlan dependent)

## Usage

    bash parseMethods.bash --workpath ./example_files/ --cfile config.conf --local

Where `--workpath` is the path to the directory that contains your program-specific output files. In the example, the files within `./example_files/` come from [PathoScope2](https://github.com/PathoScope/PathoScope)
