# parseMethods Module

Welcome to the parseMethods module! This module takes a directory containing output files from several metagenomics software and creates a standard table that can be easily imported into R for downstream analyses.
parseMethods is part of SEPA, a modular pipeline for testing metagenomic profiling software using simulations. Currently, parseMethods can process result files from pathoscope2, metamix, sigma, metaphlan2 (under development), and constrains (metaphlan dependent)

## Usage

    bash parseMethods.bash --workpath ./example_files/ --cfile config.conf --local

Where `--workpath` is the path to the directory that contains your program-specific output files. In the example, the files within `./example_files/` come from [PathoScope2](https://github.com/PathoScope/PathoScope). `--cfile` is a plain-text configuration file where you can pass options to the program as to what piece of software the files come from, whether you want proportions instead of reads, etc. Finally, the `--local` option is used when you want to run parseMethods.bash independently from SEPA, the larger modular pipeline for evaluating metagenomic profiling software.
