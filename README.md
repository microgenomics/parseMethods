![banner](https://raw.githubusercontent.com/microgenomics/tutorials/master/img/microgenomics.png)
# parseMethods Module
-------------------------

Welcome to the parseMethods module! This module takes a directory containing output files from several metagenomics software and creates a standard table (with the full lineage and read/proportions) that can be easily imported into R for downstream analyses.
parseMethods is part of SEPA, a modular pipeline for testing metagenomic profiling software using simulations. Currently, parseMethods can process result files from pathoscope2, metamix, sigma, metaphlan2 (under development), and constrains (metaphlan dependent)

## Dependencies

Bash (comes with Linux and MacOSX)
Rscript (R statistical framework) 3.0 and greater ([download here](https://www.r-project.org))
R packages rJava, methods, xlsxjars (install from R console by issuing `install.packages(c("rJava", "methods", "xlsxjars"))`)


## Usage

    bash parseMethods.bash --workpath ./example_files/ --cfile config.conf --local

Where `--workpath` is the path to the directory that contains your program-specific output files. In the example, the files within `./example_files/` come from [PathoScope2](https://github.com/PathoScope/PathoScope). `--cfile` is a plain-text configuration file where you can pass options to the program as to what piece of software the files come from, whether you want proportions instead of reads, etc. Finally, the `--local` option is used when you want to run parseMethods.bash independently from SEPA, the larger modular pipeline for evaluating metagenomic profiling software.


## Example output



## Citation

Hong C, Manimaran S, Shen Y, Perez-Rogers JF, Byrd AL, Castro-Nallar E, Crandall KA, Johnson WE.  
*PathoScope 2.0: a complete computational framework for strain identification in environmental or clinical sequencing samples* 
**Microbiome** 2014 Dec 1;2(1):1-5.
[PMID: 25225611](http://www.ncbi.nlm.nih.gov/pubmed/25225611)  

Morfopoulou S, Plagnol V.  
*Bayesian mixture analysis for metagenomic community profiling*  
**Bioinformatics** 2015 May 21:btv317.
[PMID: 26002885](http://www.ncbi.nlm.nih.gov/pubmed/26002885)

Ahn TH, Chai J, Pan C.  
*Sigma: Strain-level inference of genomes from metagenomic analysis for biosurveillance*  
**Bioinformatics** 2014 Sep 29:btu641.
[PMID: 25266224](http://www.ncbi.nlm.nih.gov/pubmed/25266224)

Truong DT, Franzosa EA, Tickle TL, Scholz M, Weingart G, Pasolli E, Tett A, Huttenhower C, Segata N.  
*MetaPhlAn2 for enhanced metagenomic taxonomic profiling* 
**Nature methods** 2015 Oct 1;12(10):902-3.
[PMID: 26418763](http://www.ncbi.nlm.nih.gov/pubmed/26418763)

Luo C, Knight R, Siljander H, Knip M, Xavier RJ, Gevers D.   
*ConStrains identifies microbial strains in metagenomic datasets*  
**Nature biotechnology** 2015 Oct 1;33(10):1045-52.
[PMID: 26344404](http://www.ncbi.nlm.nih.gov/pubmed/26344404)
