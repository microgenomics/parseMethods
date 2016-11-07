![banner](https://raw.githubusercontent.com/microgenomics/tutorials/master/img/microgenomics.png)
# parseMethods Module
-------------------------

Welcome to the parseMethods module! This module takes a directory containing output files from several metagenomics software and creates a standard table (with the full lineage and read/proportions) that can be easily imported into R for downstream analyses.
parseMethods is part of SEPA, a modular pipeline for testing metagenomic profiling software using simulations. Currently, parseMethods can process result files from pathoscope2, metamix, sigma, metaphlan2 (under development), and constrains (metaphlan dependent)

## Requirements

Bash version 4 (comes with Linux and MacOSX)
Rscript (R statistical framework) 3.0 and greater ([download here](https://www.r-project.org))
R packages rJava, methods, xlsxjars (install from R console by issuing `install.packages(c("rJava", "methods", "xlsxjars"))`)


## Usage

    bash parseMethods.bash --workpath ./example_files/ --cfile config.conf --local

Where `--workpath` is the path to the directory that contains your program-specific output files. In the example, the files within `./example_files/` come from [PathoScope2](https://github.com/PathoScope/PathoScope). `--cfile` is a plain-text configuration file where you can pass options to the program as to what piece of software the files come from, whether you want proportions instead of reads, etc. Finally, the `--local` option is used when you want to run parseMethods.bash independently from SEPA, the larger modular pipeline for evaluating metagenomic profiling software.


## Example output

| Superkingdom | Phylum         | Class          | Order               | Family               | Genus             | Specie                     | Name                                 | ti      | ID1000-sam | ID150-sam | ID300-sam | ID75-sam |
|--------------|----------------|----------------|---------------------|----------------------|-------------------|----------------------------|--------------------------------------|---------|------------|-----------|-----------|----------|
| Bacteria     | Firmicutes     | Bacilli        | Bacillales          | Bacillaceae          | Bacillus          | Bacillus anthracis         | Bacillus anthracis str. CDC 684      | 568206  | 157496     | 155096    | 156413    | 157488   |
| Bacteria     | Firmicutes     | Bacilli        | Bacillales          | Bacillaceae          | Bacillus          | Bacillus anthracis         | Bacillus anthracis str. Ames         | 198094  | 1672       | 0         | 0         | 123      |
| Bacteria     | Firmicutes     | Bacilli        | Bacillales          | Bacillaceae          | Bacillus          | Bacillus anthracis         | Bacillus anthracis str. Sterne       | 260799  | 10         | 0         | 0         | 0        |
| Bacteria     | Firmicutes     | Bacilli        | Bacillales          | Bacillaceae          | Bacillus          | Bacillus anthracis         | Bacillus anthracis str. H9401        | 768494  | 0          | 0         | 0         | 1437     |
| Bacteria     | Firmicutes     | Bacilli        | Lactobacillales     | Streptococcaceae     | Streptococcus     | Streptococcus pneumoniae   | Streptococcus pneumoniae TCH8431/19A | 525381  | 5556       | 5362      | 5244      | 5266     |
| Bacteria     | Firmicutes     | Bacilli        | Lactobacillales     | Streptococcaceae     | Streptococcus     | Streptococcus pneumoniae   | Streptococcus pneumoniae ST556       | 1130804 | 20         | 0         | 0         | 0        |
| Bacteria     | Actinobacteria | Actinobacteria | Propionibacteriales | Propionibacteriaceae | Propionibacterium | Propionibacterium acnes    | Propionibacterium acnes KPA171202    | 267747  | 5544       | 5373      | 5320      | 5504     |
| Bacteria     | Firmicutes     | Bacilli        | Lactobacillales     | Lactobacillaceae     | Lactobacillus     | Lactobacillus crispatus    | Lactobacillus crispatus ST1          | 748671  | 0          | 5689      | 5586      | 5586     |
| Bacteria     | Actinobacteria | Actinobacteria | Corynebacteriales   | Mycobacteriaceae     | Mycobacterium     | Mycobacterium tuberculosis | Mycobacterium tuberculosis CDC1551   | 83331   | 0          | 12958     | 12758     | 12946    |
| Bacteria     | Firmicutes     | Bacilli        | Lactobacillales     | Streptococcaceae     | Streptococcus     | Streptococcus anginosus    | Streptococcus anginosus C238         | 862971  | 0          | 5850      | 5648      | 5887     |
| Bacteria     | Firmicutes     | Bacilli        | Bacillales          | Staphylococcaceae    | Staphylococcus    | Staphylococcus epidermidis | Staphylococcus epidermidis RP62A     | 176279  | 0          | 4349      | 4248      | 4314     |
| Bacteria     | Firmicutes     | Bacilli        | Lactobacillales     | Lactobacillaceae     | Lactobacillus     | Lactobacillus rhamnosus    | Lactobacillus rhamnosus ATCC 8530    | 1088720 | 0          | 741       | 624       | 707      |
| Bacteria     | Actinobacteria | Actinobacteria | Propionibacteriales | Propionibacteriaceae | Propionibacterium | Propionibacterium acnes    | Propionibacterium acnes 6609         | 1031709 | 0          | 0         | 0         | 3        |
| Bacteria     | Actinobacteria | Actinobacteria | Corynebacteriales   | Mycobacteriaceae     | Mycobacterium     | Mycobacterium tuberculosis | Mycobacterium tuberculosis KZN 605   | 478435  | 0          | 0         | 0         | 164      |
| Bacteria     | Actinobacteria | Actinobacteria | Corynebacteriales   | Mycobacteriaceae     | Mycobacterium     | Mycobacterium tuberculosis | Mycobacterium tuberculosis H37Ra     | 419947  | 0          | 0         | 0         | 9        |
| Bacteria     | Actinobacteria | Actinobacteria | Corynebacteriales   | Mycobacteriaceae     | Mycobacterium     | Mycobacterium tuberculosis | Mycobacterium tuberculosis H37Rv     | 83332   | 0          | 0         | 0         | 2        |

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
