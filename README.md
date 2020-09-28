# fCAT

## Installation

To install *fCAT*, open R in your terminal

First, install devtools

`install.packages("devtools")`

Using *devtools* to install fCAT

`devtools::install_github("giangnguyen0709/fCAT")`

## Depencies

*fCAT* is depended on some tools.

### HaMStR-oneSeq (h1s)

> https://github.com/BIONF/HaMStR

### FAS

> https://github.com/BIONF/FAS

### Packages in R
> R.utils

> taxize

> EnvStats

## Usage

`fCAT::checkCompleteness(genome, fasAnno, root, coreSet, scoreMode, extend, redo, priorityList, cpu)`
* genome: The path to the genome fasta file
* fasAnno: The path to the FAS annotation file. If the user did not provide the tool will compute it.
* root: The path to the root folder
* coreSet: The name of the core set
* scoreMode: The mode to determines the method to assess the founded ortholog
* extend: A logical value to decide if the phylogenetic profile of the interested genome will be appended into the original phylogenetic profile
* redo: A logical value to decide if the tool will recheck completeness of a genome, whose ID already exists in the original phylogenetic profile
* priorityList: A vector, which contains ID of the training genomes.
* cpu: The number of the cores that h1s will use to run parallel

`fCAT::computeOriginal(root, coreSet)`

The function to compute the original phylogenetic profile. It is optional, the tool can still check the completeness of a genome, even the orginal phylogenetic profile was not computed

* root: The path to the root folder
* coreSet: The name of the core set
