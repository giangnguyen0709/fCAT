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

`fCAT::computeOriginal(root, coreSet)`
