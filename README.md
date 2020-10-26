# fCAT

fCAT is a feature-aware completeness assessment tool, which provide a solution to assess the completeness of a genome assembly or a newly sequenced genome, based on the ortholog prediction with reciprocity criterion and the feature architecture similarity between the ortholog sequences of the interested genome and the training sequences.

## Installation

To install *fCAT*, open R in your terminal

Using *devtools* to install fCAT

```
if (!requireNamespace("devtools"))
    install.packages("devtools")
devtools::install_github("giangnguyen0709/fCAT")
```

## Usage

### checkCompleteness

The function to check the completeness of an interested genome

* **genome**: The path to the genome fasta file. The name of the fasta file must be SPECIES@ncbiID@version, for example: HUMAN@9606@3
* **fasAnno**: The path to the fas annotation file (must has the same name as the genome fasta file). It can equal NULL. If fasAnno equal NULL, fCAT will compute the FAS annotation for the genome fasta file.
* **coreDir**: The path to the core directory, where the core set is stored within weight_dir, blast_dir, etc.
* **coreSet**: The name of the interested core set. The core directory can contains more than one core set and the user must specify the interested core set. The core set will be stored in the folder core_orthologs in subfolder, specify them by the name of the subfolder
* **extend**: The output folder will contain a folder named phyloprofile. The folder will contains the phyloprofiles with the name corresponding to each score mode. If extend is set to TRUE, the phylogenetic profile of the interested genome will be append to the phyloprofile file in this folder. If it exists already a genome with the same name as the interested one, fCAT will print a message to user to ask if they really want to clear the old phylogenetic profile of the genome to append the new one.
* **redo**: The score modes of fCAT reuse the fDOG's output and score mode 2 and 3 even reuses the phyloprofile output of each other to create the assessment report of an interested genome. For some reasons, the user can want to remove the old data, and recheck the completeness of the genome from the beginning step (from the ortholog prediction step). To avoid the conflict data of the old data with the new one, fCAT must remove all the old data of the genome, not only the data of the current runing score mode but the interested genome's data of all 4 modes in fCAT and also the phylogenetic profile of the genome must be removed from the extended phyloprofile file. To do this, user can set redo = TRUE. A message will be printed to make sure that the user really want to clean the old data.
* **scoreMode**: the mode determines the method to scoring the founded ortholog and how to classify them. Choices: 1, 2, 3, "len"
* **refSpecList**: A list contains one or many genome ID of the genomes, which were used to build the core set. The genome ID of this list will be stored with an priority order, the tool look at into the fasta file of each core group and uses the priority order to determine the references species for each core group. 
* **cpu**: Optional, by default is 4. determines the cores that fDOG and fdogFAS will uses to be run parallel
* **blastDir**: Optional. The user can replace the blast_dir folder in the core directory by specifying it in this argument. By default is NULL
* **weightDir**: Optional. The user can replace the weight_dir folder in the core directory by specifying it in this argument. By default is NULL
* **cleanup**: Optional, if cleanup is set to TRUE, the fDOG output and the phyloprofile output will be removed
* **output** The path to the location, where the output directory will be stored. By default it is the working directory.

The function returns two reports. A detailed report of the completeness of the interested genome and a frequency table of all taxa, which were checked completeness with fCAT with option extend = TRUE. The frequency table show how many core genes "similar", "dissimilar", "duplicated", "missing" and "ignored" in each taxon.


```
genome <- "/path/to/query/genome.fa"
fasAnno <- "/path/to/fas/annotation/genome.json"
coreDir <- "/path/to/the/core/directory"
coreSet <- "name of the core set"
extend <- TRUE #by default is FALSE
redo <- TRUE #by default is FALSE
scoreMode <- 2 #Choices: 1,2,3, "len"
priorityList <- c("HUMAN@9606@1", "ECOLI@511145@1") 
cpu <- 4
blastDir <- "/path/to/blast_dir" #Optional
weightDir <- "/path/to/weight_dir" #Optional
cleanup <- TRUE #by default is FASLE
output <- "path/to/location/to/save/output"

checkCompleteness <- function(genome, fasAnno, coreDir, coreSet, extend, redo, scoreMode, refSpecList, cpu, blastDir, weightDir, cleanup, output)

```

### computeOriginal

The function to compute the original phylogenetic profile, which will contains the phylogenetic profile of all core taxa of the core set. This phylogenetic profile can be used to assess the completeness of the core taxa and their'completeness will be reported together with the interested genome in the frequency table. It is optional, the tool can still check the completeness of a genome, even the orginal phylogenetic profile was not computed

* **coreDir**: The path to the core directory, where the core set is stored within weight_dir, blast_dir, etc.
* **coreSet**: The name of the interested core set. The core directory can contains more than one core set and the user must specify the interested core set. The core set will be stored in the folder core_orthologs in subfolder, specify them by the name of the subfolder
* **scoreMode**: the mode determines the method to scoring the founded ortholog and how to classify them. Choices: 1, 2, 3, "busco"
* **cpu**: Optional, by default is 4. Determines the cores that fDOG and fdogFAS will uses to be run parallel
* **blastDir**: Optional. The user can replace the blast_dir folder in the core directory by specifying it in this argument. By default is NULL
* **weightDir**: Optional. The user can replace the weight_dir folder in the core directory by specifying it in this argument. By default is NULL

This function will check completeness for all core genomes in the core set and store the output in the output's folder

```
coreDir <- "/path/to/the/core/directory"
coreSet <- "name of the core set"
scoreMode <- 2 #Choices: 1,2,3, "busco"
cpu <- 4

fCAT::computeOriginal(coreDir, coreSet, scoreMode, cpu)
```

### processCoreSet

The function calculate all cutoff values for all mode in the set. For score mode 1 it will calculate the avarage of all vs all FAS scores between the training sequences in the core gene. For score mode 2 it will calculate the avarage of the FAS score between each sequence against all training sequences in the core gene. The scores will be writen in a table with a column is the ID of the sequences and a column is the corresponding value. For score mode 3, the function will calculate the avarage of 1 vs all FAS scores for each training sequence in the core gene. The avarages build a distribution, the function will calculate the confidence interval of this distribution and write the upper value and the lower value of the interval in a file in the core gene folder.

* **coreDir**: The path to the core directory, where the core set is stored within weight_dir, blast_dir, etc.
* **coreSet**: The name of the interested core set. The core directory can contains more than one core set and the user must specify the interested core set. The core set will be stored in the folder core_orthologs in subfolder, specify them by the name of the subfolder

```
coreDir <- "/path/to/the/core/directory"
coreSet <- "name of the core set"

fCAT::processCoreSet(coreDir, coreSet)
```

## Examples

The test data is the eukaryota_busco set, which can be downloaded in https://applbio.biologie.uni-frankfurt.de/download/core-sets/BUSCO_Eukaryota/

The core set folder has some conflicts with the input of fCAT, which must be removed first. All the core gene folders in the folder core_orthologs must be contained in a subfolder (The name of the subfolder is the core set argument of fCAT) and this subfolder must be stored in core orthologs. In this document I will set the name of this subfolder eukaryota_busco. In the blast dir folder of CRYNE, the symbolic link of the fasta file of CRYNE was directed by a mistake to the symbolic link of CHRLE. This must be corrected before testing

The folder weight_dir of the core set contains the the xml files, which can not be run with fCAT. Please download the annotation files of the core set from https://drive.google.com/file/d/113MBwT1n7E64Xk54Ul-_r82aEK2v8jdl/view?usp=sharing and replace the xml files with the json files 

In all following examples, I assumed that I has a genome fasta file and its FAS annotation file, which named HUMAN@9606@3.fa and HUMAN@9606@3.json, the core folder named eukaryota_busco, the core set named eukaryota_busco and all this data is placed in the home folder. You can replace them by your corresponding path and names

```
genome <- "/home/user/HUMAN@9606@3.fa"
fasAnno <- "/home/user/HUMAN@9606@3.json"
coreDir <- "/home/user/eukaryota_busco"
coreSet <- "eukaryota_busco"
extend <- TRUE
refSpecList <- c("HOMSA@9606@2")
scoreMode <- 1
cpu <- 4

fCAT::checkCompleteness(genome = genome, fasAnno = fasAnno, coreDir = coreDir, coreSet = coreSet, extend = extend, refSpecList = refSpecList, scoreMode = scoreMode, cpu = cpu)
```

The report will be storede by default in /home/user/eukaryota_busco/output/eukaryota_busco/1/report

```
coreDir <- "/home/user/eukaryota_busco"
coreSet <- "eukaryota_busco"
scoreMode <- 1
cpu <- 4

fCAT::computeOriginal(coreDir = coreDir, coreSet = coreSet, scoreMode = scoreMode, cpu = cpu)
```

The phylogenetic profile of all core taxa will be computed and be stored by default in /home/user/eukaryota_busco/output/eukaryota_busco/1

```
coreDir <- "/home/user/eukaryota_busco"
coreSet <- "eukaryota_busco"

fCAT::processCoreSet(coreDir = coreDir, coreSet = coreSet)
```

The function will calculate all cutoff values for all core genes in the set and write them in a text file, which will be stored in /home/user/eukaryota_busco/core_orthologs/eukaryota_busco/core_gene/fas_dir/score_dir

## Depencies

*fCAT* is depended on some tools.

### fDOG

> https://github.com/BIONF/fDOG

### FAS

> https://github.com/BIONF/FAS

### Packages in R
> R.utils

> taxize

> EnvStats

## Bugs

Any bug reports or comments, suggestions are highly appreciated. Please open an issue on GitHub or be in touch via email.

## Contact

Thanh-Giang Nguyen > giangnguyen0709@gmail.com
