% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/placeSeed.R
\name{placeSeed}
\alias{placeSeed}
\title{This function take the path to the interested genome and its FAS 
annotation (if FAS annotation was not provided, the function will call 
annoFAS tool to compute it) as the input. The function will make a copy of 
the genome fasta file and a symbolic link (not in case annoFAS is called) 
and arrange them into the query_taxon folder and into the weight_dir folder 
of the core directory.}
\usage{
placeSeed(genome, fasAnno = NULL, root, process = FALSE, weightDir = NULL)
}
\arguments{
\item{genome}{the path to the fasta file of the genome}

\item{fasAnno}{the path to the fas annotation file of the genome. If equal
NULL, the function will compute the annotation, and arrange it into the
weight_dir folder of the core directory}

\item{root}{The path to the core directory, where the core set is stored
within weight_dir, blast_dir, etc.}

\item{process}{A logical option to determine if the function was used as a
modul in the processCoreSet function}

\item{weightDir}{The user can replace the weight_dir folder in the core
directory by specifying the path to the replacing folder in this argument}
}
\value{
none
}
\description{
This function take the path to the interested genome and its FAS 
annotation (if FAS annotation was not provided, the function will call 
annoFAS tool to compute it) as the input. The function will make a copy of 
the genome fasta file and a symbolic link (not in case annoFAS is called) 
and arrange them into the query_taxon folder and into the weight_dir folder 
of the core directory.
}
\examples{
## Take the demo data
coreFolder <- system.file("extdata", "sample", package = "fCAT")
genome <- system.file("extdata", "HUMAN@9606@3.fa", package = "fCAT")
fasAnno <- system.file("extdata", "HUMAN@9606@3.json", package = "fCAT")
## Place seed
placeSeed(genome, fasAnno, coreFolder)
## Test if the seed was place
fastaPath <- paste(coreFolder, "/query_taxon/HUMAN@9606@3/HUMAN@9606@3.fa",
    sep = ""
)
annoPath <- paste(coreFolder, "/weight_dir/HUMAN@9606@3.json", sep = "")
if (file.exists(fastaPath) && file.exists(annoPath)) {
    print("Seed is placed")
}

## Delete seed
fastaFolder <- paste(coreFolder, "/query_taxon/HUMAN@9606@3", sep = "")
annoPath <- paste(coreFolder, "/weight_dir/HUMAN@9606@3.json", sep = "")
unlink(fastaFolder, recursive = TRUE)
file.remove(annoPath)
}
