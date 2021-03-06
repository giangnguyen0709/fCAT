% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/processCoreSet.R
\name{createAnnotation}
\alias{createAnnotation}
\title{The function creates FAS annotation for the fasta file of a core gene 
in the core set. This annotation file will be saved in each core gene folder
and be used then to calculate the cut off value for each score mode.}
\usage{
createAnnotation(root, coreSet, coreGene)
}
\arguments{
\item{root}{The path to the core directory, where the core set is stored
within weight_dir, blast_dir, etc.}

\item{coreSet}{The name of the interested core set. The core directory can
contains more than one core set and the user must specify the interested
core set. The core set will be stored in the folder core_orthologs in
subfolder, specify them by the name of the subfolder}

\item{coreGene}{The ID of the core gene in the set}
}
\value{
none
}
\description{
The function creates FAS annotation for the fasta file of a core gene 
in the core set. This annotation file will be saved in each core gene folder
and be used then to calculate the cut off value for each score mode.
}
\examples{
coreFolder <- system.file("extdata", "sample", package = "fCAT")
createAnnotation(coreFolder, "test", "530670")
}
