% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getResult.R
\name{assessBusco}
\alias{assessBusco}
\title{The function to classify the founded ortholog based on its length}
\usage{
assessBusco(orthoLength, root, coreSet, coreGene, f)
}
\arguments{
\item{orthoLength}{the length of the ortholg sequence}

\item{root}{The path to the core directory, where the core set is stored
within weight_dir, blast_dir, etc.}

\item{coreSet}{The name of the interested core set. The core directory can
contains more than one core set and the user must specify the interested
core set. The core set will be stored in the folder core_orthologs in
subfolder, specify them by the name of the subfolder}

\item{coreGene}{the ID of the core gene}

\item{f}{the appearance frequence of the core gene in the phylogenetic
profile}
}
\value{
the status of the core gene, it can be "fragmented", "complete",
"duplicated, complete" or "duplicated, fragmented"
}
\description{
The function to classify the founded ortholog based on its length
}
\examples{
coreFolder <- system.file("extdata", "sample", package = "fCAT")
status <- assessBusco(600, coreFolder, "test", "530670", f = 1)
print(status)
}
