% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getResult.R
\name{translateReport}
\alias{translateReport}
\title{Translate the report table into a frequence table, which tell the user, 
how many "dissimilar", "simillar", "missing", etc.}
\usage{
translateReport(genomeID, report, scoreMode)
}
\arguments{
\item{genomeID}{the genome ID of the interested genome}

\item{report}{the report of the completeness of the interested genome based
on its phylogenetic profile}

\item{scoreMode}{the mode determines the method to scoring the founded
ortholog and how to classify them. Choices: 1, 2, 3, "busco"}
}
\value{
A frequency table of the completeness of the interested genome in
data.frame
}
\description{
Translate the report table into a frequence table, which tell the user, 
how many "dissimilar", "simillar", "missing", etc.
}
\examples{
## Create a pseudo phylogenetic profile table
geneID <- c("530670", "530730")
ncbiID <- c("ncbi9606", "ncbi9606")
orthoID <- c("530670|HUMAN@9606@3|Q16526|1", "530730|HUMAN@9606@3|P05091|1")
FAS_F <- c(1, 1)
FAS_B <- c(1, 1)
pp <- data.frame(geneID, ncbiID, orthoID, FAS_F, FAS_B)
coreFolder <- system.file("extdata", "sample", package = "fCAT")
## Translate phylogenetic profile to a detailed report
report <- reportSingle(pp, coreFolder, "test", 2, c("HUMAN@9606@1"))
## Translate the report table to a frequency table
frequencyTable <- translateReport("HUMAN@9606@3", report, 2)
print.data.frame(frequencyTable)
}
