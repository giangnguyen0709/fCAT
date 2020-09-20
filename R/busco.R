#' In the scoreMode "busco" the tool must calculate the length of the founded
#' ortholog to perfome the assessment. After HaMStR founded an ortholog the 
#' function will be called to update the length of the ortholog in phylogenetic 
#' profile
#' 
#' @param coreSet the path to the core set
#' @param coreGene the ID of the core gene
#' 
#' @return none
#' @export
updateLength <- function(coreSet, coreGene) {
  if (!endsWith(coreSet, "/")) {
    coreSet <- paste(coreSet, "/", sep="");
  }
  exFasta <- readLines(paste(coreSet, "phyloprofile", "/", "busco", "/", 
                             "hamstrout", "/", coreGene, "/",
                             coreGene, ".extended.fa", sep=""));
  i <- 1:(length(exFasta));
  i <- i[i %% 2 == 0]
  orthoLength <- unlist(lapply(i, 
                               function(i, exFasta) {
                                 return(nchar(exFasta[i]))
                               },
                               exFasta=exFasta));
  pp <- read.table(paste(coreSet, "phyloprofile", "/", "busco", "/", 
                         "hamstrout", "/", coreGene, "/",
                         coreGene, ".phyloprofile", sep=""),
                   header=TRUE,
                   sep="\t");
  pp <- cbind(pp, length=orthoLength);
  write.table(pp,
              paste(coreSet, "phyloprofile", "/", "busco", "/", 
                    "hamstrout", "/", coreGene, "/",
                    coreGene, ".phyloprofile", sep=""),
              sep="\t",
              quote=FALSE,
              row.names=FALSE);
}