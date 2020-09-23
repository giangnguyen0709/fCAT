#' In the scoreMode "busco" the tool must calculate the length of the founded
#' ortholog to perfome the assessment. After HaMStR founded an ortholog the 
#' function will be called to update the length of the ortholog in phylogenetic 
#' profile
#' 
#' @param root the path to the root folder
#' @param coreSet the core set name
#' @param coreGene the ID of the core gene
#' @param genomeName the genome ID of the interested genome
#' 
#' @return none
#' @export
updateLength <- function(root, coreSet, coreGene, genomeName) {
  if (!endsWith(root, "/")) {
    root <- paste(root, "/", sep="");
  }
  exFasta <- readLines(paste(root, "phyloprofile", "/", coreSet, "/", "busco", 
                             "/", "hamstrout", "/", coreGene, "/",
                             coreGene, ".extended.fa", sep=""));
  exFasta <- extractFasta(exFasta, genomeName)
  i <- 1:(length(exFasta));
  i <- i[i %% 2 == 0]
  orthoLength <- unlist(lapply(i, 
                               function(i, exFasta) {
                                 return(nchar(exFasta[i]))
                               },
                               exFasta=exFasta));
  pp <- read.table(paste(root, "phyloprofile", "/", coreSet, "/", "busco", "/", 
                         "hamstrout", "/", coreGene, "/",
                         coreGene, ".phyloprofile", sep=""),
                   header=TRUE,
                   sep="\t");
  pp <- cbind(pp, length=orthoLength);
  write.table(pp,
              paste(root, "phyloprofile", "/", coreSet, "/", "busco", "/", 
                    "hamstrout", "/", coreGene, "/",
                    coreGene, ".phyloprofile", sep=""),
              sep="\t",
              quote=FALSE,
              row.names=FALSE);
}