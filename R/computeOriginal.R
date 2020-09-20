#' The function compute the initial original phylogenetic profile of all the 
#' taxa in the score set
#' 
#' @param coreSet the path to the core set
#' @param scoreMode the mode determines the way to assess the founded ortholog
#' @param cpu determines the number of the cores
#' 
#' @return none
#' @export
computeOriginal <- function(coreSet, scoreMode, cpu=4) {
  if (!endsWith(coreSet, "/")) {
    coreSet <- paste(coreSet, "/", sep="");
  }
  
  blastDir <- paste(coreSet, "blast_dir", sep="");
  weightDir <- paste(coreSet, "weight_dir", sep="");
  
  for (genome in list.dirs(blastDir, full.names=FALSE, recursive=FALSE)) {
    genomeFasta <- paste(blastDir, "/", genome, "/", genome, ".fa", sep="");
    fasAnno <- paste(weightDir, "/", genome, ".json", sep="");
    computeReport(genomeFasta, fasAnno, coreSet,
                  TRUE, scoreMode, priorityList=c(genome), cpu, computeOri=TRUE);
  }
}