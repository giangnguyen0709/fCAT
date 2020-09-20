#' The function move the fasta file and the annotation file of the interested
#' genome in to the equivalent folder in the core set to perfome the assessment
#' 
#' @param genome the path to the fasta file of the genome
#' @param fasAnno the path to the fas annotation file of the genome. If equal 
#' NULL, the function will compute the annotation
#' @param coreSet the path to the core set
#' @param process A logical option to determine if the function was used as a 
#' subfunction in the processCoreSet function
#' 
#' @return none
#' @export
placeSeed <- function(genome, fasAnno=NULL, coreSet, process=FALSE) {
  if (!endsWith(coreSet,"/")) {
    coreSet <- paste(coreSet, "/", sep = "");
  }
  
  splited <- strsplit(genome, "/", fixed=TRUE)[[1]];
  splited <- strsplit(splited[length(splited)], ".", fixed=TRUE)[[1]];
  genomeName <- splited[1];
  
  genomePath <- paste(coreSet, "genome_dir", sep = "");
  weightPath <- paste(coreSet, "weight_dir", sep = "");
  
  if (process == FALSE) {
    if (!is.null(fasAnno)) {
      R.utils::createLink(paste(weightPath, "/", genomeName, ".json", sep=""), 
                          fasAnno,
                          overwrite=TRUE);
    } else {
      command <- paste("annoFAS", "-i", genome, "-o", weightPath, 
                       "-n", genomeName, sep=" ");
      system(command);
    }
  }
  
  dir.create(paste(genomePath, "/", genomeName, sep=""));
  file.copy(genome, paste(genomePath, "/", genomeName, "/", genomeName, ".fa", sep=""));
}