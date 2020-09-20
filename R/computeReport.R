#' Create the report of the completeness of the interested genome with the other
#' genome in the original phylogenetic profile
#' 
#' @param report the data frame that contains the report
#' @param coreSet the path to the coreSet
#' @param scoreMode the mode determines the way to assess the founded ortholog
#' 
#' @return none
#' @export
printReport <- function(report, coreSet, scoreMode) {
  splited <- strsplit(coreSet, "/", fixed=TRUE)[[1]];
  setName <- splited[length(splited)];
  
  reportFile <- paste(coreSet, "phyloprofile", "/", as.character(scoreMode), "/", 
                      setName, ".report", sep="");
  
  if (!file.exists(reportFile)) {
    write.table(report, 
                reportFile,
                row.names=FALSE,
                quote=FALSE,
                sep="\t");
  } else {
    write.table(report, 
                reportFile,
                row.names=FALSE,
                col.names=FALSE,
                append=TRUE,
                quote=FALSE,
                sep="\t");
  }
} 

#' Add the priority list of the interested genome to the priority file. If the 
#' file doesnt exist it will be created
#' 
#' @param genomeID the genomeID of the interested genome
#' @param priorityList the priority list
#' @param coreSet the path to the core set
#' @param scoreMode the mode determines the method to assess the founded 
#' ortholog
#' 
#' @return none
#' @export
printPriority <- function(genomeID, priorityList, coreSet, scoreMode) {
  table <- data.frame(genomeID=c(genomeID),
                      priority_list=c(paste(priorityList, collapse=",")));
  priorityFile <- paste(coreSet, "phyloprofile", "/", 
                        as.character(scoreMode), "/", "priority.list", sep="");
  
  if (!file.exists(priorityFile)) {
    write.table(table,
                paste(coreSet, "phyloprofile", "/", 
                      as.character(scoreMode), "/", "priority.list", sep=""),
                row.names=FALSE,
                quote=FALSE,
                sep="\t");
  } else {
    write.table(table,
                paste(coreSet, "phyloprofile", "/", 
                      as.character(scoreMode), "/", "priority.list", sep=""),
                row.names=FALSE,
                col.names=FALSE,
                append=TRUE,
                quote=FALSE,
                sep="\t");
  }
}

#' The function create the report of the completeness of an interested genome
#' 
#' @param genome the path to the fasta file of the interested genome
#' @param fasAnno the path to the json file of the fas annotation of the genome
#' @param coreSet the path to the core set
#' @param extend the logical argument to decision if the information of the 
#' interested genome will be appended to the original files
#' @param scoreMode the mode determines the method to assess the founded 
#' ortholog
#' @param priorityList the priority list to determines the references species
#' If equal NULL the tool will create the list based on the taxonomic distance
#' @param cpu determines the number of the cores
#' @param computeOri a logical argument to determine the special case when the 
#' function was used as a sub function for the function computeOriPP
#' 
#' @return none
#' @export
computeReport <- function(genome, fasAnno, coreSet, extend=FALSE, 
                          scoreMode, priorityList=NULL, 
                          cpu, computeOri=FALSE) {
  if (!endsWith(coreSet, "/")) {
    coreSet <- paste(coreSet, "/", sep="");
  }
  modeFolder <- paste(coreSet, "phyloprofile", "/", 
                      as.character(scoreMode), sep="");
  if (!dir.exists(modeFolder)) {
    dir.create(modeFolder, recursive=TRUE);
  }
  
  splited <- strsplit(genome, "/", fixed=TRUE)[[1]];
  splited <- splited[length(splited)];
  genomeName <- strsplit(splited, ".", fixed=TRUE)[[1]][1];
  
  
  placeSeed(genome, fasAnno, coreSet, computeOri);
  
  if (is.null(priorityList)) {
    query <- list.dirs(paste(coreSet, "genome_dir", sep=""), 
                       recursive=FALSE, 
                       full.names=FALSE)[1];
    refSet <- list.dirs(paste(coreSet, "blast_dir", sep=""),
                        recursive=FALSE,
                        full.names=FALSE);
    priorityList <- autofindPriority(query, refSet);
  }
  
  if (extend == TRUE) {
    printPriority(genomeName, priorityList, coreSet, scoreMode);
  }
  
  pp <- runHamstr(coreSet, extend, scoreMode, priorityList, cpu);
  
  report <- reportSingle(pp, coreSet, scoreMode, priorityList);
  if (extend == TRUE) {
    translated <- translateReport(genomeName, report, scoreMode);
    printReport(translated, coreSet, scoreMode);
  }
  
  unlink(paste(coreSet, "genome_dir", "/", genomeName, sep=""), recursive=TRUE);
  if (computeOri == FALSE) {
    file.remove(paste(coreSet, "weight_dir", "/", genomeName, ".json", sep=""));
  }
  return(report);
}