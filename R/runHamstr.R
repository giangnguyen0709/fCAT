#' After HaMStR searched ortholog for the interested genome, it will store the 
#' output in a temporary folder in the core set. The function will concanate all
#' the output files in the folder in a single files
#' 
#' @param directory path to the folder, that contains the output files
#' @param out the output folder, that will contain the concanated files
#' @param name the concanated files will be saved under this name
#' @return none
#' @export
concanateFiles <- function(directory, out, name) {
  if (!endsWith(directory, "/")) {
    directory <- paste(directory, "/", sep="");
  }
  if (!endsWith(out, "/")) {
    out <- paste(out, "/", sep="")
  }
  
  exFasta <- NULL;
  domain0 <- NULL;
  domain1 <- NULL;
  pp <- NULL;
  
  for (file in list.files(directory, full.names=TRUE, recursive=TRUE)) {
    if (endsWith(file, ".extended.fa")) {
      if (is.null(exFasta)) {
        exFasta <- readLines(file);
      } else {
        exFasta <- c(exFasta, readLines(file));
      }
    }
    
    if (endsWith(file, "_reverse.domains")) {
      domain <- try(read.table(file, sep="\t", comment.char=""), silent=TRUE);
      if (!inherits(domain, "try-error")) {      
        if (is.null(domain0)) {
          domain0 <- domain;
          } else {
          domain0 <- rbind(domain0, domain);
          } 
        }
      }
    
    if (endsWith(file, "_forward.domains")) {
      domain <- try(read.table(file, sep="\t", comment.char=""), silent=TRUE);
      if (!inherits(domain, "try-error")) {
        if (is.null(domain1)) {
          domain1 <- read.table(file, sep="\t", comment.char="");
        } else {
          domain1 <- rbind(domain1, read.table(file, sep="\t", comment.char=""));
        }
      }
    }
    
    if (endsWith(file, ".phyloprofile")) {
      if (is.null(pp)) {
        pp <- read.table(file, sep="\t", header=TRUE);
      } else {
        pp <- rbind(pp, read.table(file, sep="\t", header=TRUE));
      }
    }
  }
  jobname <- paste(out, name, sep="")
  write.table(pp, 
              paste(jobname, ".phyloprofile", sep=""),
              row.names=FALSE,
              sep="\t",
              quote=FALSE);
  write.table(domain0, 
              paste(jobname, "_reverse.domains", sep=""),
              row.names=FALSE,
              col.names=FALSE,
              sep="\t",
              quote=FALSE);
  write.table(domain1, 
              paste(jobname, "_forward.domains", sep=""),
              row.names=FALSE,
              col.names=FALSE,
              sep="\t",
              quote=FALSE);
  writeLines(exFasta, paste(jobname, ".extended.fa", sep=""));
}

#' This function takes a path to a core set and run HaMStR to search ortholog on
#' the genome in the folder genome_dir of core set
#' 
#' @param coreSet the path to the core set
#' @param extend if extend=TRUE the phylogenetic profile of the genome will be 
#' appended to the original phylogenetic profile
#' @param scoreMode the mode determines the way to assess the founded ortholog
#' @param priorityList the list determines the references species
#' @param cpu determines the cores that HaMStR will use
#' 
#' @return phylogenetic profile of the genome
#' @export
runHamstr <- function(coreSet, extend=FALSE, scoreMode, priorityList=NULL, cpu){
  start <- Sys.time();
  if (!endsWith(coreSet, "/")) {
    coreSet <- paste(coreSet, "/", sep="");
  }
  
  splited <- strsplit(coreSet, "/", fixed=TRUE)[[1]];
  setName <- splited[length(splited)];
  
  hmmPath <- paste(coreSet, "core_orthologs", sep="");
  blastPath <- paste(coreSet, "blast_dir", sep="");
  searchPath <- paste(coreSet, "genome_dir", sep="");
  weightPath <- paste(coreSet, "weight_dir", sep="");
  outPath <- paste(coreSet, "phyloprofile", "/", as.character(scoreMode), "/",
                   "hamstrout", sep="");
  ### - check Data - ###
  command <- paste("checkData1s",
                   "-g", searchPath,
                   "-b", blastPath,
                   "-w", weightPath);
  system(command);
  if (!dir.exists(outPath)) {
    dir.create(outPath, recursive=TRUE);
  }
  ppSet <- lapply(list.dirs(paste(coreSet, "core_orthologs", sep=""), 
                                        recursive=FALSE, full.names=FALSE),
                              function(coreGene, hmmPath, blastPath, searchPath, 
                                       outPath, weightPath, coreSet, scoreMode, 
                                       extend, priorityList) {
                                refSpec <- getSpec(paste(hmmPath, "/", coreGene,
                                                         "/", coreGene, ".fa", 
                                                         sep = ""), 
                                                   priorityList);
                                if (is.null(refSpec)) {
                                  return(NULL);
                                }
                                command <- paste("h1s",
                                                 "--seqFile", paste(coreSet,
                                                                    "core_orthologs",
                                                                    "/", coreGene,
                                                                    "/", coreGene,
                                                                    ".fa", sep=""),
                                                 "--seqName", coreGene,
                                                 "--refspec", refSpec,
                                                 "--hmmpath", hmmPath,
                                                 "--outpath", outPath,
                                                 "--blastpath", blastPath,
                                                 "--weightpath", weightPath,
                                                 "--searchpath", searchPath,
                                                 "--cleanup",
                                                 "--cpu", cpu,
                                                 "--reuseCore",
                                                 "--checkCoorthologsRef",
                                                 "--countercheck");
                                if (scoreMode == 1 || scoreMode == "busco") {
                                  command <- paste(command, "--fasoff");
                                }
                                system(command);
                                
                                if (file.exists(paste(coreGene, ".fa", sep=""))) {
                                  file.remove(paste(coreGene, ".fa", sep=""));
                                }
                                
                                if (!file.exists(paste(coreSet, "phyloprofile",
                                                       "/", as.character(scoreMode), 
                                                       "/", "hamstrout", "/",
                                                       coreGene, "/", coreGene, 
                                                       ".phyloprofile", sep=""))) {
                                  return(NULL);
                                } else {
                                  if (scoreMode == 1) {
                                    updateFasScore(coreSet, coreGene, 
                                                   extend, refSpec);
                                    pp <- read.table(paste(coreSet, "phyloprofile",
                                                           "/", as.character(scoreMode), 
                                                           "/", "hamstrout", "/",
                                                           coreGene, "/", coreGene, 
                                                           ".phyloprofile", sep=""),
                                                     header=TRUE,
                                                     sep="\t");
                                    return(pp); 
                                  } 
                                  if (scoreMode == 2 || scoreMode == 3) {
                                    pp <- read.table(paste(coreSet, "phyloprofile",
                                                           "/", as.character(scoreMode), 
                                                           "/", "hamstrout", "/",
                                                           coreGene, "/", coreGene, 
                                                           ".phyloprofile", sep=""),
                                                     header=TRUE,
                                                     sep="\t");
                                    return(pp);
                                  }
                                  if (scoreMode == "busco") {
                                    updateLength(coreSet, coreGene);
                                    pp <- read.table(paste(coreSet, "phyloprofile",
                                                           "/", as.character(scoreMode), 
                                                           "/", "hamstrout", "/",
                                                           coreGene, "/", coreGene, 
                                                           ".phyloprofile", sep=""),
                                                     header=TRUE,
                                                     sep="\t");
                                    return(pp);
                                  }
                                }
                              },
                              hmmPath=hmmPath,
                              blastPath=blastPath,
                              searchPath=searchPath, 
                              weightPath=weightPath,
                              outPath=outPath,
                              coreSet=coreSet,
                              scoreMode=scoreMode,
                              extend=extend,
                              priorityList);
  pp <- do.call("rbind", ppSet);
  if (extend == TRUE) {
    outFolder <- paste(coreSet, "phyloprofile", "/", 
                       as.character(scoreMode), sep="")
    concanateFiles(outFolder, outFolder, setName);
  }
  unlink(outPath, recursive=TRUE);
  return(pp);
}