#' After HaMStR searched ortholog for the interested genome, it will store the 
#' output in a temporary folder in the core set. The function will concanate all
#' the output files in the folder in a single files
#' 
#' @param directory path to the folder, that contains the output files
#' @param genomeName the genomeID of the genome, that need to be extracted
#' @return none
#' @export
concanateFiles <- function(directory, genomeName) {
  if (!endsWith(directory, "/")) {
    directory <- paste(directory, "/", sep="");
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
  pp <- extractPP(pp, genomeName);
  domain0 <- extractDomains(domain0, genomeName);
  domain1 <- extractDomains(domain1, genomeName);
  exFasta <- extractFasta(exFasta, genomeName);
  return(list(pp, exFasta, domain0, domain1));
}

#' Function to append the phylogenetic profile of the interested genome into the
#' original pp
#' 
#' @param root the path to the root folder
#' @param coreSet the core set name
#' @param scoreMode the mode to assess the founded ortholog
#' @param fileList a list that contains the information of the phylogenetic 
#' profile of the interested genome
#' @export
extendOriginal <- function(root, coreSet, scoreMode, fileList) {
  oriPath <- paste(root, "phyloprofile", "/", coreSet, "/", scoreMode, "/",
                   coreSet, sep="");
  if (file.exists(paste(oriPath, ".phyloprofile", sep=""))) {
    pp <- read.table(paste(oriPath, ".phyloprofile", sep=""),
                     header=TRUE,
                     sep="\t");
    pp <- rbind(pp, fileList[[1]]);
  } else {
    pp <- fileList[[1]];
  }
  
  if (file.exists(paste(oriPath, ".extended.fa", sep=""))) {
    exFasta <- readLines(paste(oriPath, ".extended.fa", sep=""));
    exFasta <- c(exFasta, fileList[[2]])
  } else {
    exFasta <- fileList[[2]];
  }
  
  if (file.exists(paste(oriPath, "_reverse.domains", sep=""))) {
    domain0 <- read.table(paste(oriPath, "_reverse.domains", sep=""),
                          sep="\t",
                          comment.char="");
    domain0 <- rbind(domain0, fileList[[3]]);
  } else {
    domain0 <- fileList[[3]];
  }
  
  if (file.exists(paste(oriPath, "_forward.domains", sep=""))) {
    domain1 <- read.table(paste(oriPath, "_forward.domains", sep=""),
                          sep="\t",
                          comment.char="");
    domain1 <- rbind(domain1, fileList[[4]]);
  } else {
    domain1 <- fileList[[4]];
  }
  write.table(pp, 
              paste(oriPath, ".phyloprofile", sep=""),
              row.names=FALSE,
              sep="\t",
              quote=FALSE);
  write.table(domain0, 
              paste(oriPath, "_reverse.domains", sep=""),
              row.names=FALSE,
              col.names=FALSE,
              sep="\t",
              quote=FALSE);
  write.table(domain1, 
              paste(oriPath, "_forward.domains", sep=""),
              row.names=FALSE,
              col.names=FALSE,
              sep="\t",
              quote=FALSE);
  writeLines(exFasta, paste(oriPath, ".extended.fa", sep=""));
}

#' This function takes a path to a core set and run HaMStR to search ortholog on
#' the genome in the folder genome_dir of core set
#' 
#' @param root Path to the root folder
#' @param coreSet The core set name
#' @param extend if extend=TRUE the phylogenetic profile of the genome will be 
#' appended to the original phylogenetic profile
#' @param scoreMode the mode determines the way to assess the founded ortholog
#' @param priorityList the list determines the references species
#' @param cpu determines the cores that HaMStR will use
#' 
#' @return phylogenetic profile of the genome
#' @export
runHamstr <- function(root, coreSet, extend=FALSE, 
                      scoreMode, priorityList=NULL, cpu){
  if (!endsWith(root, "/")) {
    root <- paste(root, "/", sep="");
  }
  
  setName <- coreSet;
  
  genomeName <- list.dirs(paste(root, "check_dir", sep=""), 
                          recursive=FALSE,
                          full.names=FALSE)[[1]];
  
  hmmPath <- paste(root, "core_orthologs", "/", coreSet, sep="");
  blastPath <- paste(root, "blast_dir", sep="");
  searchPath <- paste(root, "check_dir", sep="");
  weightPath <- paste(root, "weight_dir", sep="");
  outPath <- paste(root, "phyloprofile", "/", coreSet, "/", 
                   as.character(scoreMode), "/", "hamstrout", sep="");
  ### - check Data - ###
  command <- paste("checkData1s",
                   "-g", searchPath,
                   "-b", blastPath,
                   "-w", weightPath);
  system(command);
  if (!dir.exists(outPath)) {
    dir.create(outPath, recursive=TRUE);
  }
  ppSet <- lapply(list.dirs(paste(root, "core_orthologs", "/", 
                                  coreSet, sep=""), 
                            recursive=FALSE, full.names=FALSE),
                  function(coreGene, hmmPath, blastPath, searchPath, 
                           outPath, weightPath, root, coreSet, scoreMode, 
                           extend, priorityList, genomeName) {
                    refSpec <- getSpec(paste(hmmPath, "/", coreGene,
                                             "/", coreGene, ".fa", 
                                             sep = ""), 
                                       priorityList);
                    if (is.null(refSpec)) {
                      return(NULL);
                    }
                    command <- paste("h1s",
                                     "--seqFile", paste(root,
                                                        "core_orthologs",
                                                        "/",
                                                        coreSet,
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
                    
                    if (!file.exists(paste(root, "phyloprofile",
                                           "/", coreSet,
                                           "/", as.character(scoreMode), 
                                           "/", "hamstrout", "/",
                                           coreGene, "/", coreGene, 
                                           ".phyloprofile", sep=""))) {
                      return(NULL);
                    } else {
                      if (scoreMode == 1) {
                        updateFasScore(root, coreSet, coreGene, 
                                       extend, refSpec);
                        pp <- read.table(paste(root, "phyloprofile",
                                               "/", coreSet,
                                               "/", as.character(scoreMode), 
                                               "/", "hamstrout", "/",
                                               coreGene, "/", coreGene, 
                                               ".phyloprofile", sep=""),
                                         header=TRUE,
                                         sep="\t");
                        return(pp); 
                      } 
                      if (scoreMode == 2 || scoreMode == 3) {
                        pp <- read.table(paste(root, "phyloprofile",
                                               "/", coreSet,
                                               "/", as.character(scoreMode), 
                                               "/", "hamstrout", "/",
                                               coreGene, "/", coreGene, 
                                               ".phyloprofile", sep=""),
                                         header=TRUE,
                                         sep="\t");
                        return(pp);
                      }
                      if (scoreMode == "busco") {
                        updateLength(root, coreSet, coreGene);
                        pp <- read.table(paste(root, "phyloprofile",
                                               "/", coreSet,
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
                  root=root,
                  coreSet=coreSet,
                  scoreMode=scoreMode,
                  extend=extend,
                  priorityList);
  pp <- do.call("rbind", ppSet);
  pp <- extractPP(pp, genomeName);
  if (extend == TRUE) {
    outFolder <- paste(root, "phyloprofile", "/", coreSet, "/", 
                       as.character(scoreMode), sep="")
    fileList <- concanateFiles(outFolder, genomeName);
    extendOriginal(root, coreSet, scoreMode, fileList);
  }
  unlink(outPath, recursive=TRUE);
  return(pp);
}