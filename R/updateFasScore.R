#' The function correct information in a original domain file. This function is 
#' just applied when scoreMode=1
#' 
#' @param domain The path to the domain file
#' 
#' @return none
#' @export 
extendDomain <- function(domain) {
  domainFile <- read.table(paste(domain, sep=""), sep="\t", comment.char="");
  compareSet <- unique(as.vector(domainFile$V1));
  domainSet <- lapply(compareSet, 
                      function(compare, domainFile) {
                        sub <- subset(domainFile, V1 == compare);
                        splited <- strsplit(compare, "#", fixed=TRUE)[[1]];
                        splited2 <- strsplit(splited[1], "|", fixed=TRUE)[[1]];
                        sub$V1 <- paste(splited2[1], "#", splited[2], sep="");
                        return(sub);
                      },
                      domainFile=domainFile);
  domainFile <- do.call("rbind", domainSet);
  write.table(domainFile, domain, row.names=FALSE,
              sep="\t", quote=FALSE, col.names=FALSE);
}             

#' The function calculate the FAS score in the scoreMode 1 and update the scores
#' into phylogenetic profile of a specific core gene
#' 
#' @param coreSet the path to the core set
#' @param coreGene the ID of the core gene
#' @param extend the logical argument to decide if the information of the
#' interested genome will be appended to the original files
#' @param refSpec the references species
#' 
#' @return none
#' @export
getFasScore<- function(coreSet, coreGene, extend=FALSE, refSpec) {
  if (!endsWith(coreSet, "/")) {
    coreSet <- paste(coreSet, "/", sep="");
  }
  splited <- strsplit(coreSet, "/", fixed=TRUE)[[1]];
  setName <- splited[length(splited)];
  
  seedFasta <- paste(coreSet, "core_orthologs", "/", coreGene, "/", 
                          coreGene, ".fa", sep="");
  queryFasta <- paste(coreSet, "phyloprofile", "/", "1", "/", "hamstrout", "/", 
                      coreGene, "/", coreGene, ".extended.fa", sep="");
  
  seedFastaLines <- readLines(seedFasta);
  i <- 1:length(seedFastaLines);
  seedSet <- lapply(seedFastaLines[i[i %% 2 == 1]], 
                    function(seedID) {
                      return(substr(seedID, 2, nchar(seedID)));
                    });
  seedSet <- unlist(seedSet);
  refID <- NULL;
  for (seed in seedSet) {
    compare <- strsplit(seed, "|", fixed=TRUE)[[1]][2];
    if (compare == refSpec) {
      refID <- seed;
      break;
    }
  }
  
  annoDir <- paste(coreSet, "core_orthologs", "/", coreGene, "/", "fas_dir", "/", 
                   "annotation_dir", sep="");
  
  queryPP <- read.table(paste(coreSet, "phyloprofile", "/", "1", "/","hamstrout", "/", 
                              coreGene, "/", coreGene, ".phyloprofile", sep=""), header=TRUE, sep="\t");
  orthologNumber <- length(queryPP$orthoID);
  queryGenome <- strsplit(as.character(queryPP$orthoID[1]), "|", fixed=TRUE)[[1]][2];
  
  jobname <- paste(coreSet, "phyloprofile", "/", "1", "/", "hamstrout", "/",
                   coreGene, sep="");
  
  R.utils::createLink(paste(annoDir, "/", queryGenome, ".json", sep=""),
                      paste(coreSet, "weight_dir", "/", queryGenome, ".json", sep=""),
                      overwrite=TRUE);
  #calculate FAS scores 
  fasScores <- lapply(seedSet,
                      function(seedID, seedFasta, queryFasta, annoDir,
                               coreSet, setName, refID, jobname, 
                               queryGenome, coreGene) {
                        splited <- strsplit(seedID, "|", fixed=TRUE)[[1]];
                        R.utils::createLink(paste(annoDir, "/", splited[2], 
                                                  ".json", sep=""),
                                            paste(coreSet, "weight_dir", "/", 
                                                  splited[2], ".json", sep=""),
                                            overwrite=TRUE);
                        command <- paste("calcFAS",
                                         "-q", queryFasta,
                                         "-s", seedFasta,
                                         "-o", jobname,
                                         "-a", annoDir,
                                         "--seed_id", paste('"', seedID, '"', sep=""),
                                         "-t", 10,
                                         "--ref_2", paste(coreSet, "blast_dir", 
                                                          "/", splited[2], "/", 
                                                          splited[2], ".fa", sep=""),
                                         "-r", paste(coreSet, "genome_dir",
                                                     "/", queryGenome, "/", 
                                                     queryGenome, ".fa", sep=""),
                                         "--bidirectional",
                                         "-n", coreGene,
                                         "--raw");
                        if (seedID != refID) {
                          command <- paste(command, "--tsv");
                        }
                        lines <- system(command, intern=TRUE);
                        file.remove(paste(annoDir, "/", splited[2], 
                                           ".json", sep=""));
                        scores <- c();
                        for (line in lines) {
                          if (startsWith(line, "#")) {
                            splited <- strsplit(line, "\t", fixed=TRUE)[[1]];
                            score <- splited[length(splited)];
                            scores <- c(scores, score)
                          }
                        }
                        return(as.numeric(scores));
                      },
                      queryFasta=queryFasta,
                      seedFasta=seedFasta,
                      coreSet=coreSet,
                      setName=setName,
                      refID=refID,
                      jobname=jobname,
                      queryGenome=queryGenome,
                      annoDir=annoDir,
                      coreGene=coreGene);
  file.remove(paste(annoDir, "/", queryGenome, ".json", sep=""));
  if (extend == TRUE) {
    try(extendDomain(paste(jobname, "/", coreGene, "_forward.domains", sep="")),
        silent=TRUE);
    try(extendDomain(paste(jobname, "/", coreGene, "_reverse.domains", sep="")),
        silent=TRUE);
  }
  i <- 1:orthologNumber;
  FAS_F <- lapply(i, 
                  function(i, fasScores) {
                    scores <- lapply(fasScores, 
                                     function(scoreSet, i) {
                                       return(scoreSet[i]);
                                     }, i=i);
                    return(mean(as.numeric(scores)));
                  },
                  fasScores=fasScores);
  i <- (orthologNumber + 1):(orthologNumber * 2);
  FAS_B <- lapply(i, 
                  function(i, fasScores) {
                    scores <- lapply(fasScores, 
                                     function(scoreSet, i) {
                                       return(scoreSet[i]);
                                     }, i=i);
                    return(mean(as.numeric(scores)));
                  },
                  fasScores=fasScores);
  FAS_F <- unlist(FAS_F);
  FAS_B <- unlist(FAS_B);
  pp <- cbind(queryPP, FAS_F, FAS_B);
  write.table(pp, paste(coreSet, "phyloprofile", "/", "1", "/", 
                    "hamstrout", "/", coreGene, "/", coreGene,
                    ".phyloprofile", sep=""),
              sep="\t",
              row.names=FALSE,
              quote=FALSE);
}

#' The function calculate the FAS score in the scoreMode 1 and update the scores
#' into phylogenetic profile of a specific core gene
#' 
#' @param coreSet the path to the core set
#' @param coreGene the ID of the core gene
#' @param extend the logical argument to decide if the information of the
#' interested genome will be appended to the original files
#' @param refSpec the references species
#' 
#' @return none
#' @export
updateFasScore <- function(coreSet, coreGene, extend, refSpec) {
  if (!endsWith(coreSet, "/")) {
    coreSet <- paste(coreSet, "/", sep="");
  }
  getFasScore(coreSet, coreGene, extend, refSpec);
  file.remove(paste(coreSet, "core_orthologs", "/", coreGene, "/", 
                    "fas_dir", "/", "annotation_dir", "/", coreGene, 
                    "extended.json", sep=""));
  unlink(paste(coreSet, "core_orthologs", "/", coreGene, "/", 
               "fas_dir", "/", "annotation_dir", "/", "tmp",sep=""),
         recursive=TRUE);
}