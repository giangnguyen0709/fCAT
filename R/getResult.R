#' The function extract a phylogenetic profile of a genome from a phylogentic 
#' profile
#' 
#' @param pp the data frame of the phylogenetic profile
#' @param genome the ID of the extracted genome
#' 
#' @return the phylogenetic profile of the interested genome
#' @export
extractPP <- function(pp, genome) {
  genomeID <- unlist(lapply(pp$orthoID, 
                            function(orthoID) {
                              return(strsplit(orthoID, "|", fixed=TRUE)[[1]][2]);
                            }));
  singlePP <- cbind(pp, genomeID);
  singlePP <- subset(singlePP, genomeID == genome);
  return(singlePP[1:(ncol(singlePP)-1)]);
}

#' The function to assess the status of a founded ortholog
#' 
#' @param fasF the forward fas score of the ortholog
#' @param fasB the backward fas score of the ortholog
#' @param coreSet the path to the core set
#' @param coreGene the ID of the core gene
#' @param scoreMode the mode to determines the method to assess the ortholog
#' @param f the frequent of the core genes in the pp
#' @param priorityList the priority list to determines the references species
#' 
#' @return the status of the ortholog of the core gene
#' @export
assessStatus <- function(fasF, fasB, coreSet, coreGene, 
                         scoreMode, f, priorityList) {
  if (scoreMode == 1) {
    fas <- (fasF + fasB) / 2;
    cutoff <- read.table(paste(coreSet, "core_orthologs", "/", coreGene, "/",
                               "fas_dir", "/", "score_dir", "/", 
                               "1.cutoff", sep=""),
                         header=TRUE,
                         sep="\t");
    cutoffValue <- cutoff[1, 2];
  }
  
  if (scoreMode == 2) {
    fas <- (fasF + fasB) / 2;
    refSpec <- getSpec(paste(coreSet, "core_orthologs", "/", coreGene, "/", 
                             coreGene, ".fa", sep=""),
                       priorityList);
    cutoff <- read.table(paste(coreSet, "core_orthologs", "/", coreGene, "/",
                               "fas_dir", "/", "score_dir", "/", 
                               "2.cutoff", sep=""),
                         header=TRUE,
                         sep="\t");
    subCutoff <- subset(cutoff, taxa == refSpec);
    cutoffValue <- subCutoff[1, 2];
  }
  
  if (scoreMode == 3) {
    fas <- (fasF + fasB) / 2;
    refSpec <- getSpec(paste(coreSet, "core_orthologs", "/", coreGene, "/", 
                             coreGene, ".fa", sep=""),
                       priorityList);
    meanTable <- read.table(paste(coreSet, "core_orthologs", "/", coreGene, "/",
                                  "fas_dir", "/", "score_dir", "/",
                                  "1.cutoff", sep=""),
                            header=TRUE,
                            sep="\t");
    lcl <- meanTable[2, 2];
    ucl <- meanTable[3, 2];
    if (fas < lcl || fas > ucl ) {
      status <- "dissimilar";
    } else {
      status <- "similar";
    }
  }
  
  if (scoreMode != 3) {
    if (fas < cutoffValue) {
      status <- "dissimilar";
    } else {
      status <- "similar";
    }
  }
  
  if (f >= 2) {
    status <- paste("duplicated", ",", status);
  }
  
  return(status);
}

#' The function to calculate the cutoff value for the busco mode
#' 
#' @param coreSet the path to the core set
#' @param coreGene the ID of the core gene
#' 
#' @return a list that contains the mean length and the standard deviation of 
#' the length of the core gene
#' @export
calculateBuscoCutoff <- function(coreSet, coreGene) {
  if (!endsWith(coreSet, "/")) {
    coreSet <- paste(coreSet, "/", sep="");
  }
  
  seedFasta <- readLines(paste(coreSet, "core_orthologs", "/", coreGene, "/",
                               coreGene, ".fa", sep=""));
  i <- 1:length(seedFasta);
  seedFasta <- seedFasta[i[i %% 2 == 0]];
  lengthSet <- lapply(seedFasta, 
                      function(seq) {
                        return(nchar(seq));
                      });
  lengthSet <- unlist(lengthSet);
  
  meanLength <- mean(lengthSet);
  standardDeviation <- sqrt(mean((lengthSet - meanLength) ** 2));
  
  return(list(meanLength, standardDeviation));
}

#' The function to assess the founded ortholog with the algorithm of busco
#' 
#' @param orthoLength the length of the ortholg sequence
#' @param coreSet the path to the core set
#' @param coreGene the ID of the core gene
#' @param f the frequent of the core gene in the pp
#' 
#' @return the status of the core gene
#' @export
assessBusco <- function(orthoLength, coreSet, coreGene, f) {
  cutoff <- calculateBuscoCutoff(coreSet, coreGene);
  if (cutoff[[2]] != 0) {
    score <- (orthoLength - cutoff[[1]]) / cutoff[[2]];
    if (score > 2 || score < (-2)) {
      status <- "fragmented";
    } else {
      status <- "complete";
    }
  } else {
    score <- orthoLength - cutoff[[1]];
    if (score == 0) {
      status <- "complete";
    } else {
      status <- "fragmented";
    }
  }
  if (f >= 2) {
    status <- paste("duplicated", ",", status);
  }
  
  return(list(status, cutoff[[1]], cutoff[[2]]));
}

#' Determine if a core gene was ignored by the tool because of the unknown 
#' references species
#' 
#' @param coreSet the path to the core set
#' @param coreGene the ID of the core gene
#' @param priorityList the priority list to determine the references species
#' 
#' @return TRUE or FALSE
#' @export
filterIgnore <- function(coreSet, coreGene, priorityList) {
  if (!endsWith(coreSet, "/")) {
    coreSet <- paste(coreSet, "/", sep="");
  }
  fasta <- paste(coreSet, "core_orthologs", "/", coreGene, "/", 
                 coreGene, ".fa", sep="");
  check <- getSpec(fasta, priorityList);
  if (is.null(check)) {
    return(TRUE);
  } else {
    return(FALSE);
  }
}

#' Create the report of the completeness of a genome based on its phylogenetic 
#' profile
#' 
#' @param pp the phylogenetic profile of the genome in data frame
#' @param coreSet the path to the core set
#' @param scoreMode the mode to determined the method to assess the ortholog
#' @param priorityList the list to determinde the references species
#' 
#' @return the report in data frame
#' @export
reportSingle <- function(pp, coreSet, scoreMode, priorityList) {
  if (!endsWith(coreSet, "/")) {
    coreSet <- paste(coreSet, "/", sep="");
  }
  
  coreGeneList <- list.dirs(paste(coreSet, "core_orthologs", sep=""),
                            recursive=FALSE,
                            full.names=FALSE);
  frequency <- table(pp$geneID);
  
  if (scoreMode != "busco") {
    status <- unlist(lapply(1:nrow(pp), 
                            function(i, frequency, pp, scoreMode, 
                                     coreSet, priorityList) {
                              fasF <- pp[i, 4];
                              fasB <- pp[i, 5];
                              coreGene <- pp[i, 1];
                              f <- frequency[coreGene];
                              s <- assessStatus(fasF, fasB, coreSet, coreGene, 
                                                scoreMode, f, priorityList);
                              return(s);
                            },
                            frequency=frequency,
                            pp=pp,
                            scoreMode=scoreMode,
                            coreSet=coreSet,
                            priorityList=priorityList));
  } else {
    status <- lapply(1:nrow(pp), 
                     function(i, frequency, pp, coreSet) {
                       orthoLength <- pp[i, 4];
                       coreGene <- pp[i, 1];
                       f <- frequency[coreGene];
                       info <- assessBusco(orthoLength, coreSet, 
                                           coreGene, f);
                       r <- data.frame(status=c(info[[1]]),
                                       mean_length=c(info[[2]]),
                                       standard_deviation=c(info[[3]]));
                       return(r);
                     },
                     frequency=frequency,
                     pp=pp,
                     coreSet=coreSet);
    status <- do.call("rbind", status);
  }
  missingGene <- setdiff(coreGeneList, unique(pp$geneID));
  if (length(missingGene) != 0) {
    missingStatus <- unlist(lapply(missingGene, 
                                   function(gene, coreSet, priorityList) {
                                     check <- filterIgnore(coreSet, gene, 
                                                           priorityList);
                                     if (check == TRUE) {
                                       return("ignored");
                                     } else {
                                       return("missing")
                                     }
                                   },
                                   coreSet=coreSet,
                                   priorityList=priorityList));
    if (scoreMode != "busco") {
      report <- data.frame(geneID=pp$geneID, orthoID=pp$orthoID, status, 
                           FAS_F=pp$FAS_F, FAS_B=pp$FAS_B);
      missingTable <- data.frame(geneID=missingGene, orthoID=NA, 
                                 status=missingStatus, FAS_F=NA, FAS_B=NA);
      report <- rbind(report, missingTable);
    } else {
      report <- data.frame(geneID=pp$geneID, orthoID=pp$orthoID, 
                           status=status$status, length=pp$length, 
                           mean_length=status$mean_length, 
                           standard_deviation=status$standard_deviation);
      missingTable <- data.frame(geneID=missingGene, orthoID=NA, 
                                 status=missingStatus, length=NA,
                                 mean_length=NA,
                                 standard_deviation=NA);
      report <- rbind(report, missingTable);
    }
  }
  return(report);
}

#' Translate the report table into a frequent table, how many is complete, etc.
#' 
#' @param genomeID the genome ID of the interested genome
#' @param report the report in data frame
#' @scoreMode the mode to determine the method to assess the ortholog
#' 
#' @return A frequencies table in data frame
#' @export
translateReport <- function(genomeID, report, scoreMode) {
  if (scoreMode == "busco") {
    frequency <- table(report$status);
    
    complete <- frequency["complete"];
    fragmented <- frequency["fragmented"];
    missing <- frequency["missing"];
    ignored <- frequency["ignored"];
    
    if (is.na(complete)) {
      complete <- 0;
    }
    if (is.na(fragmented)) {
      fragmented <- 0;
    }
    if (is.na(missing)) {
      missing <- 0;
    }
    if (is.na(ignored)) {
      ignored <- 0;
    }
    
    duplicated <- length(unique(report$geneID)) - complete - 
      fragmented - missing - ignored;
    translated <- data.frame(genomeID=c(genomeID),
                             complete=c(complete),
                             fragmented=c(fragmented),
                             missing=c(missing),
                             duplicated=c(duplicated),
                             ignored=c(ignored));
    return(translated);
  } else {
    frequency <- table(report$status);
    
    similar <- frequency["similar"];
    dissimilar <- frequency["dissimilar"];
    missing <- frequency["missing"];
    ignored <- frequency["ignored"];
    if (is.na(similar)) {
      similar <- 0;
    }
    if (is.na(dissimilar)) {
      dissimilar <- 0;
    }
    if (is.na(missing)) {
      missing <- 0;
    }
    if (is.na(ignored)) {
      ignored <- 0;
    }
    
    duplicated <- length(unique(report$geneID)) - similar - 
      dissimilar - missing - ignored;
    translated <- data.frame(genomeID=c(genomeID),
                             similar=c(similar),
                             dissimilar=c(dissimilar),
                             missing=c(missing),
                             duplicated=c(duplicated),
                             ignored=c(ignored));
    return(translated);
  }
}