#' The function creates FAS annotation for a core gene in the core set
#' 
#' @param coreSet The path to the core set
#' @param coreGene The name of the core gene in the set
#' @return none
#' @export
createAnnotation <- function(coreSet, coreGene) {
  if (!endsWith(coreSet, "/")) {
    coreSet <- paste(coreSet, "/", sep="");
  }
  
  annoFolder <- paste(coreSet, "core_orthologs", "/", coreGene, "/", "fas_dir", "/", 
                      "annotation_dir", sep="");
  if (!dir.exists(annoFolder)) {
    dir.create(annoFolder, recursive=TRUE);
  }
  
  command <- paste("annoFAS", 
                   "-i", paste(coreSet, "core_orthologs", "/", coreGene, "/", coreGene, ".fa", sep=""), 
                   "-o", annoFolder,
                   "-n", paste(coreGene, sep=""));
  system(command);
  tmp <- paste(annoFolder, "/", "tmp", sep="");
  if (dir.exists(tmp)) {
    unlink(tmp, recursive=TRUE);
  }
}

#' The function creates all FAS annotation for all core genes in the set
#' 
#' @param coreSet the path to the core set
#' 
#' @return none
#' @export
createAllAnnotation <- function(coreSet) {
  startTime <- Sys.time();
  if (!endsWith(coreSet, "/")) {
    coreSet <- paste(coreSet, "/", sep="");
  }
  lapply(list.dirs(paste(coreSet, "core_orthologs", sep=""), recursive=FALSE, full.names=FALSE), 
         function(coreGene, coreSet){
           print(paste("Annotating for", coreGene, sep=" "));
           createAnnotation(coreSet, coreGene);
         }, 
         coreSet=coreSet);
  endTime <- Sys.time();
  print(endTime - startTime);
}

#' The function calculate all values, that are necessary for the assessment 
#' process for a specific core gene in the set. The function will saved the 
#' values as the file text in the folder of the core gene
#' 
#' @param coreSet the path to the core set
#' @param coreGene the ID of the core gene in the set
#' 
#' @return none
#' @export
calculateCutoff <- function(coreSet, coreGene) {
  if (!endsWith(coreSet, "/")) {
    coreSet <- paste(coreSet, "/", sep="");
  }
  fastaFile <- paste(coreSet, "core_orthologs", "/", 
                     coreGene, "/", coreGene, ".fa", sep="")
  annoDir <- paste(coreSet, "core_orthologs", "/", coreGene, "/", 
                   "fas_dir", "/", "annotation_dir", sep="")
  
  fasta <- readLines(fastaFile);
  i <- 1:length(fasta);
  i <- i[i %% 2 == 1];
  querySet <- lapply(fasta[i],
                       function(header) {
                         return(substr(header, 2, nchar(header)));
                       });
  
  genomeSet <- lapply(querySet, function(query) {
    return(strsplit(query, "|", fixed=TRUE)[[1]][2]);
  })
  
  scoreDist <- list();
  genomeScores <- list();
  for (i in 1:length(genomeSet)) {
    genomeScores[[genomeSet[[i]]]] <- 0;
    if (i != length(genomeSet)) {
      scoreDist[[genomeSet[[i]]]] <- list();
    }
  }
  
  for (i in 1:(length(genomeSet) - 1)) {
    for (j in (i+1):length(genomeSet)) {
      scoreDist[[genomeSet[[i]]]][[genomeSet[[j]]]] <- 0;
    }
  }

  scoreSet <- lapply(querySet, 
                     function(queryID, fastaFile, annoDir, coreSet) {
                       refSpec <- strsplit(queryID, "|", fixed=TRUE)[[1]][2];
                       refProteome <- paste(coreSet, "weight_dir", "/", 
                                            refSpec, ".json", sep="");
                       R.utils::createLink(paste(annoDir, "/", refSpec, ".json",
                                                 sep=""), refProteome,
                                           overwrite=TRUE);
                       command <- paste("calcFAS",
                                        "-q", fastaFile,
                                        "-s", fastaFile,
                                        "--query_id", paste('"', queryID, 
                                                            '"', sep=""),
                                        "-a", annoDir,
                                        "-o", annoDir,
                                        "--tsv",
                                        "-r", refProteome,
                                        "-t", "10",
                                        "--raw");
                       lines <- system(command, intern=TRUE);
                       scores <- list();
                       for (line in lines) {
                         if (startsWith(line, "#")) {
                           splited <- strsplit(line, "\t", fixed=TRUE)[[1]];
                           q <- strsplit(splited[3], "|", fixed=TRUE)[[1]][2];
                           s <- strsplit(splited[2], "|", fixed=TRUE)[[1]][2];
                           score <- as.numeric(splited[length(splited)]);
                           pack <- list(q, s, score);
                           scores[[length(scores) + 1]] <- pack;
                         }
                       }
                       file.remove(paste(annoDir, "/", refSpec, 
                                         ".json", sep=""));
                       return(scores);
                     },
                     fastaFile,
                     annoDir,
                     coreSet);
  for (level1 in scoreSet) {
    for (score in level1) {
      if (score[[1]] != score[[2]]) {
        a <- genomeScores[[score[[1]]]];
        b <- genomeScores[[score[[2]]]];
        genomeScores[[score[[1]]]] <- a + score[[3]];
        genomeScores[[score[[2]]]] <- b + score[[3]];
        if (!is.null(scoreDist[[score[[1]]]][[score[[2]]]])) {
          v <- scoreDist[[score[[1]]]][[score[[2]]]];
          scoreDist[[score[[1]]]][[score[[2]]]] <- v + score[[3]];
        } else {
          v <- scoreDist[[score[[2]]]][[score[[1]]]];
          scoreDist[[score[[2]]]][[score[[1]]]] <- v + score[[3]];
        }
      }
    }
  }
  scoreDist <- unlist(scoreDist) / 2;
  genomeScores <- unlist(genomeScores) / ((length(genomeSet) - 1)* 2);
  
  avaMean <- mean(genomeScores);
  
  features <- EnvStats::eexp(scoreDist, ci=TRUE);
  lcl <- 1 / (features$interval$limits[[2]]);
  ucl <- 1 / (features$interval$limits[[1]]);
  
  scoreFolder <- paste(coreSet, "core_orthologs", "/", coreGene, "/", "fas_dir",
                       "/", "score_dir", sep="");
  if (!dir.exists(scoreFolder)) {
    dir.create(scoreFolder, recursive=TRUE);
  }
  # Table 1:
  label <- c("mean", "LCL", "UCL");
  value <- c(avaMean, lcl, ucl);
  cutoffTable <- data.frame(label, value);
  
  cutoffFile <- paste(scoreFolder, "/", "1", ".cutoff", sep="");
  write.table(cutoffTable, cutoffFile, sep="\t", quote=FALSE, row.names=FALSE);
  
  # Table 2:
  cutoffTable <- data.frame(taxa=unlist(genomeSet), 
                            cutoff=genomeScores);
  cutoffFile <- paste(scoreFolder, "/", "2", ".cutoff", sep="");
  write.table(cutoffTable, cutoffFile, sep="\t", quote=FALSE, row.names=FALSE);
}

#' The function calculate cut off values for all core genes in the set
#' 
#' @param coreSet the path to the core set
#' 
#' @return none
#' @export
calculateAllCutoff <- function(coreSet) {
  startTime <- Sys.time();
  if (!endsWith(coreSet, "/")) {
    coreSet <- paste(coreSet, "/", sep="");
  }
  
  coreOrtho <- paste(coreSet, "core_orthologs", sep="");
  lapply(list.dirs(coreOrtho, full.names=FALSE, recursive=FALSE), 
         function(coreGene, coreSet, mode){
           print(paste("Starting calculate cutoff for", coreGene, sep=" "));
           calculateCutoff(coreSet, coreGene);
         }, coreSet=coreSet);
  
  endTime <- Sys.time();
  print("Done after: ");
  print(endTime - startTime);
}

#' The function compute the FAS annotations and calculate the cut off values for
#' all core genes in the set
#' 
#' @param coreSet the path to the core set
#' 
#' @return none
#' @export
processCoreSet <- function(coreSet) {
  if (!endsWith(coreSet, "/")) {
    coreSet <- paste(coreSet, "/", sep="");
  }
  
  createAllAnnotation(coreSet)
  calculateAllCutoff(coreSet)
}