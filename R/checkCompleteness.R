#' Remove the data of a failed check
#'
#' @param root the path to the root folder
#'
#' @return none
#' @export
handleError <- function(root, weightDir = NULL) {
    genomeDir <- paste(root, "query_taxon", sep = "")
    if (is.null(weightDir)) {
        weightDir <- paste(root, "weight_dir", sep = "")
    } else {
        if (endsWith(weightDir, "/")) {
            weightDir <- substr(weightDir, 1, length(weightDir) - 1)
        }
    }

    if (
        length(
            list.dirs(genomeDir, full.names = FALSE, recursive = FALSE
            )
        ) != 0
    ) {
        for (
            directory in list.dirs(
                genomeDir, full.names = FALSE, recursive = FALSE
            )
        ) {
            unlink(paste(genomeDir, "/", directory, sep = ""), recursive = TRUE)
            if (
                file.exists(
                    paste(weightDir, "/", directory, ".json", sep = "")
                )
            ) {
                file.remove(paste(weightDir, "/", directory, ".json", sep = ""))
            }
        }
    }
}

#' Check if the core set was processed
#'
#' @param root the path to the root folder
#' @param coreSet the path to the core set
#'
#' @return TRUE or FALSE
#' @export
checkPreProcess <- function(root, coreSet) {
    check <- 0
    for (
        coreGene in list.dirs(
            paste(root, "core_orthologs", "/", coreSet, sep = ""), 
            recursive = FALSE, 
            full.names = TRUE
        )
    ) {
        fasDir <- paste(coreGene, "/", "fas_dir", sep = "")
        annoDir <- paste(fasDir, "/", "annotation_dir", sep = "")
        scoreDir <- paste(fasDir, "/", "score_dir", sep = "")
        if (!dir.exists(annoDir) || !dir.exists(scoreDir)) {
            check <- 1
            break
        }
    }
    if (check == 0) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}

#' Check if the interested genome exists already in the original pp
#'
#' @param genomeName the ID of the genome
#' @param root the path to the root folder
#' @param coreSet the core set name
#' @param scoreMode the mode determines the method to assess the ortholog
#'
#' @return TRUE or FALSE
#' @export
checkExist <- function(genomeName, root, coreSet, scoreMode, ppDir) {
    setName <- coreSet
    
    if (!is.null(ppDir)) {
        reportFile <- paste(
            ppDir, setName, ".report", sep = ""
        )
    } else {
        reportFile <- paste(
            root, "output", "/", coreSet, "/",
            as.character(scoreMode), "/", setName, ".report",
            sep = ""
        )
    }
    
    if (!file.exists(reportFile)) {
        return(FALSE)
    }
    report <- read.table(reportFile,
        header = TRUE,
        sep = "\t"
    )
    if (genomeName %in% report$genomeID) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}

#' Check if all the input are correct
#'
#' @param genome The path to the genome fasta file
#' @param fasAnno The path to the fas annotation file. Can be NULL
#' @param root the path to the root folder
#' @param coreSet The core set name
#' @param extend A logical value
#' @param redo A logical value
#' @param scoreMode The mode determines the method to assess the orthologs
#' @param priorityList The list of taxa in the core set to determine the
#' references species
#' @param cpu The number of the cores that HaMStR will use
#'
#' @return A list that contains a logical value and the message to the error
#' @export
checkArguments <- function(
    genome, fasAnno = NULL, coreDir, coreSet, extend = FALSE, redo = FALSE, 
    scoreMode, priorityList = NULL, cpu = 4, blastDir = NULL, weightDir = NULL,
    outDir = NULL, cleanup = FALSE, reFdog = FALSE, fdogDir = NULL,
    ppDir = NULL
) {
    if (!endsWith(root, "/")) {
        root <- paste(root, "/", sep = "")
    }
    check <- TRUE
    status <- NULL
    if (!file.exists(genome)) {
        check <- FALSE
        status <- "Genome fasta file doesn't exist"
        return(list(check, status))
    }

    if (!is.null(fasAnno)) {
        if (!file.exists(fasAnno)) {
            check <- FALSE
            status <- "FAS annotation file doesn't exist"
            return(list(check, status))
        }
    }

    if (!dir.exists(root)) {
        check <- FALSE
        status <- "The root folder doesn't exist"
        return(list(check, status))
    } else {
        if (is.null(blastDir)) {
            blastDir <- paste(root, "blast_dir", sep = "")
            if (!dir.exists(blastDir)) {
                check <- FALSE
                status <- "blast_dir is missing"
                return(list(check, status))
            }
            
            if (length(list.dirs(blastDir)) == 0) {
                check <- FALSE
                status <- "blast_dir is empty"
                return(list(check, status))
            }
        } else {
            if (!dir.exists(blastDir)) {
                status <- "blast_dir is missing"
                check <- FALSE
                return(list(check, status))
            }
            
            if (length(list.dirs(blastDir)) == 0) {
                check <- FALSE
                status <- "blast_dir is empty"
                return(list(check, status))
            }
        }
        
        if (is.null(weightDir)) {
            weightDir <- paste(root, "weight_dir", sep = "")
            if (!dir.exists(weightDir)) {
                check <- FALSE
                status <- "weight_dir is missing"
                return(list(check, status))
            }
            
            if (length(list.dirs(weightDir)) == 0) {
                check <- FALSE
                status <- "blast_dir is empty"
                return(list(check, status))
            }
        } else {
            if (!dir.exists(weightDir)) {
                check <- FALSE
                status <- "weight_dir is missing"
                return(list(check, status))
            }
            
            if (length(list.dirs(weightDir)) == 0) {
                check <- FALSE
                status <- "blast_dir is empty"
                return(list(check, status))
            }
        }
        
        orthoDir <- paste(root, "core_orthologs", "/", coreSet, sep = "")
        if (!dir.exists(orthoDir)) {
            check <- FALSE
            status <- "The core set does not exist"
            return(list(check, status))
        } else {
            if (length(list.dirs(orthoDir)) == 0) {
                check <- FALSE
                status <- "The core_orthologs is empty"
                return(list(check, status))
            }
        }
        
    }


    modeList <- list(1, 2, 3, "busco")
    if (!(scoreMode %in% modeList)) {
        check <- FALSE
        status <- "score mode is not available"
        return(list(check, status))
    }

    if (!is.null(priorityList)) {
        if (!is.vector(priorityList)) {
            check <- FALSE
            status <- "priority list must be a list"
            return(list(check, status))
        }
        coreTaxa <- list.dirs(
            paste(root, "blast_dir", sep = ""),
            recursive = FALSE,
            full.names = FALSE
        )
        for (taxa in priorityList) {
            if (!(taxa %in% coreTaxa)) {
                check <- FALSE
                status <- 
                    "A taxa in the priority list doesn't belong to the core set"
                return(list(check, status))
            }
        }
    }

    return(list(check, status))
}

#' The function to check completeness of a new given genome
#'
#' @param genome The path to the fasta file of the genome
#' @param fasAnno the path to the json file of the annotation of the genome. If
#' equal NULL the function will create one
#' @param root the path to the root folder
#' @param coreSet the core set name
#' @param extend the logical option to decide if the information of the
#' interested will be appended to the original files and reused for the next
#' check
#' @param redo the logical option to decide if the tool should recheck for a
#' genome ID that already exists in the original pp
#' @param scoreMode the mode determines the method to assess the founded
#' ortholog
#' @param priorityList the priority list to determine the references species
#' @param cpu the number of the cores
#'
#' @return a list that contains a report table of the completeness of the
#' interested genome and a table that contains the report of the interested
#' genome within the other taxa in the original pp
#' @export
checkCompleteness <- function(
    genome, fasAnno = NULL, coreDir, coreSet, extend = FALSE, redo = FALSE, 
    scoreMode, priorityList = NULL, cpu = 4, blastDir = NULL, weightDir = NULL,
    outDir = NULL, cleanup = FALSE, reFdog = FALSE, fdogDir = NULL,
    ppDir = NULL
) {
    start <- Sys.time()
    check <- checkArguments(
        genome, fasAnno, coreDir, coreSet, extend, redo, scoreMode,
        priorityList, cpu, blastDir, weightDir, outDir, cleanup
    )
    if (check[[1]] == FALSE) {
        return(check[[2]])
    }

    if (!endsWith(coreDir, "/")) {
        coreDir <- paste(coreDir, "/", sep = "")
    }
    if (!is.null(weightDir)) {
        if (!endsWith(weightDir, "/")) {
            weightDir <- paste(weightDir, "/", sep = "")
        }   
    }
    if (!is.null(blastDir)) {
        if (!endsWith(blastDir, "/")) {
            blastDir <- paste(blastDir, "/", sep = "")
        }
    }
    if (!is.null(outDir)) {
        if (!endsWith(outDir, "/")) {
            outDir <- paste(outDir, "/", sep = "")
        }
    }
    if (!is.null(fdogDir)) {
        if (!endsWith(fdogDir, "/")) {
            fdogDir <- paste(fdogDir, "/", sep = "")
        }
    }
    if (!is.null(ppDir)) {
        if (!endsWith(ppDir, "/")) {
            ppDir <- paste(ppDir, "/", sep = "")
        }
    }
    
    if (!checkPreProcess(coreDir, coreSet)) {
        processCoreSet(coreDir, coreSet)
    }

    handleError(coreDir)

    setName <- coreSet

    splited <- strsplit(genome, "/", fixed = TRUE)[[1]]
    splited <- splited[length(splited)]
    genomeName <- strsplit(splited, ".", fixed = TRUE)[[1]][1]

    if (!checkExist(genomeName, coreDir, coreSet, scoreMode, ppDir)) {
        compute <- TRUE
    } else {
        if (redo == FALSE) {
            compute <- FALSE
        } else {
            if (extend == TRUE) {
                if (!is.null(ppDir)) {
                    correctFiles(ppDir, genomeName)
                } else {
                    correctFiles(
                        paste(
                            coreDir, "output", "/", coreSet, "/",
                            as.character(scoreMode),
                            sep = ""
                        ),
                        genomeName
                    )
                }
            }
            compute <- TRUE
        }
    }

    if (compute == TRUE) {
        singleReport <- computeReport(
            genome, fasAnno, coreDir, coreSet, extend,
            scoreMode, priorityList, cpu, FALSE, 
            blastDir, weightDir, cleanup, reFdog, fdogDir, ppDir
        )
        translated <- translateReport(genomeName, singleReport, scoreMode)
        
        if (!is.null(ppDir)) {
            reportFile <- paste(ppDir, setName, ".report", sep = "")
        } else {
            reportFile <- paste(
                coreDir, "output", "/", coreSet, "/",
                as.character(scoreMode), "/", setName, ".report",
                sep = ""
            )
        }
        if (file.exists(reportFile)) {
            allReport <- read.table(
                reportFile,
                header = TRUE,
                sep = "\t"
            )
            allReport <- subset(allReport, genomeID != genomeName)
            allReport <- rbind(translated, allReport)
        } else {
            allReport <- translated
        }
    } else {
        if (!is.null(ppDir)) {
            ppPath <- paste(ppDir, setName, ".phyloprofile", sep = "")
        } else {
            ppPath <- paste(
                coreDir, "output", "/", coreSet, "/", 
                as.character(scoreMode), "/", setName, 
                ".phyloprofile", sep = ""
            )
        }
        phyloprofile <- read.table(
            ppPath,
            header = TRUE,
            sep = "\t"
        )
        phyloprofile <- extractPP(phyloprofile, genomeName)
        if (!is.null(ppDir)) {
            prioPath <- paste(ppDir, setName, ".prioritylist", sep = "")
        } else {
            prioPath <- paste(
                coreDir, "output", "/", coreSet, "/", 
                as.character(scoreMode), "/", setName, ".prioritylist", sep = ""
            )
        }
        
        priorityTable <- read.table(
            prioPath,
            header = TRUE,
            sep = "\t"
        )
        priorityTable <- subset(priorityTable, genomeID == genomeName)
        priorityList <- strsplit(
            priorityTable[1, 2], ",", fixed = TRUE
        )[[1]]
        singleReport <- reportSingle(
            phyloprofile, coreDir, coreSet, scoreMode,
            priorityList
        )
        
        if (!is.null(ppDir)) {
            reportFile <- paste(ppDir, setName, ".report", sep = "")
        } else {
            reportFile <- paste(
                coreDir, "output", "/", coreSet, "/",
                as.character(scoreMode), "/", setName, ".report",
                sep = ""
            )
        }

        allReport <- read.table(
            reportFile,
            header = TRUE,
            sep = "\t"
        )
    }
    
    if (!is.null(outDir)) {
        write.table(
            singleReport,
            paste(outDir, setName, "_details.report", sep = ""),
            row.names = FALSE,
            quote = FALSE,
            sep = "\t"
        )
        write.table(
            allReport,
            paste(outDir, setName, ".report", sep = ""),
            row.names = FALSE,
            quote = FALSE,
            sep = "\t"
        )
    } else {
        write.table(
            singleReport,
            paste(coreDir, "output", "/", coreSet, "/", as.character(scoreMode), 
                  "/", setName, "_details.report", sep = ""),
            row.names = FALSE,
            quote = FALSE,
            sep = "\t"
        )
    }

    end <- Sys.time()
    print(paste("Running time is", as.character(end - start)))
    return(list(singleReport, allReport))
}