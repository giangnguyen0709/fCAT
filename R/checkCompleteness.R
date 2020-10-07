#' A fail running of fCAT's checkCompleteness function will leave behind a copy
#' of the fasta file and the symbolic link of the FAS annotaion of the 
#' interested genome. This can cause some warning or crash while rerun the tool.
#' This function will check if the query taxon is empty before the check, if not
#' the function will delete the fasta file in the query taxon and its 
#' corresponding annotation file in weight_dir
#' @param root The path to the core directory, where the core set is stored 
#' within weight_dir, blast_dir, etc.
#' @param weightDir The user can replace the weight_dir folder in the core 
#' directory by specifying the path to the replacing folder in this argument
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

#' Check if the core set was processed with the function processCoreSet()
#'
#' @param root The path to the core directory, where the core set is stored 
#' within weight_dir, blast_dir, etc.
#' @param coreSet The name of the interested core set. The core directory can 
#' contains more than one core set and the user must specify the interested 
#' core set. The core set will be stored in the folder core_orthologs in 
#' subfolder, specify them by the name of the subfolder
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

#' Check if ID of the interested genome exists already in the original 
#' phylogenetic profile in the output folder of the core set or in the user's
#' phyloprofile folder
#'
#' @param genomeName the ID of the genome
#' @param root The path to the core directory, where the core set is stored 
#' within weight_dir, blast_dir, etc.
#' @param coreSet The name of the interested core set. The core directory can 
#' contains more than one core set and the user must specify the interested 
#' core set. The core set will be stored in the folder core_orthologs in 
#' subfolder, specify them by the name of the subfolder
#' @param scoreMode the mode determines the method to scoring the founded 
#' ortholog and how to classify them. Choices: 1, 2, 3, "busco"
#' @param ppDir The user can replace the default folder output in the core 
#' directory, where the phylogenetic profiles are stored by his folder. The user
#' can specify the path to his folder in this argument
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

#' The function to validate the input. If a input is not validated, it will 
#' return a message to the user
#'
#' @param genome The path to the genome fasta file
#' @param fasAnno The path to the fas annotation file. It can equal NULL
#' @param root The path to the core directory, where the core set is stored 
#' within weight_dir, blast_dir, etc.
#' @param coreSet The name of the interested core set. The core directory can 
#' contains more than one core set and the user must specify the interested 
#' core set. The core set will be stored in the folder core_orthologs in 
#' subfolder, specify them by the name of the subfolder
#' @param extend The output of the function is a phylogenetic profile of the 
#' interested genome. It contains 4 files, .phyloprofile, .extended.fa, 
#' _reverse.domains and _forward.domains. If extend = TRUE, the files will be 
#' appended into the old files in the folder output of the core directory or in 
#' the inputed folder by the user with the argument ppDir. If there is no old 
#' files in the folder, the output files of the function will be writen in the 
#' new files.
#' @param redo If it exists already the genome ID of the interested genome in 
#' the old phylogenetic profile. The tool will extract direct this pp to assess
#' the completeness. If user don't want this happens, they can set redo to TRUE
#' to get a new phylogenetic profile
#' @param scoreMode the mode determines the method to scoring the founded 
#' ortholog and how to classify them. Choices: 1, 2, 3, "busco"
#' @param priorityList A list contains one or many genome ID of the genomes,
#' which were used to build the core set. The genome ID of this list will be 
#' stored with an priority order, the tool look at into the fasta file of each
#' core group and determine with the priority order to determine the references
#' species for each core group.
#' @param cpu determines the cores that fDOG and fdogFAS will uses to be run 
#' parallel
#' @param blastDir The user can replace the blast_dir folder in the core
#' directory by specifying it in this argument
#' @param weightDir The user can replace the weight_dir folder in the core
#' directory by specifying it in this argument
#' @param outDir The user can specify the directory to save the output report
#' file of the completeness of the interested genome by specifying the path
#' to the folder in this argument
#' @param cleanup The fDOG's output is a set of phylogenetic profile of each
#' core group to the interested genome. The phylogenetic profile will be stored
#' into a folder in the core set. The function will merge all the small 
#' phylogenetic profile, calculate the FAS score or length to have the whole 
#' phylogenetic profile of the interested genome to the core set. This fDOG's
#' output can be reused for all score modes. When cleanup is set to TRUE, the
#' fDOG's output will not be stored to be reused but to be removed
#' @param reFdog If it already exist a fDOG's output for a specific core group
#' the tool will skip this core group and go to the next core group. If reFdog 
#' is set to TRUE, the tool will remove all the existed fDOG's output and rerun
#' fDOG for all core groups of the set
#' @param fdogDir Normally the fDOG's output will be stored in the folder 
#' fdogout in the core directory, but the user can specify the folder for 
#' fDOG's output by specify the path to it in this argument. Notice here, is 
#' the fDOG's output folder will contains the subfolder, equivalent to the name
#' of the interested genome, for example, the folder can contain "HUMAN@9606@3" 
#' and "AMPQU@400682@2", for a completeness checking on an interested genome,
#' which has a subfolder in the fDOG's output folder with the same name, the 
#' function will look into the subfolder to find the existed fDOG's output
#' @param ppDir The user can replace the default folder output in the core 
#' directory, where the phylogenetic profiles are stored by his folder. The user
#' can specify the path to his folder in this argument
#'
#' @return A list that contains a logical value and the message to the error
#' @export
checkArguments <- function(
    genome, fasAnno = NULL, root, coreSet, extend = FALSE, redo = FALSE, 
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

#' The function to check completeness of a interested genome based on an
#' interested core set.
#'
#' @param genome The path to the genome fasta file
#' @param fasAnno The path to the fas annotation file. It can equal NULL
#' @param root The path to the core directory, where the core set is stored 
#' within weight_dir, blast_dir, etc.
#' @param coreSet The name of the interested core set. The core directory can 
#' contains more than one core set and the user must specify the interested 
#' core set. The core set will be stored in the folder core_orthologs in 
#' subfolder, specify them by the name of the subfolder
#' @param extend The output of the function is a phylogenetic profile of the 
#' interested genome. It contains 4 files, .phyloprofile, .extended.fa, 
#' _reverse.domains and _forward.domains. If extend = TRUE, the files will be 
#' appended into the old files in the folder output of the core directory or in 
#' the inputed folder by the user with the argument ppDir. If there is no old 
#' files in the folder, the output files of the function will be writen in the 
#' new files.
#' @param redo If it exists already the genome ID of the interested genome in 
#' the old phylogenetic profile. The tool will extract direct this pp to assess
#' the completeness. If user don't want this happens, they can set redo to TRUE
#' to get a new phylogenetic profile
#' @param scoreMode the mode determines the method to scoring the founded 
#' ortholog and how to classify them. Choices: 1, 2, 3, "busco"
#' @param priorityList A list contains one or many genome ID of the genomes,
#' which were used to build the core set. The genome ID of this list will be 
#' stored with an priority order, the tool look at into the fasta file of each
#' core group and determine with the priority order to determine the references
#' species for each core group.
#' @param cpu determines the cores that fDOG and fdogFAS will uses to be run 
#' parallel
#' @param blastDir The user can replace the blast_dir folder in the core
#' directory by specifying it in this argument. By default is NULL
#' @param weightDir The user can replace the weight_dir folder in the core
#' directory by specifying it in this argument. By default is NULL
#' @param outDir The user can specify the directory to save the output report
#' file of the completeness of the interested genome by specifying the path
#' to the folder in this argument. By default is NULL
#' @param cleanup The fDOG's output is a set of phylogenetic profile of each
#' core group to the interested genome. The phylogenetic profile will be stored
#' into a folder in the core set. The function will merge all the small 
#' phylogenetic profile, calculate the FAS score or length to have the whole 
#' phylogenetic profile of the interested genome to the core set. This fDOG's
#' output can be reused for all score modes. When cleanup is set to TRUE, the
#' fDOG's output will not be stored to be reused but to be removed
#' @param reFdog If it already exist a fDOG's output for a specific core group
#' the tool will skip this core group and go to the next core group. If reFdog 
#' is set to TRUE, the tool will remove all the existed fDOG's output and rerun
#' fDOG for all core groups of the set
#' @param fdogDir Normally the fDOG's output will be stored in the folder 
#' fdogout in the core directory, but the user can specify the folder for 
#' fDOG's output by specify the path to it in this argument. Notice here, is 
#' the fDOG's output folder will contains the subfolder, equivalent to the name
#' of the interested genome, for example, the folder can contain "HUMAN@9606@3" 
#' and "AMPQU@400682@2", for a completeness checking on an interested genome,
#' which has a subfolder in the fDOG's output folder with the same name, the 
#' function will look into the subfolder to find the existed fDOG's output
#' @param ppDir The user can replace the default folder output in the core 
#' directory, where the phylogenetic profiles are stored by his folder. The user
#' can specify the path to his folder in this argument. By default is NULL
#'
#' @return A list which contains 2 data.frame. The first table is the 
#' completeness report of the interested genome with details information about
#' the classification of the founded ortholog and which gene is missing. The
#' second table is the frequency table of the interested genome within the other
#' genomes, which are present in the old phylogenetic profile. The frequency
#' table give an general sight about how many "dissimilar", "similar", 
#' "duplicated" and "missing" genes founded in the interested genome.
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
        rpDir <- paste(
            coreDir, "output", "/", coreSet, "/", as.character(scoreMode), 
            "/", "report", sep = "")
        if (!dir.exists(rpDir)) {
            dir.create(rpDir)
        }
        write.table(
            singleReport,
            paste(coreDir, "output", "/", coreSet, "/", as.character(scoreMode), 
                  "/", "report", "/", setName, "_details.report", sep = ""),
            row.names = FALSE,
            quote = FALSE,
            sep = "\t"
        )
        write.table(
            allReport,
            paste(coreDir, "output", "/", coreSet, "/", as.character(scoreMode), 
                  "/", "report", "/", setName, ".report", sep = ""),
            row.names = FALSE,
            quote = FALSE,
            sep = "\t"
        )
    }

    end <- Sys.time()
    print(paste("Running time is", as.character(end - start)))
    return(list(singleReport, allReport))
}