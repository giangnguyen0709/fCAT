#' A fail running of fCAT's checkCompleteness function will leave behind a
#' copy of the fasta file and the symbolic link of the FAS annotaion of the
#' interested genome. This can cause some warning or crash while rerun the 
#' tool. This function will check if the query taxon is empty before the check,
#' if not the function will delete the fasta file in the query taxon and its 
#' corresponding annotation file in weight_dir
#' @param root The path to the core directory, where the core set is stored
#' within weight_dir, blast_dir, etc.
#' @param weightDir The user can replace the weight_dir folder in the core
#' directory by specifying the path to the replacing folder in this argument
#'
#' @return none
#' @examples
#' coreFolder <- system.file("extdata", "sample", package = "fCAT")
#' handleError(coreFolder)
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
            list.dirs(genomeDir, full.names = FALSE, recursive = FALSE)
        ) != 0
    ) {
        for (
            directory in list.dirs(
                genomeDir,
                full.names = FALSE, recursive = FALSE
            )
        ) {
            unlink(
                paste(genomeDir, "/", directory, sep = ""), recursive = TRUE
            )
            if (
                file.exists(
                    paste(weightDir, "/", directory, ".json", sep = "")
                )
            ) {
                file.remove(
                    paste(weightDir, "/", directory, ".json", sep = "")
                )
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
#' @examples
#' coreFolder <- system.file("extdata", "sample", package = "fCAT")
#' checkPreProcess(coreFolder, "test")
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

#' The function to validate the input. If a input is not validated, it will
#' return a message to the user
#' 
#' @usage 
#' checkArguments(
#'     genome, fasAnno = NULL, root, coreSet, extend = FALSE,
#'     scoreMode, refSpecList = NULL, cpu = 4, blastDir = NULL, 
#'     weightDir = NULL, output = NULL, cleanup = FALSE, redo = FALSE
#' )
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
#' @param redo when redo is set to TRUE, all old data of the interested genome 
#' include fdogOutput, phyloprofileOutput, completenessOutput and the extended
#' phyloprofile will be removed and fCAT will recheck for this interested genome
#' @param scoreMode the mode determines the method to scoring the founded
#' ortholog and how to classify them. Choices: 1, 2, 3, "len"
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
#' @param output The directory which contains the output directory. It it is 
#' equal to NULL output will be set to working directory
#' @examples
#' coreFolder <- system.file("extdata", "sample", package = "fCAT")
#' genome <- system.file("extdata", "HUMAN@9606@3.fa", package = "fCAT")
#' fasAnno <- system.file("extdata", "HUMAN@9606@3.json", package = "fCAT")
#' checkArguments(
#'     genome[1], fasAnno[1], coreFolder[1], "test",
#'     scoreMode = 2, priorityList = c("HUMAN@9606@1")
#' )
#' @return A list that contains a logical value and the message to the error
#' @export
checkArguments <- function(
    genome, fasAnno = NULL, root, coreSet, extend = FALSE,
    scoreMode, refSpecList = NULL, cpu = 4, blastDir = NULL, weightDir = NULL,
    output = NULL, cleanup = FALSE, redo = FALSE) {
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


    modeList <- list(1, 2, 3, "len")
    if (!(scoreMode %in% modeList)) {
        check <- FALSE
        status <- "score mode is not available"
        return(list(check, status))
    }

    if (!is.null(refSpecList)) {
        if (!is.vector(refSpecList)) {
            check <- FALSE
            status <- "references species list must be a vector"
            return(list(check, status))
        }
        
        if (is.null(blastDir)) {
            coreTaxa <- list.dirs(
                paste(root, "blast_dir", sep = ""),
                recursive = FALSE,
                full.names = FALSE
            )
        } else {
            coreTaxa <- list.dirs(
                blastDir,
                recursive = FALSE,
                full.names = FALSE
            )
        }
        for (taxa in refSpecList) {
            if (!(taxa %in% coreTaxa)) {
                check <- FALSE
                status <-
                    "taxa in the priority list don't belong to the core set"
                return(list(check, status))
            }
        }
    }

    return(list(check, status))
}

#' The function to check completeness of a interested genome based on an
#' interested core set.
#' 
#' @usage 
#' checkCompleteness(
#'         genome, fasAnno = NULL, coreDir, coreSet, extend = FALSE, 
#'         scoreMode, refSpecList = NULL, cpu = 4, blastDir = NULL, 
#'         weightDir = NULL, output = NULL, cleanup = FALSE, redo = FALSE
#' )
#'
#' @param genome The path to the genome fasta file
#' @param fasAnno The path to the fas annotation file. It can equal NULL
#' @param coreDir The path to the core directory, where the core set is stored
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
#' @param redo when redo is set to TRUE, all old data of the interested genome 
#' include fdogOutput, phyloprofileOutput, completenessOutput and the extended
#' phyloprofile will be removed and fCAT will recheck for this interested genome
#' @param scoreMode the mode determines the method to scoring the founded
#' ortholog and how to classify them. Choices: 1, 2, 3, "busco"
#' @param refSpecList A list contains one or many genome ID of the genomes,
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
#'
#' @return A list which contains 2 data.frame. The first table is the
#' completeness report of the interested genome with details information about
#' the classification of the founded ortholog and which gene is missing. The
#' second table is the frequency table of the interested genome within the 
#' other genomes, which are present in the old phylogenetic profile. The 
#' frequency table give an general sight about how many "dissimilar", 
#' "similar", "duplicated" and "missing" genes founded in the interested 
#' genome.
#' @examples
#' coreFolder <- system.file("extdata", "sample", package = "fCAT")
#' genome <- system.file("extdata", "HUMAN@9606@3.fa", package = "fCAT")
#' fasAnno <- system.file("extdata", "HUMAN@9606@3.json", package = "fCAT")
#' checkCompleteness(
#'     genome[1], fasAnno[1], coreFolder[1], "test",
#'     scoreMode = 2, priorityList = c("HUMAN@9606@1"), extend = TRUE,
#'
#' )
#' @export
checkCompleteness <- function(
    genome, fasAnno = NULL, coreDir, coreSet, extend = FALSE,
    scoreMode, refSpecList = NULL, cpu = 4, blastDir = NULL, weightDir = NULL,
    output = NULL, cleanup = FALSE, redo = FALSE
) {
    argumentCheck <- checkArguments(
        genome, fasAnno, coreDir, coreSet, extend,
        scoreMode, refSpecList, cpu, blastDir, weightDir,
        output, cleanup, redo
    )
    if (argumentCheck[[1]] == FALSE) {
        return(argumentCheck[[2]])
    }
    start <- Sys.time()
    if (!endsWith(coreDir, "/")) {
        coreDir <- paste(coreDir, "/", sep = "")
    }
    if (is.null(output)) {
        output <- getwd()
    }
    if (!endsWith(output, "/")) {
        output <- paste(output, "/", sep = "")
    }
    outDir <- paste(output, "fcatOutput", "/", coreSet, "/", sep = "")
    
    if (is.null(refSpecList)) {
        query <- list.dirs(paste(coreDir, "query_taxon", sep = ""),
                           recursive = FALSE,
                           full.names = FALSE
        )[1]
        refSet <- list.dirs(paste(coreDir, "blast_dir", sep = ""),
                            recursive = FALSE,
                            full.names = FALSE
        )
        refSpecList <- autofindPriority(query, refSet)
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

    if (!checkPreProcess(coreDir, coreSet)) {
        processCoreSet(coreDir, coreSet)
    }

    handleError(coreDir)

    setName <- coreSet

    splited <- strsplit(genome, "/", fixed = TRUE)[[1]]
    splited <- splited[length(splited)]
    genomeName <- strsplit(splited, ".", fixed = TRUE)[[1]][1]
    
    redoCheck <- 0
    refSpecListCheck <- 0
    refSpecListFile <- paste(outDir, genomeName, "/", "refspec.list", sep = "")
    if (file.exists(refSpecListFile)) {
        if (redo == "TRUE") {
            redoCheck <- 1
        } else {
            oriRefSpecList <- readLines(refSpecListFile)
            if (oriRefSpecList != refSpecList) {
                refSpecListCheck <- 1
            }
        }
    }
    
    if (refSpecListCheck == 1 || redoCheck == 1) {
        if (refSpecListCheck == 1) {
            message <- paste(
                "It already exist data of", genomeName, 
                "with a different references", 
                "species list. Do you want to clean all", 
                "data of this genome in the", 
                "output directory. If yes, make sure",
                "that you saved your interested", 
                "data in a safe location [y/n]:")
        }
        if (redoCheck == 1) {
            message <- paste(
                "redo is set to TRUE.", "Do you want to clean all data",
                "of this genome in the", 
                "output directory. If yes, make sure", 
                "that you saved your interested", 
                "data in a safe location [y/n]:")
        }
        while (TRUE) {
            cleanCheck <- readline(prompt = message)
            if (cleanCheck == "y" || cleanCheck == "n") {
                break
            }
        }
        
        if (cleanCheck == "y") {
            unlink(paste(outDir, genomeName, sep = ""), recursive = TRUE)
            scoreModeList <- list(1, 2, "len")
            for (mode in scoreModeList) {
                correctFiles(
                    paste(outDir, "phyloprofile", sep = ""),
                    genomeName,
                    mode
                )
            }
        }
    }
    
    if (scoreMode == "len") {
        reportFile <- paste(
            outDir, genomeName, "/", "mode_len", "/", "full_text.txt", sep = ""
        )
        summaryFile <- paste(
            outDir, genomeName, "/", "mode_len", "/", "summary.txt", sep = ""
        )
    } else {
        reportFile <- paste(
            outDir, genomeName, "/", "mode", as.character(scoreMode), 
            "/", "full_text.txt", sep = ""
        )
        summaryFile <- paste(
            outDir, genomeName, "/", "mode", as.character(scoreMode), 
            "/", "summary.txt", sep = ""
        )
    }
    if (!file.exists(reportFile)) {
        computeReport(
            genome, fasAnno, coreDir, coreSet, extend, scoreMode, refSpecList,
            cpu, FALSE, blastDir, weightDir, output)
    }
    summaryTable <- read.table(
        summaryFile,
        header = TRUE,
        sep = "\t"
    )
    
    if (cleanup == TRUE) {
        fdogFolder <- paste(
            outDir, genomeName, "/", "fdogOutput", "/", genomeName, sep = ""
        )
        if (scoreMode == 1) {
            file <- "mode1"
        }
        if (scoreMode == 2 || scoreMode == 3) {
            file <- "other"
        }
        if (scoreMode == "len") {
            file <- "len"
        }
        pathPrefix <- paste(
            outDir, genomeName, "/", "phyloprofileOutput", "/", file, 
            sep = ""
        )
        suffixList <- c(
            ".phyloprofile", "_reverse.domains", "_forward.domains")
        for (suffix in suffixList) {
            if (file.exists(paste(pathPrefix, suffix, sep = ""))) {
                file.remove(paste(pathPrefix, suffix, sep = ""))
            }
        }
    }

    end <- Sys.time()
    print(paste("Running time is", as.character(end - start)))
    return(summaryTable)
}
