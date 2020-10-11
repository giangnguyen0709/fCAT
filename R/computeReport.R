#' The function to append the frequency table of the completeness classification
#' of the interested genome into the original frequency table of the core set.
#' The original frequency table describe the frequencies of each classification
#' for all genomes, which were checked the completeness with fCAT and the extend
#' was set to TRUE while checking.(How many "complete", "fragmented", "missing",
#' and "duplicated" for busco mode or how many "similar", "dissimilar", 
#' "missing", and "duplicated" for the other modes)
#'
#' @param report A data frame that contains the frequency table of the 
#' interested genome
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
#' @return none
#' @export
printReport <- function(report, root, coreSet, scoreMode, ppDir) {
    setName <- coreSet

    if (!is.null(ppDir)) {
        reportFile <- paste(ppDir, setName, ".report", sep = "")
    } else {
        reportFile <- paste(root, "output", "/", coreSet, "/",
                            as.character(scoreMode), "/", setName, ".report",
                            sep = ""
        )
    }

    if (!file.exists(reportFile)) {
        write.table(report,
            reportFile,
            row.names = FALSE,
            quote = FALSE,
            sep = "\t"
        )
    } else {
        write.table(report,
            reportFile,
            row.names = FALSE,
            col.names = FALSE,
            append = TRUE,
            quote = FALSE,
            sep = "\t"
        )
    }
}

#' Add the priority list of the interested genome to the priority file. If the
#' file doesnt exist it will be created
#'
#' @param genomeID the genome ID of the interested genome
#' @param priorityList A list contains one or many genome ID of the genomes,
#' which were used to build the core set. The genome ID of this list will be 
#' stored with an priority order, the tool look at into the fasta file of each
#' core group and determine with the priority order to determine the references
#' species for each core group.
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
#' @return none
#' @export
printPriority <- function(
    genomeID, priorityList, root, coreSet, scoreMode, ppDir
) {
    table <- data.frame(
        genomeID = c(genomeID),
        priority_list = c(paste(priorityList, collapse = ","))
    )
    
    if (!is.null(ppDir)) {
        priorityFile <- paste(ppDir, coreSet, ".prioritylist", sep = "")
    } else {
        priorityFile <- paste(root, "output", "/", coreSet, "/",
                              as.character(scoreMode), "/", coreSet, 
                              ".prioritylist", sep = ""
        )
    }

    if (!file.exists(priorityFile)) {
        write.table(
            table,
            priorityFile,
            row.names = FALSE,
            quote = FALSE,
            sep = "\t"
        )
    } else {
        write.table(
            table,
            priorityFile,
            row.names = FALSE,
            col.names = FALSE,
            append = TRUE,
            quote = FALSE,
            sep = "\t"
        )
    }
}

#' The function take a path to a genome fasta file and its FAS annotation file
#' (if the annotation file is not provided, the funtion will compute it) and
#' compute the detailed report of the completeness of the genome
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
#' @return a report of the completeness of the interested genome with detailed 
#' information of every core genes in the core set. Which core gene is "similar"
#' , which is "dissimilar", "duplicated" and "missing"
#' @export
computeReport <- function(genome, fasAnno, root, coreSet, extend = FALSE, 
    scoreMode, priorityList = NULL, cpu, computeOri = FALSE,
    blastDir = NULL, weightDir = NULL, cleanup = FALSE, 
    reFdog = FALSE, fdogDir = NULL, ppDir = NULL) {
    if (!endsWith(root, "/")) {
        root <- paste(root, "/", sep = "")
    }
    
    if (!is.null(ppDir)) {
        modeFolder <- ppDir
    }
    modeFolder <- paste(root, "output", "/", coreSet, "/",
                        as.character(scoreMode),
                        sep = ""
    )
    
    if (!dir.exists(modeFolder)) {
        dir.create(modeFolder, recursive = TRUE)
    }

    splited <- strsplit(genome, "/", fixed = TRUE)[[1]]
    splited <- splited[length(splited)]
    genomeName <- strsplit(splited, ".", fixed = TRUE)[[1]][1]


    placeSeed(genome, fasAnno, root, computeOri, weightDir)

    if (is.null(priorityList)) {
        query <- list.dirs(paste(root, "query_taxon", sep = ""),
            recursive = FALSE,
            full.names = FALSE
        )[1]
        refSet <- list.dirs(paste(root, "blast_dir", sep = ""),
            recursive = FALSE,
            full.names = FALSE
        )
        priorityList <- autofindPriority(query, refSet)
    }

    if (extend == TRUE) {
        printPriority(genomeName, priorityList, root, coreSet, scoreMode, ppDir)
    }

    if (scoreMode == "busco") {
        pp <- runFdogBusco(
            root, coreSet, extend, scoreMode, priorityList, cpu,
            blastDir, weightDir, cleanup, reFdog, fdogDir, ppDir)
    } else {
        pp <- runFdog(
            root, coreSet, extend, scoreMode, priorityList, cpu,
            blastDir, weightDir, cleanup, reFdog, fdogDir, ppDir)
    }

    report <- reportSingle(pp, root, coreSet, scoreMode, priorityList)
    if (extend == TRUE) {
        translated <- translateReport(genomeName, report, scoreMode)
        printReport(translated, root, coreSet, scoreMode, ppDir)
    }

    unlink(
        paste(root, "query_taxon", "/", genomeName, sep = ""), 
        recursive = TRUE
    )
    if (computeOri == FALSE) {
        if (!is.null(weightDir)) {
            if (!endsWith(weightDir, "/")) {
                weightDir <- paste(weightDir, "/", sep = "")
            }
            file.remove(
                paste(weightDir, genomeName, ".json", sep = "")
            )
        } else {
            file.remove(
                paste(root, "weight_dir", "/", genomeName, ".json", sep = "")
            )
        }
    }
    return(report)
}
