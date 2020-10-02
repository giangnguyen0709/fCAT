#' Create the report of the completeness of the interested genome with the other
#' genome in the original phylogenetic profile
#'
#' @param report the data frame that contains the report
#' @param root the path to the root folder
#' @param coreSet the core set name
#' @param scoreMode the mode determines the way to assess the founded ortholog
#'
#' @return none
#' @export
printReport <- function(report, root, coreSet, scoreMode, outDir=NULL) {
    setName <- coreSet

    if (!is.null(outDir)) {
        if (!endsWith(outDir, "/")) {
            outDir <- paste(outDir, "/", sep = "")
        }
        reportFile <- paste(outDir, setName, ".report")
    }
    reportFile <- paste(root, "output", "/", coreSet, "/",
        as.character(scoreMode), "/", setName, ".report",
        sep = ""
    )

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
#' @param genomeID the genomeID of the interested genome
#' @param priorityList the priority list
#' @param root the path to the root folder
#' @param coreSet the core set name
#' @param scoreMode the mode determines the method to assess the founded
#' ortholog
#'
#' @return none
#' @export
printPriority <- function(
    genomeID, priorityList, root, coreSet, scoreMode, outDir = NULL
) {
    table <- data.frame(
        genomeID = c(genomeID),
        priority_list = c(paste(priorityList, collapse = ","))
    )
    if (!is.null(outDir)) {
        if (!endsWith(outDir, "/")) {
            outDir <- paste(outDir, "/", sep = "")
        }
        priorityFile <- paste(outDir, coreSet, ".prioritylist", sep = "")
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
computeReport <- function(genome, fasAnno, root, coreSet, extend = FALSE, 
    scoreMode, priorityList = NULL, cpu, computeOri = FALSE,
    blastDir = NULL, weightDir = NULL, outDir = NULL, cleanup = FALSE) {
    if (!endsWith(root, "/")) {
        root <- paste(root, "/", sep = "")
    }
    if (!is.null(outDir)) {
        modeFolder <- outDir
    } else {
        modeFolder <- paste(root, "output", "/", coreSet, "/",
                            as.character(scoreMode),
                            sep = ""
        )
    }
    
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
        printPriority(genomeName, priorityList, root, coreSet, scoreMode)
    }

    if (scoreMode == "busco") {
        pp <- runFdogBusco(
            root, coreSet, extend, scoreMode, priorityList, cpu,
            blastDir, weightDir, outDir, cleanup)
    } else {
        pp <- runFdog(
            root, coreSet, extend, scoreMode, priorityList, cpu,
            blastDir, weightDir, outDir, cleanup)
    }

    report <- reportSingle(pp, root, coreSet, scoreMode, priorityList)
    if (extend == TRUE) {
        translated <- translateReport(genomeName, report, scoreMode)
        printReport(translated, root, coreSet, scoreMode)
    }

    unlink(
        paste(root, "query_taxon", "/", genomeName, sep = ""), 
        recursive = TRUE
    )
    if (computeOri == FALSE) {
        if (!is.null(weightDir)) {
            if (!endsWith(weightDir)) {
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
