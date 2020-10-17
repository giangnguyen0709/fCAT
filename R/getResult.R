#' The function to classify a founded ortholog from its information like FAS
#' scores, the appearance frequence of the core gene in phylogenetic profile,
#' etc.
#'
#' @param fasF the forward fas score of the ortholog
#' @param fasB the backward fas score of the ortholog
#' @param root The path to the core directory, where the core set is stored
#' within weight_dir, blast_dir, etc.
#' @param coreSet The name of the interested core set. The core directory can
#' contains more than one core set and the user must specify the interested
#' core set. The core set will be stored in the folder core_orthologs in
#' subfolder, specify them by the name of the subfolder
#' @param coreGene the ID of the core group
#' @param scoreMode the mode determines the method to scoring the founded
#' ortholog and how to classify them. Choices: 1, 2, 3, "busco"
#' @param f the appearance frequence of the core gene in the phylogenetic
#' profile
#' @param priorityList A list contains one or many genome ID of the genomes,
#' which were used to build the core set. The genome ID of this list will be
#' stored with an priority order, the tool look at into the fasta file of each
#' core group and determine with the priority order to determine the references
#' species for each core group.
#' @return the status of the ortholog of the core gene. It can be "similar",
#' "dissimilar", "duplicated, similar" or "duplicated, dissimilar"
#' #' @examples
#' coreFolder <- system.file("extdata", "sample", package = "fCAT")
#' status <- assessStatus(
#'     0.7, 0.9, coreFolder, "test", "530670", 
#'     scoreMode = 2, f = 1, priorityList = c("HUMAN@9606@1"))
#' print(status)
#' @export
assessStatus <- function(
    fasF, fasB, root, coreSet, coreGene, scoreMode, f, priorityList) {
    if (scoreMode == 1) {
        fas <- (fasF + fasB) / 2
        cutoff <- read.table(
            paste(
                root, "core_orthologs", "/", coreSet, "/",
                coreGene, "/", "fas_dir", "/", "score_dir", "/",
                "1.cutoff",
                sep = ""
            ),
            header = TRUE,
            sep = "\t"
        )
        cutoffValue <- cutoff[1, 2]
    }

    if (scoreMode == 2) {
        fas <- (fasF + fasB) / 2
        refSpec <- getSpec(
            paste(root, "core_orthologs", "/", coreSet, "/",
                coreGene, "/", coreGene, ".fa",
                sep = ""
            ),
            priorityList
        )
        cutoff <- read.table(
            paste(
                root, "core_orthologs", "/", coreSet, "/", coreGene, "/",
                "fas_dir", "/", "score_dir", "/", "2.cutoff",
                sep = ""
            ),
            header = TRUE,
            sep = "\t"
        )
        subCutoff <- subset(cutoff, taxa == refSpec)
        cutoffValue <- subCutoff[1, 2]
    }

    if (scoreMode == 3) {
        fas <- (fasF + fasB) / 2
        refSpec <- getSpec(
            paste(root, "core_orthologs", "/", coreSet, "/",
                coreGene, "/", coreGene, ".fa",
                sep = ""
            ),
            priorityList
        )
        meanTable <- read.table(
            paste(
                root, "core_orthologs", "/", coreSet, "/", coreGene, "/",
                "fas_dir", "/", "score_dir", "/", "1.cutoff",
                sep = ""
            ),
            header = TRUE,
            sep = "\t"
        )
        lcl <- meanTable[2, 2]
        ucl <- meanTable[3, 2]
        if (fas < lcl || fas > ucl) {
            status <- "dissimilar"
        } else {
            status <- "similar"
        }
    }

    if (scoreMode != 3) {
        if (fas < cutoffValue) {
            status <- "dissimilar"
        } else {
            status <- "similar"
        }
    }

    if (f >= 2) {
        status <- paste("duplicated", ",", status)
    }

    return(status)
}

#' The function to calculate the cutoff value based on standard deviation of
#' the length for each core group. This function will be used in score mode
#' "busco"
#'
#' @param root The path to the core directory, where the core set is stored
#' within weight_dir, blast_dir, etc.
#' @param coreSet The name of the interested core set. The core directory can
#' contains more than one core set and the user must specify the interested
#' core set. The core set will be stored in the folder core_orthologs in
#' subfolder, specify them by the name of the subfolder
#' @param coreGene the ID of the core gene
#'
#' @return a list that contains the mean length and the standard deviation of
#' the length of the core gene
#' @examples
#' coreFolder <- system.file("extdata", "sample", package = "fCAT")
#' cutoff <- calculateBuscoCutoff(coreFolder, "test", "530670")
#' print(cutoff)
#' @export
calculateBuscoCutoff <- function(root, coreSet, coreGene) {
    if (!endsWith(root, "/")) {
        root <- paste(root, "/", sep = "")
    }

    seedFasta <- readLines(
        paste(
            root, "core_orthologs", "/", coreSet, "/",
            coreGene, "/", coreGene, ".fa",
            sep = ""
        )
    )
    i <- 1:length(seedFasta)
    seedFasta <- seedFasta[i[i %% 2 == 0]]
    lengthSet <- lapply(
        seedFasta,
        function(seq) {
            return(nchar(seq))
        }
    )
    lengthSet <- unlist(lengthSet)

    meanLength <- mean(lengthSet)
    standardDeviation <- sqrt(mean((lengthSet - meanLength)**2))

    return(list(meanLength, standardDeviation))
}

#' The function to classify the founded ortholog based on its length
#'
#' @param orthoLength the length of the ortholg sequence
#' @param root The path to the core directory, where the core set is stored
#' within weight_dir, blast_dir, etc.
#' @param coreSet The name of the interested core set. The core directory can
#' contains more than one core set and the user must specify the interested
#' core set. The core set will be stored in the folder core_orthologs in
#' subfolder, specify them by the name of the subfolder
#' @param coreGene the ID of the core gene
#' @param f the appearance frequence of the core gene in the phylogenetic
#' profile
#'
#' @return the status of the core gene, it can be "fragmented", "complete",
#' "duplicated, complete" or "duplicated, fragmented"
#' @examples
#' coreFolder <- system.file("extdata", "sample", package = "fCAT")
#' status <- assessBusco(600, coreFolder, "test", "530670", f = 1)
#' print(status)
#' @export
assessBusco <- function(orthoLength, root, coreSet, coreGene, f) {
    cutoff <- calculateBuscoCutoff(root, coreSet, coreGene)
    if (cutoff[[2]] != 0) {
        score <- (orthoLength - cutoff[[1]]) / cutoff[[2]]
        if (score > 2 || score < (-2)) {
            status <- "fragmented"
        } else {
            status <- "complete"
        }
    } else {
        score <- orthoLength - cutoff[[1]]
        if (score == 0) {
            status <- "complete"
        } else {
            status <- "fragmented"
        }
    }
    if (f >= 2) {
        status <- paste("duplicated", ",", status)
    }

    return(list(status, cutoff[[1]], cutoff[[2]]))
}

#' Determine if a core gene was ignored by the tool because of the unknown
#' references species
#'
#' @param root The path to the core directory, where the core set is stored
#' within weight_dir, blast_dir, etc.
#' @param coreSet The name of the interested core set. The core directory can
#' contains more than one core set and the user must specify the interested
#' core set. The core set will be stored in the folder core_orthologs in
#' subfolder, specify them by the name of the subfolder
#' @param coreGene the ID of the core gene
#' @param priorityList A list contains one or many genome ID of the genomes,
#' which were used to build the core set. The genome ID of this list will be
#' stored with an priority order, the tool look at into the fasta file of each
#' core group and determine with the priority order to determine the references
#' species for each core group.
#'
#' @return TRUE or FALSE
#' @examples
#' coreFolder <- system.file("extdata", "sample", package = "fCAT")
#' check <- filterIgnore(coreFolder, "test", "530670", c("HUMAN@9606@1"))
#' print(check)
#' @export
filterIgnore <- function(root, coreSet, coreGene, priorityList) {
    if (!endsWith(root, "/")) {
        root <- paste(root, "/", sep = "")
    }
    fasta <- paste(root, "core_orthologs", "/", coreSet, "/", coreGene, "/",
        coreGene, ".fa",
        sep = ""
    )
    check <- getSpec(fasta, priorityList)
    if (is.null(check)) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}

#' Create the report of the completeness of a genome based on its 
#' phylogenetic profile to the interested core set
#'
#' @param pp the phylogenetic profile of the genome in data.frame
#' @param root The path to the core directory, where the core set is stored
#' within weight_dir, blast_dir, etc.
#' @param coreSet The name of the interested core set. The core directory can
#' contains more than one core set and the user must specify the interested
#' core set. The core set will be stored in the folder core_orthologs in
#' subfolder, specify them by the name of the subfolder
#' @param scoreMode the mode determines the method to scoring the founded
#' ortholog and how to classify them. Choices: 1, 2, 3, "busco"
#' @param priorityList A list contains one or many genome ID of the genomes,
#' which were used to build the core set. The genome ID of this list will be
#' stored with an priority order, the tool look at into the fasta file of each
#' core group and determine with the priority order to determine the references
#' species for each core group.
#'
#' @return the report in data frame
#' @examples
#' ## Create a pseudo phylogenetic profile table
#' geneID <- c("530670", "530730")
#' ncbiID <- c("ncbi9606", "ncbi9606")
#' orthoID <- c("530670|HUMAN@9606@3|Q16526|1", "530730|HUMAN@9606@3|P05091|1")
#' FAS_F <- c(1, 1)
#' FAS_B <- c(1, 1)
#' pp <- data.frame(geneID, ncbiID, orthoID, FAS_F, FAS_B)
#' coreFolder <- system.file("extdata", "sample", package = "fCAT")
#' ## Translate phylogenetic profile to a detailed report
#' report <- reportSingle(pp, coreFolder, "test", 2, c("HUMAN@9606@1"))
#' print.data.frame(report)
#' @export
reportSingle <- function(pp, root, coreSet, scoreMode, priorityList) {
    if (!endsWith(root, "/")) {
        root <- paste(root, "/", sep = "")
    }

    coreGeneList <- list.dirs(
        paste(root, "core_orthologs", "/", coreSet, sep = ""),
        recursive = FALSE,
        full.names = FALSE
    )
    frequency <- table(pp$geneID)
    if (scoreMode != "len") {
        status <- unlist(lapply(1:nrow(pp),
            function(
                i, frequency, pp, scoreMode, root, coreSet, priorityList
            ) {
                fasF <- pp[i, 4]
                fasB <- pp[i, 5]
                coreGene <- pp[i, 1]
                f <- frequency[as.character(coreGene)]
                s <- assessStatus(
                    fasF, fasB, root, coreSet,
                    coreGene, scoreMode, f,
                    priorityList
                )
                return(s)
            },
            frequency = frequency,
            pp = pp,
            scoreMode = scoreMode,
            root = root,
            coreSet = coreSet,
            priorityList = priorityList
        ))
    } else {
        status <- lapply(1:nrow(pp),
            function(i, frequency, pp, root, coreSet) {
                orthoLength <- pp[i, 4]
                coreGene <- pp[i, 1]
                f <- frequency[as.character(coreGene)]
                info <- assessBusco(
                    orthoLength, root, coreSet,
                    coreGene, f
                )
                r <- data.frame(
                    status = c(info[[1]]),
                    mean_length = c(info[[2]]),
                    standard_deviation = c(info[[3]])
                )
                return(r)
            },
            frequency = frequency,
            pp = pp,
            root = root,
            coreSet = coreSet
        )
        status <- do.call("rbind", status)
    }
    missingGene <- setdiff(coreGeneList, unique(pp$geneID))
    if (length(missingGene) != 0) {
        missingStatus <- unlist(lapply(missingGene,
            function(gene, root, coreSet, priorityList) {
                check <- filterIgnore(root, coreSet, gene, priorityList)
                if (check == TRUE) {
                    return("ignored")
                } else {
                    return("missing")
                }
            },
            root = root,
            coreSet = coreSet,
            priorityList = priorityList
        ))
        if (scoreMode != "len") {
            report <- data.frame(
                geneID = pp$geneID, orthoID = pp$orthoID, status,
                FAS_F = pp$FAS_F, FAS_B = pp$FAS_B
            )
            missingTable <- data.frame(
                geneID = missingGene, orthoID = NA,
                status = missingStatus, FAS_F = NA, FAS_B = NA
            )
            report <- rbind(report, missingTable)
        } else {
            report <- data.frame(
                geneID = pp$geneID, orthoID = pp$orthoID,
                status = status$status, length = pp$length,
                mean_length = status$mean_length,
                standard_deviation = status$standard_deviation
            )
            missingTable <- data.frame(
                geneID = missingGene, orthoID = NA,
                status = missingStatus, length = NA,
                mean_length = NA,
                standard_deviation = NA
            )
            report <- rbind(report, missingTable)
        }
    } else {
        missingTable <- NULL
        if (scoreMode != "len") {
            report <- data.frame(
                geneID = pp$geneID, orthoID = pp$orthoID, status,
                FAS_F = pp$FAS_F, FAS_B = pp$FAS_B
            )
        } else {
            report <- data.frame(
                geneID = pp$geneID, orthoID = pp$orthoID,
                status = status$status, length = pp$length,
                mean_length = status$mean_length,
                standard_deviation = status$standard_deviation
            )
        }
    }
    return(list(report, missingTable))
}

#' Translate the report table into a frequence table, which tell the user, 
#' how many "dissimilar", "simillar", "missing", etc.
#'
#' @param genomeID the genome ID of the interested genome
#' @param report the report of the completeness of the interested genome based
#' on its phylogenetic profile
#' @param scoreMode the mode determines the method to scoring the founded
#' ortholog and how to classify them. Choices: 1, 2, 3, "busco"
#'
#' @return A frequency table of the completeness of the interested genome in
#' data.frame
#' @examples
#' ## Create a pseudo phylogenetic profile table
#' geneID <- c("530670", "530730")
#' ncbiID <- c("ncbi9606", "ncbi9606")
#' orthoID <- c("530670|HUMAN@9606@3|Q16526|1", "530730|HUMAN@9606@3|P05091|1")
#' FAS_F <- c(1, 1)
#' FAS_B <- c(1, 1)
#' pp <- data.frame(geneID, ncbiID, orthoID, FAS_F, FAS_B)
#' coreFolder <- system.file("extdata", "sample", package = "fCAT")
#' ## Translate phylogenetic profile to a detailed report
#' report <- reportSingle(pp, coreFolder, "test", 2, c("HUMAN@9606@1"))
#' ## Translate the report table to a frequency table
#' frequencyTable <- translateReport("HUMAN@9606@3", report, 2)
#' print.data.frame(frequencyTable)
#' @export
translateReport <- function(genomeID, report, scoreMode) {
    if (scoreMode == "len") {
        frequency <- table(report$status)

        complete <- frequency["complete"]
        fragmented <- frequency["fragmented"]
        missing <- frequency["missing"]
        ignored <- frequency["ignored"]

        if (is.na(complete)) {
            complete <- 0
        }
        if (is.na(fragmented)) {
            fragmented <- 0
        }
        if (is.na(missing)) {
            missing <- 0
        }
        if (is.na(ignored)) {
            ignored <- 0
        }

        duplicated <- length(unique(report$geneID)) -
            complete - fragmented - missing - ignored
        translated <- data.frame(
            genomeID = c(genomeID),
            complete = c(complete),
            fragmented = c(fragmented),
            missing = c(missing),
            duplicated = c(duplicated),
            ignored = c(ignored)
        )
        return(translated)
    } else {
        frequency <- table(report$status)

        similar <- frequency["similar"]
        dissimilar <- frequency["dissimilar"]
        missing <- frequency["missing"]
        ignored <- frequency["ignored"]
        if (is.na(similar)) {
            similar <- 0
        }
        if (is.na(dissimilar)) {
            dissimilar <- 0
        }
        if (is.na(missing)) {
            missing <- 0
        }
        if (is.na(ignored)) {
            ignored <- 0
        }

        duplicated <- length(unique(report$geneID)) - similar -
            dissimilar - missing - ignored
        translated <- data.frame(
            genomeID = c(genomeID),
            similar = c(similar),
            dissimilar = c(dissimilar),
            missing = c(missing),
            duplicated = c(duplicated),
            ignored = c(ignored)
        )
        return(translated)
    }
}
