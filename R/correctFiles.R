#' If user choose the option to extend the phylogenetic profile of the
#' interested genome into the original phylogenetic profile, the tool must save
#' the priority list that the user entered to use this one to get the result for
#' the next calling. When the user want to redo the check, the function will
#' remove the old taxa in teh priorityFile
#'
#' @param priorityFile The path to the text file, that contains the information
#' about the priority list of all searching
#' @param genome the genome ID that need to be remove
#'
#' @return none
#' @export
correctPriority <- function(priorityFile, genome) {
    priorityTable <- read.table(priorityFile, header = TRUE, sep = "\t")
    priorityTable <- subset(priorityTable, genomeID != genome)
    write.table(priorityTable,
        priorityFile,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
    )
}

#' Remove a specific genome from the report file
#'
#' @param reportFile the path to the report file
#' @param genome the genome ID that need to be removed
#'
#' @return none
#' @export
correctReport <- function(reportFile, genome) {
    reportTable <- read.table(reportFile,
        header = TRUE,
        sep = "\t"
    )
    reportTable <- subset(reportTable, genomeID != genome)
    write.table(reportTable,
        reportFile,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
    )
}

#' remove lines of a specific genome from a domains file
#'
#' @param domainsFile path to the file
#' @param genome name of the genome in the file, that should be remove
#'
#' @return none
#' @export
correctDomains <- function(domainsFile, genome) {
    checkGenome <- function(compareLine) {
        characterized <- as.character(compareLine)
        splited <- strsplit(characterized, "#", fixed = TRUE)[[1]]
        splited <- strsplit(splited[length(splited)], "|", fixed = TRUE)[[1]]
        genomeName <- splited[2]
        return(genomeName)
    }

    domains <- read.table(domainsFile, sep = "\t", comment.char = "")
    genomeName <- unlist(lapply(domains$V1, checkGenome))
    domains <- cbind(domains, genomeName)
    domains <- subset(domains, genomeName != genome)
    domains <- domains[c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8")]
    write.table(
        domains, domainsFile, quote = FALSE, sep = "\t", 
        row.names = FALSE, col.names = FALSE)
}

#' remove the lines of a specific genome from a phylogenetic profile
#'
#' @param PPFile path to phylogenetic profile
#' @param genome name of the genome, that should be removed
#' @return none
#' @export
correctPP <- function(PPFile, genome) {
    checkGenome <- function(orthoID) {
        characterized <- as.character(orthoID)
        splited <- strsplit(characterized, "|", fixed = TRUE)[[1]]
        genomeName <- splited[2]
        return(genomeName)
    }
    PP <- read.table(PPFile, header = TRUE, sep = "\t")
    genomeName <- unlist(lapply(PP$orthoID, checkGenome))
    PP <- cbind(PP, genomeName)
    PP <- subset(PP, genomeName != genome)
    PP <- PP[1:(ncol(PP) - 1)]
    write.table(PP, PPFile, row.names = FALSE, quote = FALSE, sep = "\t")
}

#' remove all the lines of a specific genome from a fasta file
#'
#' @param fasta the path to the fasta file
#' @param genome the name of the genome, that should be removed
#' @return none
#' @export
correctFasta <- function(fasta, genome) {
    lines <- readLines(fasta)

    newLines <- c()
    for (i in c(1:length(lines))) {
        if (i %% 2 == 1) {
            splited <- strsplit(lines[i], "|", fixed = TRUE)[[1]]
            spec <- splited[2]

            if (spec != genome) {
                newLines <- c(newLines, lines[i], lines[i + 1])
            }
        }
    }
    if (length(newLines) != 0) {
        writeLines(text = newLines, con = fasta)
    } else {
        file.remove(fasta)
    }
}

#' remove all the lines of a specific genome from all of the phyloprofile,
#' domains, fasta files, priority file, and report file in a folder
#' @param directory the path to the directory, that contains the files
#' @param genome name of the genome, that should be removed
#' @return none
#' @export
correctFiles <- function(directory, genome) {
    files <- list.files(directory, full.names = TRUE, recursive = TRUE)

    for (file in files) {
        if (endsWith(file, ".phyloprofile")) {
            correctPP(file, genome)
        }

        if (endsWith(file, ".extended.fa")) {
            correctFasta(file, genome)
        }

        if (endsWith(file, ".domains")) {
            correctDomains(file, genome)
        }

        if (endsWith(file, "priority.list")) {
            correctPriority(file, genome)
        }

        if (endsWith(file, ".report")) {
            correctReport(file, genome)
        }
    }
}
