#' If user choose the option to extend the phylogenetic profile of the
#' interested genome into the original phylogenetic profile, the tool must save
#' the priority list in a file text in the core directory, so that the tool can
#' recall it to assess the completeness of the genome, which is appended in the 
#' original phylogenetic profile. When the user want to redo the check, it can 
#' be that the user can change the priority list, because of that, anytime a
#' redo is on, the tool must first remove the old priority list of the genome in
#' the file.
#'
#' @param priorityFile The path to the prioritylist file, which contains the
#' priority lists, which were used for the completeness checking of the genomes,
#' that were checked with the tool over the time.
#' @param genome The genome ID of the genome, that need to be removed from the
#' prioritylist file
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

#' For every with the tool checked genomes, fCAT will save the frequency table 
#' of each genome into a file with ending .report in the phyloprofile folder, 
#' which by default is the output folder in core directory or can be the
#' folder specified by the user with the argument ppDir in the function 
#' checkCompleteness. When redo is on to recheck for a genome, whose pp exists
#' already in the original pp, the tool must remove the lines of the
#' genome from the original frequency table. That is the purpose of this 
#' function
#'
#' @param reportFile The path to the reprot file (or the frequency table), which
#' contains the priority lists, which were used for the completeness checking of
#' the genomes, that were checked with the tool over the time.
#' @param genome the genome ID that need to be removed
#'
#' @return none
#' @export
correctReport <- function(reportFile, genome) {
    reportTable <- read.table(
        reportFile,
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

#' For all with the tool checked genomes, fCAT will save the domains
#' of each genome into a domains file in the phyloprofile folder, 
#' which by default is the output folder in core directory or can be the
#' folder specified by the user with the argument ppDir in the function 
#' checkCompleteness. When redo is on to recheck for a genome, whose pp exists
#' already in the original pp, the tool must remove the lines of the
#' genome from the original domains file. That is the purpose of this 
#' function
#'
#' @param domainsFile path to the domains file
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

##' For every with the tool checked genomes, fCAT will save the phylogenetic 
#' profile of each genome into a extended fasta file in the phyloprofile folder, 
#' which by default is the output folder in core directory or can be the
#' folder specified by the user with the argument ppDir in the function 
#' checkCompleteness. When redo is on to recheck for a genome, whose pp exists
#' already in the original pp, the tool must remove the lines of the
#' genome from the original phyloprofile file. That is the purpose of this 
#' function
#'
#' @param PPFile path to phylogenetic profile file
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

#' For every with the tool checked genomes, fCAT will save the sequences
#' of each genome into a extended fasta file in the phyloprofile folder, 
#' which by default is the output folder in core directory or can be the
#' folder specified by the user with the argument ppDir in the function 
#' checkCompleteness. When redo is on to recheck for a genome, whose pp exists
#' already in the original pp, the tool must remove the lines of the
#' genome from the original extended fasta file. That is the purpose of this 
#' function
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

#' This function run correctPP, correctFasta, correctDomains, correctPriority 
#' and correctReport in a loop to correct all phyloprofile, extended fasta, 
#' domains, prioritylist and report file in a folder
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

        if (endsWith(file, ".prioritylist")) {
            correctPriority(file, genome)
        }

        if (endsWith(file, ".report") && !endsWith(file, "_details.report")) {
            correctReport(file, genome)
        }
    }
}
