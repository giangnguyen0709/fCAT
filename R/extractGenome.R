#' This function takes a phylogenetic profile, which contains many different
#' genomes as an input and extract the phylogenetic profile of a specific genome
#'
#' @param pp A phylogenetic profile in data.frame
#' @param genome The genome ID of the genome, whose pp need to be extracted.
#' Exp: genome = HUMAN@9606@3
#'
#' @return The phylogenetic profile of the interested genome in data.frame
#' @export
extractPP <- function(pp, genome) {
    genomeID <- unlist(lapply(
        pp$orthoID,
        function(orthoID) {
            return(strsplit(orthoID, "|", fixed = TRUE)[[1]][2])
        }
    ))
    singlePP <- cbind(pp, genomeID)
    singlePP <- subset(singlePP, genomeID == genome)
    return(singlePP[1:(ncol(singlePP) - 1)])
}

#' A domains file can contains lines of many different genomes. This function
#' extract the domains of a specific genome
#'
#' @param domains The domains in data.frame
#' @param genome The genome ID of the genome, whose domains need to be extracted
#' @param reverse if reverse = TRUE then the function will not extract the lines
#' of the genome but remove it
#'
#' @return the domains of the interested genome in data.frame
#' @export
extractDomains <- function(domains, genome, reverse = FALSE) {
    checkGenome <- function(compareLine) {
        characterized <- as.character(compareLine)
        splited <- strsplit(characterized, "#", fixed = TRUE)[[1]]
        splited <- strsplit(splited[length(splited)], "|", fixed = TRUE)[[1]]
        genomeName <- splited[2]
        return(genomeName)
    }

    genomeName <- unlist(lapply(domains$V1, checkGenome))
    domains <- cbind(domains, genomeName)
    if (reverse == TRUE) {
        domains <- subset(domains, genomeName != genome)
    } else {
        domains <- subset(domains, genomeName == genome)
    }
    domains <- domains[c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8")]
    return(domains)
}

#' This function extract the sequences of a specific genome from a mulitple 
#' fasta file
#'
#' @param lines a vector that contains lines of an extended.fasta file
#' @param genome the name of the genome, that should be extracted
#' @return the vector, which contains the sequences of the interested genome
#' @export
extractFasta <- function(lines, genome) {
    newLines <- c()
    for (i in c(1:length(lines))) {
        if (i %% 2 == 1) {
            splited <- strsplit(lines[i], "|", fixed = TRUE)[[1]]
            spec <- splited[2]

            if (spec == genome) {
                newLines <- c(newLines, lines[i], lines[i + 1])
            }
        }
    }
    return(newLines)
}
