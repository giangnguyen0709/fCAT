#' The function extract a phylogenetic profile of a genome from a phylogentic
#' profile
#'
#' @param pp the data frame of the phylogenetic profile
#' @param genome the ID of the extracted genome
#'
#' @return the phylogenetic profile of the interested genome
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

#' extract from a domains table a specific genome
#'
#' @param domains the domains table
#' @param genome name of the genome in the file, that should be remove
#' @param reverse if reverse = TRUE then the function will not extract the 
#' lines of the genome but remove it
#'
#' @return the new domains table
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

#' extract from an extended.fasta file a specific genome
#'
#' @param lines a vector that contains lines of an extended.fasta file
#' @param genome the name of the genome, that should be removed
#' @return the new vector that contains the lines of the genome
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
