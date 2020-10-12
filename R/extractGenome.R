#' This function takes a phylogenetic profile, which contains many different
#' genomes as an input and extract the phylogenetic profile of a specific genome
#'
#' @param pp A phylogenetic profile in data.frame
#' @param genome The genome ID of the genome, whose pp need to be extracted.
#' Exp: genome = HUMAN@9606@3
#'
#' @return The phylogenetic profile of the interested genome in data.frame
#' @examples 
#' ## Create pseudo phylogenetic profile, which contains two different genomes
#' geneID <- c("530670", "530730", "603043")
#' ncbiID <- c("ncbi9606", "ncbi9606", "ncbi9606")
#' orthoID <- c("530670|HUMAN@9606@3|Q16526|1", "530730|HUMAN@9606@3|P05091|1",
#' "603043|HOMSA@9606@2|Q39233|1")
#' FAS_F <- c(1, 1, 1)
#' FAS_B <- c(1, 1, 1)
#' pp <- data.frame(geneID, ncbiID, orthoID, FAS_F, FAS_B)
#' newPP <- extractPP(pp, "HUMAN@9606@3")
#' print.data.frame(newPP)
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
#' @examples
#' ## Create a pseudo domains table
#' V1 <- c("1001705at2759#1001705at2759|HUMAN@9606@3|Q15291|1",
#' "551907at2759#551907at2759|AMPQU@400682@2|400682_0:001490|1")
#' V2 <- c("1001705at2759|HOMSA@9606@2|9606_0:00004c", 
#' "551907at2759|AMPQU@400682@2|400682_0:001490|1")
#' V3 <- c(538, 1087)
#' V4 <- c("pfam_WD40", "flps_SINGLE_{G}")
#' V5 <- c(17, 1030)
#' V6 <- c(52, 1058)
#' V7 <- c(0.2265, 0.0936)
#' V8 <- c("Y", "Y")
#' domains <- data.frame(V1, V2, V3, V4, V5, V6, V7, V8)
#' ## extract domains table of HUMAN@9606@3
#' newDomains <- extractDomains(domains, "HUMAN@9606@3")
#' print.data.frame(newDomains)
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
#' @examples 
#' ## Create pseudo extended fasta file
#' fasta <- c(">1001705at2759|HUMAN@9606@3|Q15291|1",
#' "MNLELLESFGQNYPEEADGTLDCISMALTCTFNRWGT", 
#' ">1489230at2759|ARATH@3702@2|3702_0:005bc2|1",
#' "MAGRATIPARNSALIAMIADEDTVVGFLMAGVGNVDIRRKTNYLIVDS")
#' ## Extract extended fasta of HUMAN@96063
#' newFasta <- extractFasta(fasta, "HUMAN@9606@3")
#' print(newFasta)
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