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
#' @examples
#' #' ## Create a pseudo domains table
#' V1 <- c(
#'     "1001705at2759#1001705at2759|HUMAN@9606@3|Q15291|1",
#'     "551907at2759#551907at2759|AMPQU@400682@2|400682_0:001490|1"
#' )
#' V2 <- c(
#'     "1001705at2759|HOMSA@9606@2|9606_0:00004c",
#'     "551907at2759|AMPQU@400682@2|400682_0:001490|1"
#' )
#' V3 <- c(538, 1087)
#' V4 <- c("pfam_WD40", "flps_SINGLE_{G}")
#' V5 <- c(17, 1030)
#' V6 <- c(52, 1058)
#' V7 <- c(0.2265, 0.0936)
#' V8 <- c("Y", "Y")
#' domains <- data.frame(V1, V2, V3, V4, V5, V6, V7, V8)
#'
#' ## Write it in a file
#' wd <- getwd()
#' filePath <- paste(wd, "/fCAT_functiontest.domains", sep = "")
#' write.table(domains, filePath,
#'     sep = "\t", row.names = FALSE,
#'     quote = FALSE, col.names = FALSE
#' )
#'
#' ## correcting
#' correctDomains(filePath, "HUMAN@9606@3")
#'
#' ## Check if the file is corrected
#' read.table(filePath, sep = "\t", header = FALSE, comment.char = "")
#'
#' ## delete file
#' file.remove(filePath)
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
    if (nrow(domains) == 0) {
        file.remove(domainsFile)
    } else {
        write.table(
            domains, domainsFile,
            quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = FALSE
        )
    }
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
#' @examples
#' ## Create pseudo phylogenetic profile, which contains two different genomes
#' geneID <- c("530670", "530730", "603043")
#' ncbiID <- c("ncbi9606", "ncbi9606", "ncbi9606")
#' orthoID <- c(
#'     "530670|HUMAN@9606@3|Q16526|1", "530730|HUMAN@9606@3|P05091|1",
#'     "603043|HOMSA@9606@2|Q39233|1"
#' )
#' FAS_F <- c(1, 1, 1)
#' FAS_B <- c(1, 1, 1)
#' pp <- data.frame(geneID, ncbiID, orthoID, FAS_F, FAS_B)
#'
#' ## Write the pp file in a file
#' wd <- getwd()
#' filePath <- paste(wd, "/fCAT_functiontest.phyloprofile", sep = "")
#' write.table(pp, filePath, sep = "\t", row.names = FALSE, quote = FALSE)
#'
#' ## Correcting
#' correctPP(filePath, "HUMAN@9606@3")
#'
#' ## Check if the file is corrected
#' read.table(filePath, header = TRUE, sep = "\t")
#'
#' ## Delete the file
#' file.remove(filePath)
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
    if (nrow(PP) == 0) {
        file.remove(PPFile)
    } else {
        write.table(PP, PPFile, row.names = FALSE, quote = FALSE, sep = "\t")
    }
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
#' @examples
#' ## Create pseudo extended fasta file
#' fasta <- c(
#'     ">1001705at2759|HUMAN@9606@3|Q15291|1",
#'     "MNLELLESFGQNYPEEADGTLDCISMALTCTFNRWGT",
#'     ">1489230at2759|ARATH@3702@2|3702_0:005bc2|1",
#'     "MAGRATIPARNSALIAMIADEDTVVGFLMAGVGNVDIRRKTNYLIVDS"
#' )
#'
#' ## Write it in a file in the current working directory
#' wd <- getwd()
#' filePath <- paste(wd, "/fCAT_functiontest.fa", sep = "")
#' writeLines(fasta, filePath)
#'
#' ## Correcting
#' correctFasta(filePath, "HUMAN@9606@3")
#'
#' ## Check if the file is corrected
#' fasta <- readLines(filePath)
#' print(fasta)
#'
#' ## Delete the file
#' file.remove(filePath)
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

#' This function run correctPP, correctFasta, correctDomains, 
#' correctPriority and correctReport in a loop to correct all phyloprofile, 
#' extended fasta, domains, prioritylist and report file in a folder
#' @param directory the path to the directory, that contains the files
#' @param genome name of the genome, that should be removed
#' @return none
#' @examples
#' ## Create some pseudo data
#' genomeID <- c("HUMAN@9606@3", "AMPQU@400682@2")
#' priority_list <- c("HUMAN@9606@3", "AMPQU@400682@2")
#' table1 <- data.frame(genomeID, priority_list)
#'
#' genomeID <- c("HUMAN@9606@3", "AMPQU@400682")
#' similar <- c(330, 313)
#' dissimilar <- c(3, 0)
#' missing <- c(4, 11)
#' duplicated <- c(1, 0)
#' ignored <- c(8, 22)
#'
#' table2 <- data.frame(
#'     genomeID, similar, dissimilar, missing, duplicated,
#'     ignored
#' )
#'
#' ## Write them in files in a subfolder in the current working directory
#' wd <- getwd()
#' testFolder <- paste(wd, "/fCAT_functiontest", sep = "")
#' dir.create(testFolder, recursive = TRUE)
#' filePath1 <- paste(testFolder, "/table1.prioritylist", sep = "")
#' filePath2 <- paste(testFolder, "/table2.report", sep = "")
#' write.table(table1, filePath1, sep = "\t", row.names = FALSE, quote = FALSE)
#' write.table(table2, filePath2, sep = "\t", row.names = FALSE, quote = FALSE)
#'
#' ## Correcting
#' correctFiles(testFolder, "HUMAN@9606@3")
#'
#' ## Check if the file in the test folder are corrected
#' table1 <- read.table(filePath1, sep = "\t", header = TRUE)
#' print.data.frame(table1)
#' table2 <- read.table(filePath2, sep = "\t", header = TRUE)
#' print.data.frame(table2)
#'
#' ## Delete the folder
#' unlink(testFolder, recursive = TRUE)
#' @export
correctFiles <- function(directory, genome) {
    files <- list.files(directory, full.names = TRUE)

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
    }
}
