#' This function take the path to the interested genome and its FAS 
#' annotation (if FAS annotation was not provided, the function will call 
#' annoFAS tool to compute it) as the input. The function will make a copy of 
#' the genome fasta file and a symbolic link (not in case annoFAS is called) 
#' and arrange them into the query_taxon folder and into the weight_dir folder 
#' of the core directory.
#'
#' @param genome the path to the fasta file of the genome
#' @param fasAnno the path to the fas annotation file of the genome. If equal
#' NULL, the function will compute the annotation, and arrange it into the
#' weight_dir folder of the core directory
#' @param root The path to the core directory, where the core set is stored
#' within weight_dir, blast_dir, etc.
#' @param process A logical option to determine if the function was used as a
#' modul in the processCoreSet function
#' @param weightDir The user can replace the weight_dir folder in the core
#' directory by specifying the path to the replacing folder in this argument
#'
#' @return none
#' @examples
#' ## Take the demo data
#' coreFolder <- system.file("extdata", "sample", package = "fCAT")
#' genome <- system.file("extdata", "HUMAN@9606@3.fa", package = "fCAT")
#' fasAnno <- system.file("extdata", "HUMAN@9606@3.json", package = "fCAT")
#' ## Place seed
#' placeSeed(genome, fasAnno, coreFolder)
#' ## Test if the seed was place
#' fastaPath <- paste(coreFolder, "/query_taxon/HUMAN@9606@3/HUMAN@9606@3.fa",
#'     sep = ""
#' )
#' annoPath <- paste(coreFolder, "/weight_dir/HUMAN@9606@3.json", sep = "")
#' if (file.exists(fastaPath) && file.exists(annoPath)) {
#'     print("Seed is placed")
#' }
#'
#' ## Delete seed
#' fastaFolder <- paste(coreFolder, "/query_taxon/HUMAN@9606@3", sep = "")
#' annoPath <- paste(coreFolder, "/weight_dir/HUMAN@9606@3.json", sep = "")
#' unlink(fastaFolder, recursive = TRUE)
#' file.remove(annoPath)
#' @export
placeSeed <- function(
    genome, fasAnno = NULL, root, process = FALSE, weightDir = NULL) {
    if (!endsWith(root, "/")) {
        root <- paste(root, "/", sep = "")
    }

    splited <- strsplit(genome, "/", fixed = TRUE)[[1]]
    splited <- strsplit(splited[length(splited)], ".", fixed = TRUE)[[1]]
    genomeName <- splited[1]

    genomePath <- paste(root, "query_taxon", sep = "")
    if (!dir.exists(genomePath)) {
        dir.create(genomePath)
    }

    if (!is.null(weightDir)) {
        weightPath <- weightDir
    } else {
        weightPath <- paste(root, "weight_dir", sep = "")
    }

    if (process == FALSE) {
        if (!is.null(fasAnno)) {
            R.utils::createLink(
                paste(weightPath, "/", genomeName, ".json", sep = ""),
                fasAnno,
                overwrite = TRUE
            )
        } else {
            command <- paste(
                "annoFAS", "-i", genome, "-o", weightPath,
                "-n", genomeName,
                sep = " "
            )
            system(command)
        }
    }

    dir.create(paste(genomePath, "/", genomeName, sep = ""))
    file.copy(
        genome,
        paste(genomePath, "/", genomeName, "/", genomeName, ".fa", sep = "")
    )
}
