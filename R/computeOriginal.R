#' The function compute the initial original phylogenetic profile of all the
#' taxa in the score set
#'
#' @param root the path to the root folder
#' @param coreSet the core set name
#' @param scoreMode the mode determines the way to assess the founded ortholog
#' @param cpu determines the number of the cores
#'
#' @return none
#' @export
computeOriginal <- function(root, coreSet, scoreMode, cpu = 4) {
    if (!endsWith(root, "/")) {
        root <- paste(root, "/", sep = "")
    }

    genomeDir <- paste(root, "genome_dir", sep = "")
    weightDir <- paste(root, "weight_dir", sep = "")

    for (
        genome in list.dirs(genomeDir, full.names = FALSE, recursive = FALSE)
    ) {
        genomeFasta <- paste(
            genomeDir, "/", genome, "/", genome, ".fa", sep = ""
        )
        fasAnno <- paste(weightDir, "/", genome, ".json", sep = "")
        computeReport(genomeFasta, fasAnno, root, coreSet,
            TRUE, scoreMode,
            priorityList = c(genome), cpu, computeOri = TRUE
        )
    }
}
