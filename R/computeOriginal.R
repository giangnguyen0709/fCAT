#' The function compute the initial original phylogenetic profile of all the
#' taxa in the score set
#' 
#' @usage 
#' computeOriginal(
#'     coreDir, coreSet, scoreMode, cpu = 4, cleanup = FALSE,
#'     ppDir = NULL
#' ) 
#'
#' @param coreDir The path to the core directory, where the core set is stored
#' within weight_dir, blast_dir, etc.
#' @param coreSet The name of the interested core set. The core directory can
#' contains more than one core set and the user must specify the interested
#' core set. The core set will be stored in the folder core_orthologs in
#' subfolder, specify them by the name of the subfolder
#' @param scoreMode the mode determines the method to scoring the founded
#' ortholog and how to classify them. Choices: 1, 2, 3, "busco"
#' @param cpu determines the cores that fDOG and fdogFAS will uses to be run
#' parallel
#' @param cleanup The fDOG's output is a set of phylogenetic profile of each
#' core group to the interested genome. The phylogenetic profile will be stored
#' into a folder in the core set. The function will merge all the small
#' phylogenetic profile, calculate the FAS score or length to have the whole
#' phylogenetic profile of the interested genome to the core set. This fDOG's
#' output can be reused for all score modes. When cleanup is set to TRUE, the
#' fDOG's output will not be stored to be reused but to be removed
#' @param ppDir The user can replace the default folder output in the core
#' directory, where the phylogenetic profiles are stored by his folder. The user
#' can specify the path to his folder in this argument
#'
#' @return none
#' @examples
#' coreFolder <- system.file("extdata", "sample", package = "fCAT")
#' computeOriginal(coreFolder, "test",
#'     scoreMode = 2, cleanup = TRUE,
#'     ppDir = NULL
#' )
#' @export
computeOriginal <- function(
    coreDir, coreSet, scoreMode, cpu = 4, cleanup = FALSE,
    ppDir = NULL) {
    if (!endsWith(coreDir, "/")) {
        coreDir <- paste(coreDir, "/", sep = "")
    }
    if (!checkPreProcess(coreDir, coreSet)) {
        processCoreSet(coreDir, coreSet)
    }

    genomeDir <- paste(coreDir, "genome_dir", sep = "")
    weightDir <- paste(coreDir, "weight_dir", sep = "")

    for (
        genome in list.dirs(genomeDir, full.names = FALSE, recursive = FALSE)
    ) {
        if (!checkExist(genome, coreDir, coreSet, scoreMode, ppDir)) {
            genomeFasta <- paste(
                genomeDir, "/", genome, "/", genome, ".fa",
                sep = ""
            )
            fasAnno <- paste(weightDir, "/", genome, ".json", sep = "")
            computeReport(genomeFasta, fasAnno, coreDir, coreSet,
                extend = TRUE, scoreMode,
                priorityList = c(genome), cpu, computeOri = TRUE,
                cleanup = cleanup,
                ppDir = ppDir
            )
        }
    }
}
