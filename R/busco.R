#' In the scoreMode "busco" the tool must calculate the length of the founded
#' ortholog to perfome the assessment. After HaMStR founded an ortholog the
#' function will be called to update the length of the ortholog in phylogenetic
#' profile
#'
#' @param root the path to the root folder
#' @param coreSet the core set name
#' @param coreGene the ID of the core gene
#'
#' @return none
#' @export
updateLength <- function(root, coreSet, coreGene) {
    if (!endsWith(root, "/")) {
        root <- paste(root, "/", sep = "")
    }
    exFasta <- readLines(
        paste(
            root, "phyloprofile", "/", coreSet, "/", "busco", "/", "hamstrout", 
            "/", coreGene, "/", coreGene, ".extended.fa", sep = ""
        )
    )

    pp <- read.table(
        paste(
            root, "phyloprofile", "/", coreSet, "/", "busco", "/",
            "hamstrout", "/", coreGene, "/", coreGene, ".phyloprofile",
            sep = ""
        ),
        header = TRUE,
        sep = "\t"
    )
    i <- (1:nrow(pp))

    orthoLength <- unlist(
        lapply(
            i, 
            function(i, exFasta, pp) {
                orthoID <- pp$orthoID[i]
                for (j in 1:length(exFasta)) {
                    if (exFasta[j] == paste(">", orthoID, sep = "")) {
                        return(nchar(exFasta[j + 1]))
                    }
                }
            },
            exFasta = exFasta,
            pp = pp
    ))
    pp <- cbind(pp, length = orthoLength)
    write.table(
        pp,
        paste(
            root, "phyloprofile", "/", coreSet, "/", "busco", "/",
            "hamstrout", "/", coreGene, "/",
            coreGene, ".phyloprofile",
            sep = ""
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
}
