#' The function read the fasta file in each core gene and extract the references
#' species based on the priority order in the priorityList
#'
#' @param fasta the path to the fasta file
#' @param priorityList a vector that contains the taxa of the core set (it does
#' not need to contain all the taxa, can be all or fewer) in a priority list
#'
#' @return The function will return the references species. If the function did
#' not find any references species, it will return the value NULL
#' @export
getSpec <- function(fasta, priorityList) {
    lines <- readLines(fasta)
    check <- 0
    for (genomeID in priorityList) {
        for (i in 1:length(lines)) {
            if (i %% 2 == 1) {
                splited <- strsplit(lines[i], "|", fixed = TRUE)[[1]]
                refGenome <- splited[2]
                if (refGenome == genomeID) {
                    check <- 1
                    break
                }
            }
        }
        if (check == 1) {
            break
        }
    }
    if (check == 0) {
        return(NULL)
    }
    return(refGenome)
}

#' The function determines the priority order of a list of taxa based on the
#' taxonomic distance of the taxa in the list with a given taxon
#'
#' @param query The genome ID of the interested genome
#' @param refSet the list of taxa in the core set
#'
#' @return The list of taxa with a priority order based on the taxanomic
#' distance
#' @export
autofindPriority <- function(query, refSet) {
    position <- function(i, j, n) {
        return(n * (i - 1) - i * (i - 1) / 2 + j - i)
    }

    swap <- function(i, j, set) {
        x <- set[i]
        y <- set[j]
        set[i] <- y
        set[j] <- x
        return(set)
    }
    queryNcbi <- strsplit(query, "@", fixed = TRUE)[[1]][2]

    refNcbiSet <- lapply(refSet, function(genomeID) {
        return(strsplit(genomeID, "@", fixed = TRUE)[[1]][2])
    })
    refNcbiSet <- unlist(refNcbiSet)

    if (queryNcbi %in% refNcbiSet) {
        check <- 1
        ncbiSet <- c()
        for (ncbi in refNcbiSet) {
            if (ncbi != queryNcbi) {
                ncbiSet <- c(ncbiSet, ncbi)
            }
        }
        refNcbiSet <- c(queryNcbi, ncbiSet)
    } else {
        check <- 0
        refNcbiSet <- c(queryNcbi, refNcbiSet)
    }
    cl <- taxize::classification(refNcbiSet, db = "ncbi")
    tree <- taxize::class2tree(cl)
    print(tree$distmat)
    n <- length(refNcbiSet)
    distanceSet <- unlist(lapply(2:n, function(i, distmat, n) {
        return(distmat[position(1, i, n)])
    },
    distmat = tree$distmat,
    n = n
    ))
    print(distanceSet)

    for (i in 1:(n - 2)) {
        min <- distanceSet[i]
        index <- i
        for (j in (i + 1):(n - 1)) {
            if (distanceSet[j] < min) {
                min <- distanceSet[j]
                index <- j
            }
        }
        distanceSet <- swap(i, index, distanceSet)
        refNcbiSet <- swap(i + 1, index + 1, refNcbiSet)
    }
    if (check == 0) {
        refNcbiSet <- refNcbiSet[2:length(refNcbiSet)]
    }

    refSet <- unlist(lapply(refNcbiSet,
        function(ncbi, refSet) {
            r <- NULL
            for (take in refSet) {
                n <- strsplit(take, "@", fixed = TRUE)[[1]][2]
                if (n == ncbi) {
                    r <- take
                    break
                }
            }
            return(r)
        },
        refSet = refSet
    ))
    print(refSet)
    return(refSet)
}
